#include "../include/shapefile_io_operations.h"

#include <map>

using namespace boost::geometry;

//helper function to convert OGR polygons to CGAL polygons
Polygon_wh OGRPolygonToCGAL(OGRPolygon* poPolygon) {
	Polygon outer;
	// Exterior ring
	OGRLinearRing* exterior = poPolygon->getExteriorRing();
	if (exterior) {
		int n = exterior->getNumPoints();
		if (n > 0) {
			// Add first point explicitly
			outer.push_back(K::Point_2(exterior->getX(0), exterior->getY(0)));

			// Add subsequent points only if different from previous one
			for (int i = 1; i < n - 1; ++i) { // skip last point since it's a duplicate of the first
				double x = exterior->getX(i);
				double y = exterior->getY(i);
				K::Point_2 curr(x, y);
				if (outer.size() == 0 || curr != outer[outer.size() - 1]) {
					outer.push_back(curr);
				}
			}
		}
	}
	if (!outer.is_counterclockwise_oriented()) outer.reverse_orientation();

	std::vector<Polygon> holes;

	// Interior rings
	for (int i = 0; i < poPolygon->getNumInteriorRings(); ++i) {
		OGRLinearRing* ring = poPolygon->getInteriorRing(i);
		Polygon hole;
		int n = ring->getNumPoints();
		if (n > 0) {
			hole.push_back(K::Point_2(ring->getX(0), ring->getY(0)));

			for (int j = 1; j < n - 1; ++j) { // skip last point (duplicate of first)
				double x = ring->getX(j);
				double y = ring->getY(j);
				K::Point_2 curr(x, y);
				if (hole.size() == 0 || curr != hole[hole.size() - 1]) {
					hole.push_back(curr);
				}
			}
		}
		if (!hole.is_empty() && hole.is_simple()) {
			if (!hole.is_clockwise_oriented()) hole.reverse_orientation();
			holes.push_back(hole);
		}
	}

	//Error message on invalid polygons
	Polygon_wh poly(outer, holes.begin(), holes.end());
	if (!poly.is_valid()) {
		std::cerr << "invalid polygon: " << std::fixed << std::setprecision(8)<< endl;
		std::cerr << "POLYGON(("; for (const auto& p : poly.outer_boundary()) std::cerr << to_double(p.x()) << " " << to_double(p.y()) << ","; cout << "))" << endl;
		std::cerr << "holes: " << endl;
		for (const auto& hole: poly.hole_polygons) {
			std::cerr << "POLYGON(("; for (const auto& p : hole) std::cerr << to_double(p.x()) << " " << to_double(p.y()) << ","; cout << "))" << endl;
		}
		//abort execution as behavior is undefined with invalid polygons
		abort();
	}

	return poly;
}

bool ReadGeoPackage(const std::string& filename, std::vector<Polygon_wh>& polys) {
    // Register all drivers
    GDALAllRegister();

    // Open the dataset
    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(
        filename.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));

    if (poDS == nullptr) {
        std::cerr << "Failed to open GPKG: " << filename << std::endl;
        return 0;
    }

    // Loop over layers
    for (int i = 0; i < poDS->GetLayerCount(); ++i) {
        OGRLayer* poLayer = poDS->GetLayer(i);
        if (poLayer == nullptr) continue;

        poLayer->ResetReading();
        OGRFeature* poFeature;
        while ((poFeature = poLayer->GetNextFeature()) != nullptr) {
            OGRGeometry* poGeometry = poFeature->GetGeometryRef();
            if (poGeometry != nullptr &&
                (wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon ||
                 wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon)) {

                std::vector<Polygon_wh> feature_polys;

                // Convert geometry to CGAL polygons
                if (wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon) {
                    Polygon_wh poly = OGRPolygonToCGAL(
                        dynamic_cast<OGRPolygon*>(poGeometry));
                    feature_polys.push_back(poly);
                } else if (wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon) {
                    OGRMultiPolygon* poMulti = dynamic_cast<OGRMultiPolygon*>(poGeometry);
                    for (auto&& poPart : *poMulti) {
                        Polygon_wh poly = OGRPolygonToCGAL(
                            dynamic_cast<OGRPolygon*>(poPart));
                        feature_polys.push_back(poly);
                    }
                }

                // Add polygons to output vector
                polys.insert(polys.end(), feature_polys.begin(), feature_polys.end());
            }
            OGRFeature::DestroyFeature(poFeature);
        }
    }

    GDALClose(poDS);

	//assign ids to polygons
	for (int i=0; i<polys.size(); i++) {
		polys[i].global_id = i;
	}

    return 0;
}

bool ReadGeoPackageWithAttributes(
    const std::string& filename,
    std::vector<Polygon_wh>& polys,
    std::vector<std::map<std::string, std::string>>& attributes_per_feature)
{
    GDALAllRegister();

    GDALDataset *poDS = static_cast<GDALDataset*>(GDALOpenEx(
        filename.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));

    if (poDS == nullptr) {
        std::cerr << "Failed to open GPKG: " << filename << std::endl;
        return false;
    }

    // Loop over layers
    for (int i = 0; i < poDS->GetLayerCount(); ++i) {
        OGRLayer* poLayer = poDS->GetLayer(i);
        if (poLayer == nullptr) continue;

        // === get CRS and store it globally ===
        OGRSpatialReference* srs = poLayer->GetSpatialRef();
        if (srs != nullptr) {
            char* auth_name = nullptr;
            char* auth_code = nullptr;
            if (srs->AutoIdentifyEPSG() == OGRERR_NONE) {
                auth_name = const_cast<char*>(srs->GetAuthorityName(nullptr));
                auth_code = const_cast<char*>(srs->GetAuthorityCode(nullptr));
                if (auth_name && auth_code) {
                    global_crs = std::string(auth_name) + ":" + std::string(auth_code); // e.g., "EPSG:4326"
                }
            }
            if (global_crs.empty()) {
                char* wkt = nullptr;
                srs->exportToWkt(&wkt);
                if (wkt) {
                    global_crs = wkt;
                    CPLFree(wkt);
                }
            }
        }

        poLayer->ResetReading();
        OGRFeature* poFeature;
        while ((poFeature = poLayer->GetNextFeature()) != nullptr) {
            OGRGeometry* poGeometry = poFeature->GetGeometryRef();
            if (poGeometry != nullptr &&
                (wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon ||
                 wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon)) {

                std::vector<Polygon_wh> feature_polys;

                // Convert geometry to CGAL polygons
                if (wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon) {
                    Polygon_wh poly = OGRPolygonToCGAL(
                        dynamic_cast<OGRPolygon*>(poGeometry));
                    feature_polys.push_back(poly);
                } else if (wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon) {
                    OGRMultiPolygon* poMulti = dynamic_cast<OGRMultiPolygon*>(poGeometry);
                    for (auto&& poPart : *poMulti) {
                        Polygon_wh poly = OGRPolygonToCGAL(
                            dynamic_cast<OGRPolygon*>(poPart));
                        feature_polys.push_back(poly);
                    }
                }

                // === read all attributes into a map ===
                std::map<std::string, std::string> attr_map;
                int field_count = poFeature->GetFieldCount();
                for (int fi = 0; fi < field_count; ++fi) {
                    OGRFieldDefn* field_def = poFeature->GetDefnRef()->GetFieldDefn(fi);
                    std::string field_name = field_def->GetNameRef();

                    if (poFeature->IsFieldSetAndNotNull(fi)) {
                        switch (field_def->GetType()) {
                            case OFTInteger:
                                attr_map[field_name] = std::to_string(poFeature->GetFieldAsInteger(fi));
                                break;
                            case OFTInteger64:
                                attr_map[field_name] = std::to_string(poFeature->GetFieldAsInteger64(fi));
                                break;
                            case OFTReal:
                                attr_map[field_name] = std::to_string(poFeature->GetFieldAsDouble(fi));
                                break;
                            case OFTString:
                                attr_map[field_name] = poFeature->GetFieldAsString(fi);
                                break;
                            default:
                                attr_map[field_name] = poFeature->GetFieldAsString(fi); // fallback
                                break;
                        }
                    }
                }

                // === add each polygon and attribute map ===
                for (auto& poly : feature_polys) {
                    polys.push_back(poly);
                    attributes_per_feature.push_back(attr_map);
                }
            }
            OGRFeature::DestroyFeature(poFeature);
        }
    }

    GDALClose(poDS);

	//assign IDs to polys
	for (int i=0; i<polys.size(); i++) {
		polys[i].global_id = i;
	}

    return true;
}

void writeToGeoPackage(const std::vector<Polygon_wh>& polys, Solution& sol, bool set_index, const std::string& path) {
    GDALAllRegister();

    // Create data source (GeoPackage)
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GPKG");
    if (!driver) {
        throw std::runtime_error("GPKG driver not available");
    }

	std::string file = path + ".gpkg";
    GDALDataset* dataset = driver->Create(file.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!dataset) {
        throw std::runtime_error("Failed to create GPKG file");
    }

	// Create spatial reference from CRS string
	OGRSpatialReference srs;
	if (srs.SetFromUserInput(global_crs.c_str()) != OGRERR_NONE) {
		GDALClose(dataset);
		throw std::runtime_error("Invalid CRS: " + global_crs);
	}

	// Create layer with CRS
	OGRLayer* layer = dataset->CreateLayer("polygons", &srs, wkbPolygon, nullptr);
	if (!layer) {
		GDALClose(dataset);
		throw std::runtime_error("Failed to create layer");
	}

    // Create fields
    OGRFieldDefn field_match("match", OFTInteger);
    layer->CreateField(&field_match);

    OGRFieldDefn field_number("number", OFTInteger);
    layer->CreateField(&field_number);

    OGRFieldDefn field_weight("weight", OFTReal);
    field_weight.SetWidth(10);
    field_weight.SetPrecision(3);
    layer->CreateField(&field_weight);

    OGRFieldDefn field_mult_this("mult_this", OFTInteger);
    layer->CreateField(&field_mult_this);

    OGRFieldDefn field_mult_other("mult_other", OFTInteger);
    layer->CreateField(&field_mult_other);

    // Collect attribute data
    std::vector<int> group_index = sol.getMatchIndices(set_index);
    std::vector<double> match_weight = sol.getMatchWeights(set_index);
    auto matching_multiplicities = sol.getMatchMultiplicitiesPerPolygon();
    std::vector<std::pair<int,int>> multiplicities = !set_index ? matching_multiplicities.first : matching_multiplicities.second;

    int p = 0;
    for (const auto& poly : polys) {
        // Build geometry
        OGRPolygon ogr_poly;

        // Outer boundary
        OGRLinearRing outer_ring;
        for (auto vit = poly.outer_boundary().vertices_begin(); vit != poly.outer_boundary().vertices_end(); ++vit) {
            outer_ring.addPoint(to_double(vit->x()), to_double(vit->y()));
        }
        outer_ring.closeRings();
        ogr_poly.addRing(&outer_ring);

        // Holes
        for (auto hole: poly.hole_polygons) {
            OGRLinearRing hole_ring;
            // Need to reverse order for clockwise orientation
            std::vector<Point> vertices;
            for (auto v : hole.vertices()) {
                vertices.push_back(v);
            }
            for (auto it = vertices.rbegin(); it != vertices.rend(); ++it) {
                hole_ring.addPoint(to_double(it->x()), to_double(it->y()));
            }
            hole_ring.closeRings();
            ogr_poly.addRing(&hole_ring);
        }

        // Create feature
        OGRFeature* feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
        feature->SetGeometry(&ogr_poly);
        feature->SetField("match", group_index[p]);
        feature->SetField("number", p);
        feature->SetField("weight", match_weight[p]);
        feature->SetField("mult_this", multiplicities[p].first);
        feature->SetField("mult_other", multiplicities[p].second);

        if (layer->CreateFeature(feature) != OGRERR_NONE) {
            OGRFeature::DestroyFeature(feature);
            GDALClose(dataset);
            throw std::runtime_error("Failed to create feature");
        }

        OGRFeature::DestroyFeature(feature);
        p++;
    }

    GDALClose(dataset);
}

void writeToGeoPackageWithAttributes(
    const std::vector<Polygon_wh>& polys,
    Solution& sol,
    bool set_index,
    const std::vector<std::map<std::string, std::string>>& attributes,
    const std::string& path)
{
    if (polys.size() != attributes.size()) {
        throw std::runtime_error("Vector sizes of polys and attributes must match!");
    }

    GDALAllRegister();

    // Create data source (GeoPackage)
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GPKG");
    if (!driver) {
        throw std::runtime_error("GPKG driver not available");
    }

    std::string file = path + ".gpkg";
    GDALDataset* dataset = driver->Create(file.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!dataset) {
        throw std::runtime_error("Failed to create GPKG file");
    }

    // Create spatial reference from CRS string
    OGRSpatialReference srs;
    if (srs.SetFromUserInput(global_crs.c_str()) != OGRERR_NONE) {
        GDALClose(dataset);
        throw std::runtime_error("Invalid CRS: " + global_crs);
    }

    // Create layer
    OGRLayer* layer = dataset->CreateLayer("polygons", &srs, wkbPolygon, nullptr);
    if (!layer) {
        GDALClose(dataset);
        throw std::runtime_error("Failed to create layer");
    }

    // === Dynamically add fields from attributes ===
    if (!attributes.empty()) {
        const auto& first_attr = attributes[0];
        for (const auto& kv : first_attr) {
            OGRFieldDefn fld(kv.first.c_str(), OFTString); // store as string, safe choice
            fld.SetWidth(254);
            if (layer->CreateField(&fld) != OGRERR_NONE) {
                GDALClose(dataset);
                throw std::runtime_error("Failed to create attribute field: " + kv.first);
            }
        }
    }

    // === Add your existing fields ===
    OGRFieldDefn field_match("match", OFTInteger);
    layer->CreateField(&field_match);

    OGRFieldDefn field_number("number", OFTInteger);
    layer->CreateField(&field_number);

    OGRFieldDefn field_weight("weight", OFTReal);
    field_weight.SetWidth(10);
    field_weight.SetPrecision(3);
    layer->CreateField(&field_weight);

    OGRFieldDefn field_mult_this("mult_this", OFTInteger);
    layer->CreateField(&field_mult_this);

    OGRFieldDefn field_mult_other("mult_other", OFTInteger);
    layer->CreateField(&field_mult_other);

    // === Collect your attribute data from Solution ===
    std::vector<int> group_index = sol.getMatchIndices(set_index);
    std::vector<double> match_weight = sol.getMatchWeights(set_index);
    auto matching_multiplicities = sol.getMatchMultiplicitiesPerPolygon();
    std::vector<std::pair<int,int>> multiplicities = !set_index ? matching_multiplicities.first : matching_multiplicities.second;

    if (polys.size() != group_index.size() || polys.size() != match_weight.size() || polys.size() != multiplicities.size()) {
        GDALClose(dataset);
        throw std::runtime_error("Vector sizes from Solution do not match polys!");
    }

    // === Write features ===
    for (size_t p = 0; p < polys.size(); ++p) {
        const auto& poly = polys[p];
        const auto& attr = attributes[p];

        // Build geometry
        OGRPolygon ogr_poly;

        // Outer boundary
        OGRLinearRing outer_ring;
        for (auto vit = poly.outer_boundary().vertices_begin(); vit != poly.outer_boundary().vertices_end(); ++vit) {
            outer_ring.addPoint(to_double(vit->x()), to_double(vit->y()));
        }
        outer_ring.closeRings();
        ogr_poly.addRing(&outer_ring);

        // Holes
        for (auto hole: poly.hole_polygons) {
            OGRLinearRing hole_ring;
            // Reverse to get clockwise orientation
            std::vector<Point> vertices;
            for (auto v : hole.vertices()) vertices.push_back(v);
            for (auto it = vertices.rbegin(); it != vertices.rend(); ++it) {
                hole_ring.addPoint(to_double(it->x()), to_double(it->y()));
            }
            hole_ring.closeRings();
            ogr_poly.addRing(&hole_ring);
        }

        // Create feature
        OGRFeature* feature = OGRFeature::CreateFeature(layer->GetLayerDefn());
        feature->SetGeometry(&ogr_poly);

        // Set dynamic attributes
        for (const auto& kv : attr) {
            feature->SetField(kv.first.c_str(), kv.second.c_str());
        }

        // Set your existing attributes
        feature->SetField("match", group_index[p]);
        feature->SetField("number", static_cast<int>(p));
        feature->SetField("weight", match_weight[p]);
        feature->SetField("mult_this", multiplicities[p].first);
        feature->SetField("mult_other", multiplicities[p].second);

        if (layer->CreateFeature(feature) != OGRERR_NONE) {
            OGRFeature::DestroyFeature(feature);
            GDALClose(dataset);
            throw std::runtime_error("Failed to create feature");
        }

        OGRFeature::DestroyFeature(feature);
    }

    GDALClose(dataset);
}

void updateGeoPackageWithAttributes(
    const std::string& path,
    const std::vector<Polygon_wh>& polys,
    Solution& sol,
    bool set_index)
{
    GDALAllRegister();

    GDALDataset* dataset = static_cast<GDALDataset*>(
        GDALOpenEx(path.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE, nullptr, nullptr, nullptr));
    if (!dataset) {
        throw std::runtime_error("Failed to open GeoPackage for update: " + path);
    }

    std::vector<int> group_index = sol.getMatchIndices(set_index);
    std::vector<double> match_weight = sol.getMatchWeights(set_index);
    auto matching_multiplicities = sol.getMatchMultiplicitiesPerPolygon();
    std::vector<std::pair<int,int>> multiplicities = !set_index ?
        matching_multiplicities.first : matching_multiplicities.second;

    // Loop over layers
    for (int i = 0; i < dataset->GetLayerCount(); ++i) {
        OGRLayer* layer = dataset->GetLayer(i);
        if (!layer) continue;

        // BEGIN TRANSACTION: speeds up massively
    	if (layer->StartTransaction() != OGRERR_NONE) {
    		GDALClose(dataset);
    		throw std::runtime_error("Failed to start transaction on layer.");
    	}

        // Add fields if missing
        if (layer->FindFieldIndex("match", TRUE) == -1) {
            layer->CreateField(new OGRFieldDefn("match", OFTInteger));
        }
        if (layer->FindFieldIndex("number", TRUE) == -1) {
            layer->CreateField(new OGRFieldDefn("number", OFTInteger));
        }
        if (layer->FindFieldIndex("weight", TRUE) == -1) {
            OGRFieldDefn fld("weight", OFTReal);
            fld.SetWidth(10);
            fld.SetPrecision(3);
            layer->CreateField(&fld);
        }
        if (layer->FindFieldIndex("mult_this", TRUE) == -1) {
            layer->CreateField(new OGRFieldDefn("mult_this", OFTInteger));
        }
        if (layer->FindFieldIndex("mult_other", TRUE) == -1) {
            layer->CreateField(new OGRFieldDefn("mult_other", OFTInteger));
        }

        layer->ResetReading();

        OGRFeature* feature;
        int p = 0;
        while ((feature = layer->GetNextFeature()) != nullptr) {
            if (p >= static_cast<int>(polys.size())) {
                OGRFeature::DestroyFeature(feature);
                break;
            }

            feature->SetField("match", group_index[p]);
            feature->SetField("number", p);
            feature->SetField("weight", match_weight[p]);
            feature->SetField("mult_this", multiplicities[p].first);
            feature->SetField("mult_other", multiplicities[p].second);

            if (layer->SetFeature(feature) != OGRERR_NONE) {
                OGRFeature::DestroyFeature(feature);
                layer->RollbackTransaction();
                GDALClose(dataset);
                throw std::runtime_error("Failed to update feature in layer.");
            }

            OGRFeature::DestroyFeature(feature);
            ++p;
        }

        // COMMIT TRANSACTION
    	if (layer->CommitTransaction() != OGRERR_NONE) {
    		GDALClose(dataset);
    		throw std::runtime_error("Failed to commit transaction on layer.");
    	}
    }

    GDALClose(dataset);
}


void writeToWKT(std::vector<Polygon_wh> polys, std::string path) {

	std::ofstream poly_file;
	poly_file.open(path);

	//insert polygons
	poly_file << "MULTIPOLYGON(";
	for (const auto& poly : polys) {
		poly_file << "((";
		for (const Point v : poly.outer_boundary().vertices()) {
			poly_file << std::fixed << v.x() << " " << v.y() << ",";

		}
		poly_file << ")),\n";
	}
	poly_file << ")" << endl;
	poly_file.close();

}

//writes Segments into File
void writeToWKT(std::vector<Segment> segments, std::string path) {

	std::ofstream poly_file;
	poly_file.open(path);

	//insert polygons
	poly_file << "MULTILINESTRING(";
	for (const auto& s : segments) {
		poly_file << "((";
		for (const Point v :  {s.source(),s.target()}) {
			poly_file << std::fixed << v.x() << " " << v.y() << ",";

		}
		poly_file << ")),\n";
	}
	poly_file << ")" << endl;
	poly_file.close();

}

//writes analysis data to csv file
void writeToCSV(Solution& s, std::string path) {
	std::ofstream file;
	file.open(path + ".csv");

	//write header 
	file << "map, poly_id, match_id, match_weight" << endl;

	//insert data per polygon
	for (bool map : {0, 1}) {
		std::vector<double> match_weights = s.getMatchWeights(map);
		std::vector<int> match_indices = s.getMatchIndices(map);
		for (int id = 0; id < match_weights.size(); id++) {
			file << map << "," << id << "," << match_indices[id] << "," << match_weights[id] << endl;
		}
	}
	
	file.close();

	//additionally write a file containing the included polygons PER MATCH for easier analysis
	file.open(path + "_per_match.csv");

	//write header
	file << "match_id,polys1,polys2,match_weight" << endl;

	int offset = abs(s.match_ids.first) - 1;

	std::vector<double> match_weights1 = s.getMatchWeights(0);
	std::vector<int> match_indices1 = s.getMatchIndices(0), match_indices2 = s.getMatchIndices(1);

	std::vector<std::stringstream> sstreams1(s.match_count), sstreams2(s.match_count);
	std::vector<double> weights(s.match_count,0.0);

	for (int i=0; i<match_indices1.size(); i++) {
		sstreams1[match_indices1[i] + offset] << i << " ";
		weights[match_indices1[i] + offset] = match_weights1[i];
	}
	for (int j=0; j<match_indices2.size(); j++) {
		sstreams2[match_indices2[j] + offset] << j << " ";
	}

	for (int id = s.match_ids.first + 1; id < s.match_ids.second; id++) {
		file << id << "," << sstreams1[id+offset].str() << "," << sstreams2[id+offset].str() << "," << weights[id+offset] << endl;
	}

	file.close();

}

//writes analysis data to csv file
void writeToCSV(std::vector<double>& obj_of_set, std::string path) {
	std::ofstream file;
	file.open(path + ".csv");

	//write header
	file << "set_id,objective" << endl;

	//insert data per polygon
	for (int i=0; i<obj_of_set.size(); i++) {
		file << i << "," << obj_of_set[i] << endl;
	}

	file.close();

}

void writeToCSV(std::vector<int> set_sizes, std::vector<std::pair<int,int>> set_sizes_after_precomp, std::vector<double> execution_times, std::vector<bool> completed_exploration, std::string path) {
	std::ofstream file;
	file.open(path + ".csv");

	//write header 
	file << "set_id,set_size,set_size after pregrouping,set_size after simple matches,execution_time,completed_exploration" << endl;

	//insert data 
	for (int id = 0; id < set_sizes.size(); id++) {
		file << id << "," << set_sizes[id] << "," << set_sizes_after_precomp[id].first << "," << set_sizes_after_precomp[id].second << "," << std::fixed << std::setprecision(3) << execution_times[id] << "," << completed_exploration[id] << endl;
	}

	file.close();

}

//reads out two csv files and puts out indices of matches in fileA, which are different to the matches of fileB
void findDifferingMatches(std::string fileA, std::string fileB) {
	std::ifstream csvA, csvB;
	csvA.open(fileA); csvB.open(fileB);
	//read both files and store as vectors
	std::vector<int> matchesA, matchesB;
    std::vector<double> weightsA, weightsB;

	//skip headers
	std::string line; 
	getline(csvA, line);

	while (getline(csvA, line)) {
		//parse line, we want the third line, as this stores the match id
		std::vector<std::string> columns;
		std::stringstream ss(line);
		while(ss.good()) {
			std::string substr;
			getline(ss, substr, ',');
			columns.push_back(substr);
		}


		matchesA.push_back(stoi(columns[2]));
        weightsA.push_back(stof(columns[3]));
	}

	getline(csvB, line);

	while (getline(csvB, line)) {
		//parse line, we want the third line, as this stores the match id
		std::vector<std::string> columns;
		std::stringstream ss(line);
		while (ss.good()) {
			std::string substr;
			getline(ss, substr, ',');
			columns.push_back(substr);
		}


		matchesB.push_back(stoi(columns[2]));
        weightsB.push_back(stof(columns[3]));
	}

	//remember if any differences have been found
	bool found_difference = false;

	int max_match_id_A = *std::max_element(matchesA.begin(), matchesA.end());
	std::vector<bool> visited_matches_A(max_match_id_A,false);

	//now the matches per polygon are stores in the matches-vectors
	//iterate over every polygon in fileA 
	for (int i = 0; i < matchesA.size(); i++) {
		int match_id_A = matchesA[i];
        double weight_match_A = weightsA[i];

		//check if match was already visited
		if (visited_matches_A[match_id_A]) continue;

		//has not been visited yet, add to visited memory
		visited_matches_A[match_id_A] = true;

		std::vector<int> polys_in_match_A = { i };

		//find every polygon with the same match ID
		for (int j = i+1; j < matchesA.size(); j++) {
			if (matchesA[j] == match_id_A) polys_in_match_A.push_back(j);
		}

		//check the same for the opposing map
		int match_id_B = matchesB[i];
        double weight_match_B = weightsB[i];
		std::vector<int> polys_in_match_B = { i };
		for (int j = 0; j < matchesB.size(); j++) {
			if (j == i) continue;
			if (matchesB[j] == match_id_B) polys_in_match_B.push_back(j);
		}

		// make sure both vectors are sorted
		std::sort(polys_in_match_A.begin(), polys_in_match_A.end());
		std::sort(polys_in_match_B.begin(), polys_in_match_B.end());

		//if the vectors are not equal, a differing match has been found, output this match
		if (polys_in_match_A != polys_in_match_B) {
			found_difference = true;
			cout << "Found differing match! ID in first matching: " << match_id_A << " | ID in second matching: " << match_id_B << endl;
			cout << "Matching A: "; for (const auto& p : polys_in_match_A) cout << p << ", "; cout << " | Obj.: " << weight_match_A << endl;
			cout << "Matching B: "; for (const auto& p : polys_in_match_B) cout << p << ", "; cout << " | Obj.: " << weight_match_B << endl;

		}
        //also check for plausability of the solutions, if match is equal, the quality should also be equal
        else if(weight_match_A != weight_match_B) {
            cout << "Found inconsistency in match " << match_id_A << " | " << match_id_B << ", even though matches are identical, their qualities differ: " << weight_match_A << " | " << weight_match_B << endl;
            //weight_match_A << " | " << weight_match_B << endl;
        }
	}	

	if (!found_difference) cout << "No differences have been found, the two matchings are equal!" << endl;
}