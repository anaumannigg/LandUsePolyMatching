#ifndef _shapefile_io_operations_included_
#define _shapefile_io_operations_included_
#include "cgal_includes.h"
#include "polygon_wh.h"
#include "solution.h"
#include "polygon_operations.h"
#include "ogrsf_frmts.h"
#include "gdal_priv.h"

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <filesystem>
namespace fs = std::filesystem;

bool ReadGeoPackage(const std::string& filename, std::vector<Polygon_wh>& polys);

bool ReadGeoPackageWithAttributes(
    const std::string& filename,
    std::vector<Polygon_wh>& polys,
    std::vector<std::map<std::string, std::string>>& attributes_per_feature);

//does the same as writeToShapeFile, but writes to a gpgk file instead, used for larger file sizes
void writeToGeoPackage(const std::vector<Polygon_wh>& polys, Solution& sol, bool set_index, const std::string& path);

void writeToGeoPackageWithAttributes(
    const std::vector<Polygon_wh>& polys,
    Solution& sol,
    bool set_index,
    const std::vector<std::map<std::string, std::string>>& attributes,
    const std::string& path);

void updateGeoPackageWithAttributes(
    const std::string& path,                    // path to copied GPKG (e.g., export/input.gpkg)
    const std::vector<Polygon_wh>& polys,       // used to know feature count/order
    Solution& sol,                              // your data source
    bool set_index);

//writes polygons into wkt-file (.txt)
void writeToWKT(std::vector<Polygon_wh> polys, std::string path);

//writes Segmets into wkt-file (.txt)
void writeToWKT(std::vector<Segment> segments, std::string path);

//writes solution data into csv
void writeToCSV(Solution& s, std::string path);

//write singular objectives of each connected component to csv
void writeToCSV(std::vector<double>& obj_of_set, std::string path);

//writes analysis data into csv
void writeToCSV(std::vector<int> set_sizes, std::vector<std::pair<int, int>> set_sizes_after_precomp, std::vector<double> execution_times, std::vector<bool> completed_exploration, std::string path);

//helper function to find the matches, in which two solutions differ
void findDifferingMatches(std::string fileA, std::string fileB);


#endif