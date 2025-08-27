import os
import re
import geopandas as gpd
import pandas as pd
import time
import shapely
from shapely.geometry import Polygon, MultiPolygon, Point
import numpy as np
from scipy.ndimage import gaussian_filter
from pathlib import Path
from dataclasses import dataclass
import yaml
from collections import Counter
import warnings
from src.coloring import apply_visualization
import src.loader as loader

#class to define granularity configuration
@dataclass
class GranularityConfig:
    cell_size : int
    raster_resolution : int

class Compute :

    def __init__(self,dir_path, subdir_path):
        """Init the Computer Object."""
        self.__current_dir = dir_path
        print(f"Loading files from {subdir_path}...")

        # Dataset name = last directory in dir_path
        self.__dataset_name = os.path.basename(os.path.normpath(dir_path))

        #find matching suffixes
        base = os.path.basename(dir_path)
        # List all files ending with _matched.gpkg in dir_path
        files = [f for f in os.listdir(subdir_path) if f.endswith("_matched.gpkg")]

        # Pattern to extract the 'X' part: e.g. 'basename_X_matched.gpkg'
        pattern = re.compile(rf"^{re.escape(base)}_(.+)_matched\.gpkg$")

        self.__geopkg_suffixes = []  # will store X, Y

        for f in files:
            match = pattern.match(f)
            if match:
                suffix = match.group(1)
                self.__geopkg_suffixes.append(suffix)

        if len(self.__geopkg_suffixes) != 2:
            raise ValueError(f"Expected 2 geopackages with suffixes, found {len(self.__geopkg_suffixes)}")

        # Reorder so that 'iacs' is always first if it exists
        if 'iacs' in self.__geopkg_suffixes and self.__geopkg_suffixes[0] != 'iacs':
            self.__geopkg_suffixes.sort(key=lambda x: 0 if x == 'iacs' else 1)

        # assign to named instance vars
        self.__suffix1, self.__suffix2 = self.__geopkg_suffixes

        self.__data1_path = os.path.join(subdir_path, f"{os.path.basename(dir_path)}_" + self.__suffix1 + "_matched.gpkg")
        self.__data2_path = os.path.join(subdir_path, f"{os.path.basename(dir_path)}_" + self.__suffix2 + "_matched.gpkg")
        self.__csv_path = os.path.join(subdir_path, f"{os.path.basename(dir_path)}_data_per_match.csv")

        # initialize granularity levels
        self.__granularity_levels = ["fine", "medium", "coarse"]

        # note : finest_resolution / ratio_for_coarser_resolutions^(len(granularity_levels)) should be integer!
        self.__finest_cell_size = 2500
        self.__finest_resolution = 800
        self.__ratio_for_coarser_resolutions = 2

        # Create aggregated file paths
        data1_path = Path(self.__data1_path)
        data2_path = Path(self.__data2_path)

        # Create dictionaries for aggregated file paths
        self.__aggregated_path1 = data1_path.with_name(f"{data1_path.stem}_aggregated{data1_path.suffix}")
        self.__aggregated_path2 = data2_path.with_name(f"{data2_path.stem}_aggregated{data2_path.suffix}")

        self.__aggregated_gdf1 = None
        self.__aggregated_gdf2 = None

        self.__gdf1 : gpd.GeoDataFrame = None
        self.__gdf2 : gpd.GeoDataFrame = None
        self.__df_matching : pd.DataFrame = None

        #filter mean of empty slice warnings for coarser heatmap computation, as this is intended
        warnings.filterwarnings("ignore", category=RuntimeWarning, message="Mean of empty slice")

    def load_files(self):
        """Try to load the already aggregated files if possible, else load the two GeoPackages and CSV from the selected subdir."""
        self.__aggregated_data_is_set = False  # Flag to know if we need to load raw gdfs

        if self.__aggregated_path1.exists() and self.__aggregated_path2.exists():
            self.__aggregated_data_is_set = True
            self.__aggregated_gdf1 = gpd.read_file(self.__aggregated_path1)
            self.__aggregated_gdf2 = gpd.read_file(self.__aggregated_path2)

        if not self.__aggregated_data_is_set:
            print("Some aggregated files missing, loading raw dataframes...")
            self.__gdf1 = gpd.read_file(self.__data1_path)
            self.__gdf2 = gpd.read_file(self.__data2_path)

        self.__df_matching = pd.read_csv(self.__csv_path, dtype={
            "match_id": int,
            "polys1": str,
            "polys2": str,
            "match_weight": float
        })

        self.__modified_original_gdfs = False

        print("Files loaded successfully.")

    def mark_unmatched_polygons(self):
        if(self.__aggregated_data_is_set):
            return
        if("is_matched" in self.__gdf1.columns and "is_matched" in self.__gdf2.columns):
            return

        self.__gdf1['is_matched'] = (self.__gdf1["match"] >= 0).astype(float)
        self.__gdf2['is_matched'] = (self.__gdf2["match"] >= 0).astype(float)
        self.__modified_original_gdfs = True

    def compute_moeller_score_of_polygons(self,poly1, poly2):
        """
        Compute match quality score based on:
        - Overlap ratios o_1 and o_2
        - Centroid proximity scores p_1 and p_2
        Returns: geometric mean of [o_1, o_2, p_1, p_2]
        """
        # Defensive checks
        if poly1 is None or poly2 is None:
            return 0.0
        if poly1.is_empty or poly2.is_empty:
            return 0.0

        # Step 1: intersection
        intersection = poly1.intersection(poly2)
        if intersection.is_empty:
            return 0.0

        # Step 2: overlap ratios
        area_intersection = intersection.area
        area_poly1 = poly1.area
        area_poly2 = poly2.area

        o_1 = area_intersection / area_poly1 if area_poly1 > 0 else 0
        o_2 = area_intersection / area_poly2 if area_poly2 > 0 else 0

        # Step 3: centroids
        centroid_intersection = intersection.centroid
        centroid_poly1 = poly1.centroid
        centroid_poly2 = poly2.centroid

        # Step 4: compute differences
        diff1 = poly1.difference(poly2)
        diff2 = poly2.difference(poly1)

        # Step 5: get furthest centroid distance for each difference
        def get_furthest_distance(diff_geom, reference_centroid):
            if diff_geom.is_empty:
                return 0
            # Ensure always iterable
            if isinstance(diff_geom, (Polygon)):
                diff_parts = [diff_geom]
            elif isinstance(diff_geom, MultiPolygon):
                diff_parts = list(diff_geom.geoms)
            else:
                # unexpected geometry
                return 0

            # Compute max distance
            return max(reference_centroid.distance(part.centroid) for part in diff_parts)

        furthest_d1 = get_furthest_distance(diff1, centroid_intersection)
        furthest_d2 = get_furthest_distance(diff2, centroid_intersection)

        # Distances
        d1 = centroid_intersection.distance(centroid_poly1)
        d2 = centroid_intersection.distance(centroid_poly2)

        # Step 6: compute p_i
        p_1 = 1 - (d1 / furthest_d1) if furthest_d1 > 0 else 0
        p_2 = 1 - (d2 / furthest_d2) if furthest_d2 > 0 else 0

        # Avoid negatives in geometric mean
        o_1 = max(o_1, 0)
        o_2 = max(o_2, 0)
        p_1 = max(p_1, 0)
        p_2 = max(p_2, 0)

        geo_mean = (o_1 * o_2 * p_1 * p_2) ** 0.25

        return geo_mean

    def compute_moeller_score(self):
        if(self.__aggregated_data_is_set):
            return

        """Compute new score using self.__gdf_iacs, self.__gdf_sat, and df_match."""
        print("Computing new score...")

        score_string = "moeller_score"

        #first check if score is already computed, skip in this case
        if score_string in self.__gdf1.columns and score_string in self.__gdf2.columns:
            return

        #score is not computed yet, compute it
        # Build lookup tables for fast access
        iacs_geom_map = self.__gdf1.set_index('number').geometry.to_dict()
        sat_geom_map = self.__gdf2.set_index('number').geometry.to_dict()
        new_scores = []

        for idx, row in self.__df_matching.iterrows():

            ids_iacs = str(row["polys1"]).split() if pd.notna(row["polys1"]) else []
            ids_sat  = str(row["polys2"]).split() if pd.notna(row["polys2"]) else []

            if ids_iacs==[] or ids_sat ==[]:
                new_scores.append(0)
                continue

            # Collect geometries
            iacs_geoms = [iacs_geom_map.get(int(id_)) for id_ in ids_iacs if int(id_) in iacs_geom_map]
            sat_geoms  = [sat_geom_map.get(int(id_)) for id_ in ids_sat  if int(id_) in sat_geom_map]

            # Filter out Nones (missing geometries)
            iacs_geoms = [geom for geom in iacs_geoms if geom is not None]
            sat_geoms  = [geom for geom in sat_geoms if geom is not None]

            # Dissolve to single MultiPolygon (or Polygon)
            poly1 = shapely.union_all(iacs_geoms) if iacs_geoms else None
            poly2 = shapely.union_all(sat_geoms)  if sat_geoms  else None


            # Compute custom score
            score = self.compute_moeller_score_of_polygons(poly1, poly2)

            new_scores.append(score)

        # === Add / overwrite column in dataframe ===
        self.__df_matching['moeller_score'] = new_scores

        # Update self.__gdf_sat
        self.__gdf2 = self.__gdf2.drop(columns=["moeller_score"], errors="ignore")  # remove old score if exists
        self.__gdf2 = self.__gdf2.merge(
            self.__df_matching[["match_id", "moeller_score"]],
            left_on="match",
            right_on="match_id",
            how="left"
        ).drop(columns=["match_id"])

        # Update self.__gdf_iacs
        self.__gdf1 = self.__gdf1.drop(columns=["moeller_score"], errors="ignore")
        self.__gdf1 = self.__gdf1.merge(
            self.__df_matching[["match_id", "moeller_score"]],
            left_on="match",
            right_on="match_id",
            how="left"
        ).drop(columns=["match_id"])

        self.__modified_original_gdfs = True

    def compute_consistencies(self, config_path="config.yaml"):
        """
        Marks polygons in both gdf1 and gdf2 with 'is-matched-correctly'
        based on majority category decision per match, using YAML config.
        Only runs if dataset name exists in config.
        """
        if self.__aggregated_data_is_set:
            return
        if "is_matched_correctly" in self.__gdf1.columns and "is_matched_correctly" in self.__gdf2.columns:
            return

        # Load config
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        dataset_name = self.__dataset_name  # Assumes you store the current dataset name here
        if dataset_name not in config:
            print(f"Dataset '{dataset_name}' not in config — skipping consistency computation.")
            return

        ds_cfg = config[dataset_name]

        # Helper to map categories for each gdf
        def map_categories(gdf, gdf_key):
            col = ds_cfg[gdf_key]["category_column"]
            mapping = ds_cfg[gdf_key]["mapping"]
            return gdf[col].map(mapping)

        self.__gdf1["_category"] = map_categories(self.__gdf1, "gdf1")
        self.__gdf2["_category"] = map_categories(self.__gdf2, "gdf2")

        # Init correctness column
        self.__gdf1["is_matched_correctly"] = False
        self.__gdf2["is_matched_correctly"] = False

        # Lookup for polygon numbers
        gdf1_lookup = {num: idx for idx, num in self.__gdf1["number"].items()}
        gdf2_lookup = {num: idx for idx, num in self.__gdf2["number"].items()}

        # Iterate matches
        for _, row in self.__df_matching.iterrows():
            match_id = row["match_id"]
            if match_id <= 0:
                continue

            polys1_nums = [int(x) for x in row["polys1"].split() if x.strip().isdigit()]
            polys2_nums = [int(x) for x in row["polys2"].split() if x.strip().isdigit()]

            # Gather categories
            cats = []
            for num in polys1_nums:
                idx = gdf1_lookup.get(num)
                if idx is not None:
                    cats.append(self.__gdf1.at[idx, "_category"])
            for num in polys2_nums:
                idx = gdf2_lookup.get(num)
                if idx is not None:
                    cats.append(self.__gdf2.at[idx, "_category"])

            if not cats:
                continue

            # Majority vote
            majority_cat, _ = Counter(cats).most_common(1)[0]

            #in a 1:1, if both polys do not agree, we do not know the correct category
            if len(polys1_nums) == 1 and len(polys2_nums) == 1 and len(cats) == 2 and cats[0] == cats[1]:
                majority_cat = None

            # Mark correctness
            for num in polys1_nums:
                idx = gdf1_lookup.get(num)
                if idx is not None and self.__gdf1.at[idx, "_category"] == majority_cat:
                    self.__gdf1.at[idx, "is_matched_correctly"] = True
            for num in polys2_nums:
                idx = gdf2_lookup.get(num)
                if idx is not None and self.__gdf2.at[idx, "_category"] == majority_cat:
                    self.__gdf2.at[idx, "is_matched_correctly"] = True

        # Drop temp
        self.__gdf1.drop(columns="_category", inplace=True)
        self.__gdf2.drop(columns="_category", inplace=True)

        self.__modified_original_gdfs = True

    def overwrite_geopackages(self):
        """Overwrite the existing GeoPackages with updated data."""
        print("Overwriting GeoPackages...")

        #skip if files have not been modified
        if(not self.__modified_original_gdfs):
            return

        if os.path.exists(self.__data1_path):
            os.remove(self.__data1_path)

        if os.path.exists(self.__data2_path):
            os.remove(self.__data2_path)

        self.__gdf1.to_file(self.__data1_path, driver="GPKG")
        self.__gdf2.to_file(self.__data2_path, driver="GPKG")

    def write_heatmap_to_file(self, arr, x_limits,y_limits, filename):
        with open(filename, 'w') as f:
            # Write limits in the first row
            f.write(f"{x_limits[0]},{x_limits[1]},{y_limits[0]},{y_limits[1]}\n")
            # Write each row of arr
            for row in arr:
                f.write(','.join(map(str, row)) + '\n')

    def generate_heatmaps(self):
        """Generate and save heatmap based on current data."""
        print("Generating heatmaps...")

        #map modes to output titles
        mode_titles = {
            "quality" : "quality-IoU",
            "mult" : "multiplicity",
            "moeller_score" : "quality-combined",
            "unmatched" : "is-matched",
            "is_matched_correctly" : "is-matched-correctly"
        }

        #make sure output directory exists
        if not os.path.exists("heatmaps"):
            os.makedirs("heatmaps")


        # load country polygon
        country_name = loader.detect_country_from_filename(self.__current_dir, "region_boundaries")
        country_gdf : gpd.GeoDataFrame = None
        if country_name:
            path = f"region_boundaries/{country_name}.gpkg"
            country_gdf = gpd.read_file(path).to_crs("EPSG:3857")
        country_polygon = country_gdf.geometry.values[0] if country_gdf is not None else None



        # Check if aggregated files already exist
        if not self.__aggregated_data_is_set:
            # make sure all geodataframes are in the same CRS
            self.__gdf2.to_crs("EPSG:3857", inplace=True)
            self.__gdf1.to_crs("EPSG:3857", inplace=True)

            # generate grid of country
            grid = loader.generate_grid_over_polygon(country_gdf, self.__finest_cell_size)
            self.__aggregated_gdf1 = loader.aggregate_attributes_to_grid(self.__gdf1, grid)
            self.__aggregated_gdf2 = loader.aggregate_attributes_to_grid(self.__gdf2, grid)

            # save aggregated results
            self.__aggregated_gdf1.to_file(self.__aggregated_path1, driver="GPKG")
            self.__aggregated_gdf2.to_file(self.__aggregated_path2, driver="GPKG")

        # extract bounds
        min_x1, min_y1, max_x1, max_y1 = self.__aggregated_gdf1.total_bounds
        min_x2, min_y2, max_x2, max_y2 = self.__aggregated_gdf2.total_bounds
        min_x = min(min_x1, min_x2)
        min_y = min(min_y1, min_y2)
        max_x = max(max_x1, max_x2)
        max_y = max(max_y1, max_y2)


        for perspective in [self.__suffix1, self.__suffix2]:
            for mode in ["quality","mult","moeller_score","unmatched","is_matched_correctly"]:
                #skip the is_matched_correctly_mode, it may not be set for all datasets
                if(mode=="is_matched_correctly" and self.__aggregated_gdf1.get(mode) is None):
                    continue

                #load styled gdf
                styled_gdf = apply_visualization(self.__aggregated_gdf1, mode) if perspective == self.__suffix1 else apply_visualization(self.__aggregated_gdf2, mode)

                # 1. Define raster resolution — choose something reasonable for performance
                raster_size = self.__finest_resolution

                # 2. Prepare a grid for rasterization
                x_grid = np.linspace(min_x, max_x, raster_size)
                y_grid = np.linspace(min_y, max_y, raster_size)

                heatmap = np.full((raster_size, raster_size), np.nan, dtype=float)

                # 3. For each polygon, fill the raster cells that intersect it with its value
                # This is a simple rasterization: assign polygon value to raster pixels it covers

                # Convert polygons to shapely geometries
                polygons = styled_gdf.geometry.values
                values = styled_gdf['color_value'].values  # We'll need numeric values for filtering

                # Build 2D grid of cell centers
                xv, yv = np.meshgrid((x_grid[:-1] + x_grid[1:]) / 2, (y_grid[:-1] + y_grid[1:]) / 2)

                for poly, val in zip(polygons, values):
                    # Get bounding box indices of polygon in raster
                    minx, miny, maxx, maxy = poly.bounds

                    x_inds = np.where((xv[0, :] >= minx) & (xv[0, :] <= maxx))[0]
                    y_inds = np.where((yv[:, 0] >= miny) & (yv[:, 0] <= maxy))[0]

                    for i in y_inds:
                        for j in x_inds:
                            pt = Point(xv[i, j], yv[i, j])
                            if poly.contains(pt):
                                heatmap[i, j] = val

                # now we want to create one heatmap per granularity level
                for granularity_level in self.__granularity_levels:
                    # generalize if not the finest granularity level
                    if granularity_level != self.__granularity_levels[0]:
                        heatmap= self.make_heatmap_coarser(heatmap)

                    # blur heatmap
                    heatmap = self.blur_heatmap(heatmap)

                    # mask with country border
                    heatmap_masked = self.mask_heatmap_with_country_border(heatmap, country_polygon, min_x, min_y, max_x, max_y)

                    # Save Heatmap to file
                    self.write_heatmap_to_file(heatmap_masked, [min_x,max_x], [min_y,max_y], "../vis_tool/data/" + country_name + "_" + perspective + "_heatmap_" + mode_titles[mode] + "_" + granularity_level + ".csv")

    def blur_heatmap(self, heatmap, sigma=1):
        mask = ~np.isnan(heatmap)
        data_filled = np.where(mask, heatmap, 0)
        blurred = gaussian_filter(data_filled, sigma=sigma, mode='reflect')
        mask_blurred = gaussian_filter(mask.astype(float), sigma=sigma, mode='reflect')
        # Normalize, keep NaNs outside
        with np.errstate(invalid='ignore', divide='ignore'):
            heatmap_blurred = np.where(mask_blurred > 0, blurred / mask_blurred, np.nan)

        return heatmap_blurred

    def mask_heatmap_with_country_border(self, heatmap, country_polygon, min_x, min_y, max_x, max_y):
        heatmap_local = heatmap.copy()
        ny, nx = heatmap_local.shape
        xs = np.linspace(min_x, max_x, nx)
        ys = np.linspace(min_y, max_y, ny)
        xv, yv = np.meshgrid(xs, ys)

        pts = [Point(x, y) for x, y in zip(xv.ravel(), yv.ravel())]
        inside_mask = np.array([country_polygon.contains(pt) for pt in pts], dtype=bool).reshape(ny, nx)

        heatmap_local[~inside_mask] = np.nan
        return heatmap_local

    def make_heatmap_coarser(self, heatmap):
        n = heatmap.shape[0]
        if(n%self.__ratio_for_coarser_resolutions!=0):
            raise Exception("Heatmap size is not divisible by ratio for coarser resolutions.")
        # reshape into (n//ratio, ratio, n//ratio, ratio)
        coarse = heatmap.reshape(n // self.__ratio_for_coarser_resolutions, self.__ratio_for_coarser_resolutions, n // self.__ratio_for_coarser_resolutions, self.__ratio_for_coarser_resolutions)
        # average over the small blocks (axes 1 and 3)
        return np.nanmean(np.nanmean(coarse, axis=3), axis=1)
