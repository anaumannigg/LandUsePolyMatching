import geopandas as gpd
import pandas as pd
import time
from shapely.geometry import box
import os

def generate_grid_over_polygon(polygon_gdf, cell_size):
    """
    Creates a grid (GeoDataFrame) of square cells over the bounding box of the input polygon_gdf,
    clipped to the actual boundary.
    `cell_size` is in map units (e.g., meters for EPSG:3857)
    """
    # Use EPSG:3857 for projected coordinates
    if polygon_gdf.crs != "EPSG:3857":
        polygon_gdf = polygon_gdf.to_crs("EPSG:3857")

    bounds = polygon_gdf.total_bounds  # [minx, miny, maxx, maxy]
    xmin, ymin, xmax, ymax = bounds

    rows = int((ymax - ymin) // cell_size) + 1
    cols = int((xmax - xmin) // cell_size) + 1

    cells = []
    for i in range(cols):
        for j in range(rows):
            x0 = xmin + i * cell_size
            y0 = ymin + j * cell_size
            x1 = x0 + cell_size
            y1 = y0 + cell_size
            cell = box(x0, y0, x1, y1)
            cells.append(cell)

    grid = gpd.GeoDataFrame(geometry=cells, crs="EPSG:3857")
    clipped = gpd.overlay(grid, polygon_gdf, how="intersection")
    clipped["cell_id"] = clipped.index
    return clipped


def aggregate_attributes_to_grid(data_gdf, grid_gdf):
    """
    Aggregates average values of specified attributes from data_gdf into each grid cell.
    """
    joined = gpd.sjoin(grid_gdf, data_gdf, how="left", predicate="intersects")

    def mean_exclude_zeros(series):
        return series[series != 0].mean()

    # Base aggregation
    agg_dict = {
        "geometry": "first",  # grid cell geometry
        "weight": "mean",
        "mult_this": mean_exclude_zeros,
        "mult_other": mean_exclude_zeros,
        "moeller_score": "mean",
        "is_matched": "mean"
    }

    # Conditionally include 'is_matched_correctly' if it exists
    if "is_matched_correctly" in joined.columns:
        agg_dict["is_matched_correctly"] = "mean"

    aggregated = joined.groupby("cell_id").agg(agg_dict).reset_index()

    result = gpd.GeoDataFrame(aggregated, geometry="geometry", crs=grid_gdf.crs)
    return result


def detect_country_from_filename(filepath, country_boundary_dir):
    """
    Detects the country name based on the filename.
    Matches if a country boundary file name is contained in the filename.

    Parameters:
    - filepath: Path to the loaded GPKG file
    - country_boundary_dir: Directory containing country boundary files (e.g. GeoJSON)

    Returns:
    - country name (without extension) if matched, else None
    """

    filename = os.path.basename(filepath).lower()

    # List all country boundary files
    for fname in os.listdir(country_boundary_dir):
        if not fname.lower().endswith((".geojson", ".gpkg", ".shp")):
            continue  # skip irrelevant files

        country_name = os.path.splitext(fname)[0].lower()

        # If country name is a substring in the filename
        if country_name in filename:
            return country_name  # Match found

    return None  # No match


def load_file_and_grid(filepath, country_boundary_dir, cell_size=20000):
    gdf = gpd.read_file(filepath)
    gdf = gdf.to_crs("EPSG:3857")  # for performance/grid

    # Step 1: detect country from filename
    country_name = detect_country_from_filename(filepath, country_boundary_dir)

    if country_name:
        boundary_path = f"{country_boundary_dir}/{country_name}.gpkg"
        country_gdf = gpd.read_file(boundary_path).to_crs("EPSG:3857")

        # Step 2: generate grid
        grid = generate_grid_over_polygon(country_gdf, cell_size)

        # Step 3: aggregate
        aggregated = aggregate_attributes_to_grid(gdf, grid)
    else:
        aggregated = None  # fallback if no country match

    return aggregated, country_gdf
