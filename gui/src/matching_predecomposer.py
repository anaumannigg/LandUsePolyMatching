import geopandas as gpd
import networkx as nx
import shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import os

# Parameters
OVERLAP_EPSILON = 1e-8  # minimum area to consider polygons as overlapping

def select_vector_files():
    root = tk.Tk()
    root.withdraw()
    file_paths = filedialog.askopenfilenames(title='Select TWO polygon shapefiles',
                                             filetypes=[("Vector files", "*.shp *.gpkg"), 
                                                        ("Shapefiles", "*.shp"),
                                                        ("GeoPackage", "*.gpkg")])
    return list(file_paths)

def select_output_file():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.asksaveasfilename(title='Save output shapefile as',
                                              defaultextension=".shp",
                                              filetypes=[("Shapefiles", "*.shp")])
    return file_path

def build_intersection_graph(gdf_a, gdf_b, epsilon=OVERLAP_EPSILON):
    G = nx.Graph()

    # Add nodes
    for i, geom in enumerate(gdf_a.geometry):
        G.add_node(('A', i), geometry=geom)

    for j, geom in enumerate(gdf_b.geometry):
        G.add_node(('B', j), geometry=geom)

    # Spatial index
    sindex_b = gdf_b.sindex

    for i, geom_a in enumerate(gdf_a.geometry):
        if geom_a is None or geom_a.is_empty:
            continue
        if not geom_a.is_valid:
            geom_a = geom_a.buffer(0)  # attempt fix
            if not geom_a.is_valid:
                continue  # skip if still invalid

        possible_matches_index = list(sindex_b.intersection(geom_a.bounds))
        for j in possible_matches_index:
            geom_b = gdf_b.geometry.iloc[j]
            if geom_b is None or geom_b.is_empty:
                continue
            if not geom_b.is_valid:
                geom_b = geom_b.buffer(0)
                if not geom_b.is_valid:
                    continue

            try:
                intersection = geom_a.intersection(geom_b)
                if intersection.is_empty or intersection.area < epsilon:
                    continue
                G.add_edge(('A', i), ('B', j))
            except shapely.errors.GEOSException as e:
                print(f"Skipped invalid intersection between A[{i}] and B[{j}]: {e}")
                continue

    return G

def compute_union_components(G):
    union_geometries = []
    for component in nx.connected_components(G):
        polygons = [G.nodes[node]['geometry'] for node in component]
        merged = unary_union(polygons)
        union_geometries.append(merged)
    return union_geometries

def write_components_txt(G, output_path):
    with open(output_path, 'w') as f:
        for component in nx.connected_components(G):
            a_ids = sorted(str(i) for node_type, i in component if node_type == 'A')
            b_ids = sorted(str(i) for node_type, i in component if node_type == 'B')
            line = ",".join(a_ids) + " " + ",".join(b_ids)
            f.write(line + "\n")

def compute_connected_components(fileA,fileB):
    #shapefiles = select_vector_files()
    #if len(shapefiles) != 2:
    #    print("Please select exactly two shapefiles.")
    #    return

    gdf_a = gpd.read_file(fileA)
    gdf_b = gpd.read_file(fileB)

    # Ensure same CRS
    if gdf_a.crs != gdf_b.crs:
        gdf_b = gdf_b.to_crs(gdf_a.crs)

    # Build intersection graph
    G = build_intersection_graph(gdf_a, gdf_b)

    # Write connected components to TXT
    common_dir = os.path.dirname(fileA)
    output_txt_path = os.path.join(common_dir, "connected_components.txt")
    write_components_txt(G, output_txt_path)
    print(f"Connected components written to: {output_txt_path}")

