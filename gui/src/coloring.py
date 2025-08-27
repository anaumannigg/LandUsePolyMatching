import numpy as np
import pandas as pd

MAX_MULTIPLICITY = 5

def fill_empty_grid_cells(gdf, value_col):
    """
    For rows in gdf where value_col is NaN, fill with the average
    of neighboring cells' values (cells whose geometries touch this cell).
    """
    gdf = gdf.copy()  # work on a copy to avoid warnings

    # Identify which rows need filling
    missing_mask = gdf[value_col].isna()

    # For performance: build spatial index on geometries
    geoms = list(gdf.geometry)
    from shapely.strtree import STRtree
    tree = STRtree(geoms)

    # Map geometry id to index for quick lookup
    geom_id_to_idx = {id(geom): idx for idx, geom in enumerate(geoms)}

    # Loop over cells with missing value
    for idx in gdf[missing_mask].index:
        geom = gdf.at[idx, 'geometry']

        # Query spatial index for candidates that intersect bbox
        candidates = tree.query(geom)

        # Find neighbors that actually touch this geometry (excluding self)
        neighbor_indices = [
            i for i in candidates
            if i != idx and geoms[i].touches(geom)
        ]

        # Get neighbor values, drop missing neighbors
        neighbor_vals = gdf.loc[neighbor_indices, value_col].dropna()

        if not neighbor_vals.empty:
            # Set missing value to average of neighbors
            gdf.at[idx, value_col] = neighbor_vals.mean()

    return gdf

def apply_visualization(gdf, mode):
    gdf = gdf.copy()

    if mode == "quality":
        attribute = 'weight'
        invert = False

        # Normalize between 0-1
        values = gdf[attribute]
        norm = (values - values.min()) / (values.max() - values.min()) if values.max() != values.min() else values * 0


    elif mode == "moeller_score":
        attribute = 'moeller_score'
        invert = False

        # Normalize between 0-1
        values = gdf[attribute]
        norm = (values - values.min()) / (values.max() - values.min()) if values.max() != values.min() else values * 0


    elif mode in ["mult","mult_this", "mult_other"]:
        attribute = "mult_other"
        invert = True  # if you want to invert these

        values = gdf[attribute]
        norm = values.copy()  # use raw values (not normalized)

        #clip to maximum multiplicity
        norm[norm > MAX_MULTIPLICITY] = MAX_MULTIPLICITY

        if invert:
            norm = norm.max() - norm  # simple inversion of raw values
            norm = (norm - norm.min()) / (
                        norm.max() - norm.min()) if norm.max() != norm.min() else norm * 0 # normalize to values between 0 and 1

    elif mode == "unmatched":
        attribute = "is_matched"

        # Normalize between 0-1
        values = gdf[attribute]
        norm = (values - values.min()) / (values.max() - values.min()) if values.max() != values.min() else values * 0

    elif mode == "is_matched_correctly":
        attribute = "is_matched_correctly"

        values = pd.to_numeric(gdf[attribute], errors="coerce")
        norm = (values - values.min()) / (values.max() - values.min()) if values.max() != values.min() else values * 0

    else:
        raise ValueError(f"Unknown visualization mode: {mode}")

    # Optional: fill missing norm values before coloring
    gdf['__norm'] = norm
    gdf = fill_empty_grid_cells(gdf, '__norm')

    #for heatmap visualization
    gdf['color_value'] = gdf['__norm']

    return gdf.drop(columns=['__norm'])

