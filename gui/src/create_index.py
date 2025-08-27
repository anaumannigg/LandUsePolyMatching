import os
import json


def generate_index(base_dir):
    DATA_FOLDER = os.path.join(base_dir, "data")  # data folder inside the passed dir
    INDEX_FILE = os.path.join(DATA_FOLDER, "index.json")
    index = {}

    for filename in os.listdir(DATA_FOLDER):
        if filename.endswith(".csv"):
            name = filename[:-4]  # remove .csv
            parts = name.split("_")

            if len(parts) < 4:
                print("Skipping unexpected filename: {}".format(filename))
                continue

            # Only take the first four parts, ignore additional granularity
            dataset, perspective, vis_type, measure = parts[0], parts[1], parts[2], parts[3]

            if dataset not in index:
                index[dataset] = {}

            if perspective not in index[dataset]:
                index[dataset][perspective] = {
                    "visType": vis_type,
                    "measures": []
                }

            if measure not in index[dataset][perspective]["measures"]:
                index[dataset][perspective]["measures"].append(measure)

    with open(INDEX_FILE, "w", encoding="utf-8") as f:
        json.dump(index, f, indent=2)

    print("Index.json generated with {} datasets".format(len(index)))

