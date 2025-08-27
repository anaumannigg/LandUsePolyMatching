import tkinter as tk
from tkinter import ttk, filedialog
import os
import threading
from src.computing import Compute
import subprocess
import itertools
from src.matching_predecomposer import compute_connected_components
from src.create_index import generate_index
import webbrowser
import socketserver
import http.server
import re

class Application:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("MatchingVis Data Preparation")
        self.root.geometry("1080x900")

        # Directory dropdown for "../input/"
        self.valid_dirs = self.find_valid_input_dirs()
        self.selected_input_dir = tk.StringVar()
        if self.valid_dirs:
            self.selected_input_dir.set(self.valid_dirs[0])

        self.selected_dir = tk.StringVar()
        self.subdirs = []
        self.subdir_var = tk.StringVar()

        self.status_message = tk.StringVar()
        self.status_message.set("Ready.")

        self.vis_server = None
        self.vis_thread = None
        self.vis_running = False

        self.create_widgets()

    def find_valid_input_dirs(self):
        """Scan ../input/ for directories with at least 2 gpkg files
        named <dirname>_xyz.gpkg"""
        base_path = os.path.join("..", "input")
        valid = []
        if os.path.isdir(base_path):
            for d in os.listdir(base_path):
                dir_path = os.path.join(base_path, d)
                if os.path.isdir(dir_path):
                    gpkg_files = [f for f in os.listdir(dir_path)
                                  if f.startswith(d + "_") and f.endswith(".gpkg")]
                    if len(gpkg_files) >= 2:
                        valid.append(d)
        return sorted(valid)

    def get_descriptor_pairs(self, dataset):
        """Return all possible descriptor pairs for a dataset, alphabetically sorted"""
        dataset_dir = os.path.join("..", "input", dataset)
        if not os.path.isdir(dataset_dir):
            return []

        # Extract descriptors from filenames of form dataset_DESCRIPTOR.gpkg
        descriptors = []
        for f in os.listdir(dataset_dir):
            if f.startswith(dataset + "_") and f.endswith(".gpkg"):
                desc = f[len(dataset) + 1: -5]  # strip "dataset_" and ".gpkg"
                descriptors.append(desc)

        descriptors = sorted(set(descriptors))
        # Build all unique pairs (combinations of 2 descriptors)
        pairs = [f"{a} - {b}" for a, b in itertools.combinations(descriptors, 2)]
        return pairs

    def on_dataset_selected(self, event=None):
        """Update descriptor pairs when dataset changes"""
        dataset = self.selected_input_dir.get()
        pairs = self.get_descriptor_pairs(dataset)
        self.pair_menu["values"] = pairs
        if pairs:
            self.pair_var.set(pairs[0])
        else:
            self.pair_var.set("")

    def create_widgets(self):
        # --- Matching section ---
        tk.Label(self.root, text="--- Matching ---").pack(pady=(10, 2))
        tk.Label(self.root, text="Select Input Dataset:").pack(pady=(10, 2))
        self.input_menu = ttk.Combobox(
            self.root, textvariable=self.selected_input_dir,
            values=self.valid_dirs, state="readonly", width=50
        )
        self.input_menu.pack(pady=5)
        self.input_menu.bind("<<ComboboxSelected>>", self.on_dataset_selected)

        tk.Label(self.root, text="Select Descriptor Pair:").pack(pady=(10, 2))
        self.pair_var = tk.StringVar()
        self.pair_menu = ttk.Combobox(
            self.root, textvariable=self.pair_var,
            values=[], state="readonly", width=50
        )
        self.pair_menu.pack(pady=5)

        if self.valid_dirs:
            self.selected_input_dir.set(self.valid_dirs[0])
            self.on_dataset_selected(None)

        compute_btn = tk.Button(self.root, text="Compute Matching", command=self.compute_matching_thread)
        compute_btn.pack(pady=5)

        separator = ttk.Separator(self.root, orient="horizontal")
        separator.pack(fill="x", pady=10)

        # --- Prepare Heatmaps section ---
        tk.Label(self.root, text="--- Prepare Heatmaps ---").pack(pady=(10, 2))

        tk.Label(self.root, text="Select Export Dataset:").pack(pady=(10, 2))
        self.valid_export_dirs = self.find_valid_export_dirs()
        self.selected_export_dir = tk.StringVar()

        self.export_menu = ttk.Combobox(
            self.root, textvariable=self.selected_export_dir,
            values=self.valid_export_dirs, state="readonly", width=50
        )
        self.export_menu.pack(pady=5)
        self.export_menu.bind("<<ComboboxSelected>>", self.on_export_dataset_selected)

        tk.Label(self.root, text="Select Lambda Subdir:").pack(pady=(10, 2))
        self.lambda_var = tk.StringVar()
        self.lambda_menu = ttk.Combobox(
            self.root, textvariable=self.lambda_var,
            values=[], state="readonly", width=50
        )
        self.lambda_menu.pack(pady=5)

        if self.valid_export_dirs:
            self.selected_export_dir.set(self.valid_export_dirs[0])
            self.on_export_dataset_selected(None)

        run_btn = tk.Button(self.root, text="Run", command=self.run_tasks_thread)
        run_btn.pack(pady=5)

        separator = ttk.Separator(self.root, orient="horizontal")
        separator.pack(fill="x", pady=10)

        # --- Launch Visualization Tool section ---
        tk.Label(self.root, text="--- Launch Visualization Tool ---").pack(pady=(10, 2))
        self.vis_btn = tk.Button(self.root, text="Launch Visualization Tool", command=self.toggle_vis_tool)
        self.vis_btn.pack(pady=5)

        separator = ttk.Separator(self.root, orient="horizontal")
        separator.pack(fill="x", pady=10)


        # --- Progress section ---
        tk.Label(self.root, text="--- Progress of your current operation ---").pack(pady=(10, 2))
        self.progress = ttk.Progressbar(self.root, orient="horizontal", length=400, mode="determinate")
        self.progress.pack(pady=10)

        status_label = tk.Label(self.root, textvariable=self.status_message)
        status_label.pack(pady=5)

        separator = ttk.Separator(self.root, orient="horizontal")
        separator.pack(fill="x", pady=10)

        # Console
        tk.Label(self.root, text="--- Console Output of C++ Tool ---").pack(pady=(5, 2))
        self.console = tk.Text(self.root, height=3, width=100, bg="black", fg="white")
        self.console.pack(pady=2, fill="x")
        self.console.config(state="disabled")  # start read-only

        close_btn = tk.Button(self.root, text="Close", command=self.close)
        close_btn.pack(pady=5)

    def toggle_vis_tool(self):
        if not self.vis_running:
            self.start_vis_tool()
        else:
            self.stop_vis_tool()

    def start_vis_tool(self):
        vis_dir = "../vis_tool/"
        generate_index(vis_dir)

        PORT = 8000
        Handler = http.server.SimpleHTTPRequestHandler
        os.chdir(vis_dir)

        class ReusableTCPServer(socketserver.TCPServer):
            allow_reuse_address = True

        self.vis_server = ReusableTCPServer(("", PORT), Handler)

        self.vis_thread = threading.Thread(target=self.vis_server.serve_forever, daemon=True)
        self.vis_thread.start()

        self.status_message.set(f"Visualization running on http://localhost:{PORT}")
        webbrowser.open(f"http://localhost:{PORT}")

        self.vis_running = True
        self.vis_btn.config(text="Stop Visualization Tool")

    def stop_vis_tool(self):
        if self.vis_server:
            self.vis_server.shutdown()
            self.vis_server.server_close()
            self.vis_thread.join()  # wait for thread to finish
            self.vis_server = None
            self.vis_thread = None

        self.status_message.set("Visualization server stopped.")
        self.vis_running = False
        self.vis_btn.config(text="Run Visualization Tool")

    def compute_matching_thread(self):
        thread = threading.Thread(target=self.compute_matching)
        thread.start()

    def compute_matching(self):
        dataset = self.selected_input_dir.get()
        if not dataset:
            self.status_message.set("No valid dataset selected.")
            return

        pair = self.pair_var.get()
        if not pair or " - " not in pair:
            self.status_message.set("No valid descriptor pair selected.")
            return

        # Split "DESCRIPTOR1 - DESCRIPTOR2"
        dd1, dd2 = pair.split(" - ")

        build_dir = os.path.join("..", "build")
        input_dir = os.path.join("..", "input", dataset)
        connected_file = os.path.join(input_dir, "connected_components.txt")

        # Build step list
        steps = []

        # Step 0: optional connected components
        if not os.path.exists(connected_file):
            fileA = os.path.join(input_dir, f"{dataset}_{dd1}.gpkg")
            fileB = os.path.join(input_dir, f"{dataset}_{dd2}.gpkg")
            steps.append((
                None,
                f"Computing connected components for {dd1} and {dd2}...",
                lambda: compute_connected_components(fileA, fileB)
            ))

        # Normal C++ build + run pipeline
        steps.extend([
            (["mkdir"], "Creating build directory...", lambda: os.makedirs(build_dir, exist_ok=True)),
            (["cmake", "-DCMAKE_BUILD_TYPE=Release", ".."], "Running CMake...", None),
            (["make", "-j", "8"], "Building project...", None),
            (
                ["./TCPolygonMatching", "-d", dataset, "-t", "12", "-l", "0.3", "-dd1", dd1, "-dd2", dd2],
                f"Running Matching on {dd1} vs {dd2}...",
                None
            ),
        ])

        # Configure progress bar
        self.progress["maximum"] = len(steps)
        self.progress["value"] = 0

        try:
            for i, (cmd, msg, py_action) in enumerate(steps, start=1):
                self.status_message.set(msg)
                self.root.update_idletasks()

                if py_action:  # Python action (mkdir or compute_connected_components)
                    py_action()
                else:
                    # Run subprocess and capture output
                    proc = subprocess.Popen(
                        cmd, cwd=build_dir,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True
                    )
                    for line in proc.stdout:
                        self.console_write(line.rstrip())
                    proc.wait()
                    if proc.returncode != 0:
                        raise subprocess.CalledProcessError(proc.returncode, cmd)

                self.progress["value"] = i
                self.root.update_idletasks()

            self.status_message.set("Matching completed successfully.")

        except subprocess.CalledProcessError as e:
            self.status_message.set(f"Command failed: {e}")
        except Exception as e:
            self.status_message.set(f"Error during matching: {e}")

        #refresh heatmap directories
        self.valid_export_dirs = self.find_valid_export_dirs()
        if self.valid_export_dirs:
            self.export_menu["values"] = self.valid_export_dirs
            self.selected_export_dir.set(self.valid_export_dirs[0])
            self.on_export_dataset_selected(None)

    def load_directory(self):
        folder = filedialog.askdirectory()
        if folder:
            self.selected_dir.set(folder)
            self.update_subdirs()

    def find_valid_export_dirs(self):
        """Scan ../export/ for datasets with a csv + at least one lambda subdir containing 2 gpkg files"""
        base_path = os.path.join("..", "export")
        valid = []
        if os.path.isdir(base_path):
            for d in os.listdir(base_path):
                dir_path = os.path.join(base_path, d)
                if os.path.isdir(dir_path):
                    # must contain at least one lambda* subdir with 2 gpkg files
                    lambda_dirs = [
                        sub for sub in os.listdir(dir_path)
                        if sub.startswith("lambda") and os.path.isdir(os.path.join(dir_path, sub))
                    ]
                    valid_lambda = False
                    for ld in lambda_dirs:
                        gpkg_files = [f for f in os.listdir(os.path.join(dir_path, ld)) if f.endswith(".gpkg")]
                        # must contain [dataset]_data_per_match.csv
                        csv_file = os.path.join(dir_path, ld, f"{d}_data_per_match.csv")
                        if len(gpkg_files) >= 2 and os.path.isfile(csv_file):
                            valid_lambda = True
                            break
                    if valid_lambda:
                        valid.append(d)
        return sorted(valid)

    def get_lambda_subdirs(self, dataset):
        """Return all lambda* subdirs of a dataset that contain 2 gpkg files"""
        dataset_dir = os.path.join("..", "export", dataset)
        lambda_dirs = []
        if os.path.isdir(dataset_dir):
            for sub in os.listdir(dataset_dir):
                subdir_path = os.path.join(dataset_dir, sub)
                if sub.startswith("lambda") and os.path.isdir(subdir_path):
                    gpkg_files = [f for f in os.listdir(subdir_path) if f.endswith(".gpkg")]
                    if len(gpkg_files) >= 2:
                        lambda_dirs.append(sub)
        return sorted(lambda_dirs)

    def on_export_dataset_selected(self, event=None):
        """Update lambda subdirs when export dataset changes"""
        dataset = self.selected_export_dir.get()
        lambdas = self.get_lambda_subdirs(dataset)
        self.lambda_menu["values"] = lambdas
        if lambdas:
            self.lambda_var.set(lambdas[0])
        else:
            self.lambda_var.set("")

    def run_tasks_thread(self):
        """Run tasks in a separate thread so UI doesn't freeze."""
        thread = threading.Thread(target=self.run_tasks)
        thread.start()

    def run_tasks(self):
        dir_path = "../export/"  + self.selected_export_dir.get()
        subdir_name = self.lambda_var.get()
        subdir_path = os.path.join(dir_path, subdir_name)
        print("dir path: ",dir_path, " | subdir path: ",subdir_path)
        comp = Compute(dir_path, subdir_path)

        steps = [
            (lambda: comp.load_files(),"Loading files..."),
            (lambda: comp.mark_unmatched_polygons(),"Marking unmatched polygons..."),
            (lambda: comp.compute_moeller_score(),"Computing Moeller Score..."),
            (lambda: comp.compute_consistencies(),"Checking for consistency..."),
            (lambda: comp.overwrite_geopackages(),"Saving updated files..."),
            (lambda: comp.generate_heatmaps(),"Generating heatmaps...")
        ]

        self.progress["maximum"] = len(steps)
        self.progress["value"] = 0

        for i, (step,message) in enumerate(steps, start=1):
            self.status_message.set(message)
            step()
            self.progress["value"] = i
            self.root.update_idletasks()

        self.status_message.set("All tasks completed.")

    def run(self):
        self.root.mainloop()

    def console_write(self, msg):
        # Remove all ANSI escape sequences first
        ansi_escape = re.compile(r'\x1b\[[0-9;]*[A-Za-z]')
        msg_clean = ansi_escape.sub('', msg)

        if "\x1b[1A" in msg or "\x1b[2K" in msg:
            # Delete the last line and replace it
            self.console.config(state="normal")
            last_line_start = self.console.index("end-2l linestart")
            last_line_end = self.console.index("end-2l lineend+1c")
            self.console.delete(last_line_start, last_line_end)
            self.console.insert(last_line_start, msg_clean + "\n")
            self.console.see("end")
            self.console.update_idletasks()
            self.console.config(state="disabled")
        else:
            # Normal append
            self.console.config(state="normal")
            self.console.insert("end", msg_clean + "\n")
            self.console.see("end")
            self.console.update_idletasks()
            self.console.config(state="disabled")

    def close(self):
        if self.vis_running:
            self.stop_vis_tool()
        self.root.quit()
