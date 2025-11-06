import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
import requests
from io import StringIO
from Bio import SeqIO

INFLUENZA_ACCESSIONS = [
    "CY121680", "CY121681", "CY121682", "CY121683", "CY121684",
    "CY121685", "CY121686", "CY121687", "CY013986", "CY013987"
]

ECO_RI = "GAATTC" 
def download_fasta_by_accession(acc):
    """
    Download a sequence in FASTA form using NCBI viewer.
    Returns a plain A/C/G/T string (uppercase).
    """
    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={acc}&db=nuccore&report=fasta"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    fasta_text = r.text
    record = SeqIO.read(StringIO(fasta_text), "fasta")
    # Convert RNA to DNA if needed (U->T)
    return str(record.seq).upper().replace("U", "T")

def digest_ecori(seq):
    """
    Cut DNA sequence at EcoRI (GAATTC).
    Return list of fragment lengths (bp).
    """
    cut = ECO_RI
    cut_len = len(cut)
    positions = []
    i = seq.find(cut)
    while i != -1:
        positions.append(i)
        i = seq.find(cut, i + 1)

    if not positions:
        return [len(seq)]

    frags = []
    start = 0
    for pos in positions:
        frags.append(pos - start)          
        start = pos + cut_len             
    frags.append(len(seq) - start)        
    return [f for f in frags if f > 0]

def migration_from_sizes(bp_lengths, a=1.0, b=0.32, noise_sd=0.012):
    """
    Very simple agarose migration model:
    distance = a - b*log10(size_bp) + small noise
    Returned values are in [~0..1] top->bottom.
    """
    sizes = np.array(bp_lengths, dtype=float)
    sizes = sizes[sizes > 0]
    if sizes.size == 0:
        return np.array([])
    y = a - b * np.log10(sizes)
    y += np.random.normal(0, noise_sd, size=y.shape[0])
    return y

class InfluenzaEcoRIGUI:
    def __init__(self, root):
        self.root = root
        root.title("Influenza EcoRI Digest — Gel Simulation")
        root.geometry("1000x800")

        # Top controls
        title = tk.Label(root, text="Influenza Genomes → EcoRI Digest → Gel Simulation", font=("Segoe UI", 16, "bold"))
        title.pack(pady=10)

        btn_frame = tk.Frame(root)
        btn_frame.pack(pady=6)

        self.download_btn = tk.Button(
            btn_frame, text="Download 10 genomes & Digest (EcoRI)",
            command=self.download_and_digest, font=("Segoe UI", 11)
        )
        self.download_btn.grid(row=0, column=0, padx=5)

        self.combined_btn = tk.Button(
            btn_frame, text="Show Combined Gel",
            command=self.show_combined_gel, state="disabled", font=("Segoe UI", 11)
        )
        self.combined_btn.grid(row=0, column=1, padx=5)

        # Select genome
        self.selected_acc = tk.StringVar()
        self.acc_combo = ttk.Combobox(
            btn_frame, textvariable=self.selected_acc, state="readonly", width=20, font=("Segoe UI", 10)
        )
        self.acc_combo.grid(row=0, column=2, padx=10)
        self.acc_combo.bind("<<ComboboxSelected>>", lambda e: self.show_single_gel())

        self.single_btn = tk.Button(
            btn_frame, text="Show Selected Gel", command=self.show_single_gel,
            state="disabled", font=("Segoe UI", 11)
        )
        self.single_btn.grid(row=0, column=3, padx=5)

        # Info labels
        self.status = tk.Label(root, text="Status: idle", anchor="w")
        self.status.pack(fill="x", padx=10, pady=5)

        self.most_label = tk.Label(root, text="Most fragments: —", anchor="w", font=("Segoe UI", 10, "bold"))
        self.most_label.pack(fill="x", padx=10)

        # Table of fragment counts
        table_frame = tk.Frame(root)
        table_frame.pack(pady=6)
        self.tree = ttk.Treeview(
            table_frame, columns=("Fragments", "Length (bp, summary)"), show="headings", height=8
        )
        self.tree.column("Fragments", width=120, anchor="center")
        self.tree.column("Length (bp, summary)", width=420, anchor="w")
        self.tree.heading("Fragments", text="# Fragments")
        self.tree.heading("Length (bp, summary)", text="Fragment sizes (first 10 shown)")
        self.tree.pack()

        # Matplotlib figure
        self.fig, self.ax = plt.subplots(figsize=(6, 7), dpi=120)
        self.canvas_frame = tk.Frame(root)
        self.canvas_frame.pack(fill="both", expand=True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # Data containers
        self.sequences = {}        
        self.fragments = {}        
        self.migrations = {}       

    def download_and_digest(self):
        self.status.config(text="Status: downloading sequences from NCBI…")
        self.root.update_idletasks()

        self.sequences.clear()
        self.fragments.clear()
        self.migrations.clear()

        errors = []
        for acc in INFLUENZA_ACCESSIONS:
            try:
                seq = download_fasta_by_accession(acc)
                self.sequences[acc] = seq
                frags = digest_ecori(seq)
                self.fragments[acc] = frags
                self.migrations[acc] = migration_from_sizes(frags)
                self.status.config(text=f"Downloaded + digested {acc} — {len(frags)} fragments")
                self.root.update_idletasks()
            except Exception as e:
                errors.append(f"{acc}: {e}")

        if errors:
            messagebox.showwarning("Download issues", "Some accessions failed:\n\n" + "\n".join(errors))

        if not self.fragments:
            self.status.config(text="No genomes processed.")
            return

        # Update UI elements
        self.acc_combo["values"] = list(self.fragments.keys())
        if not self.selected_acc.get() and self.acc_combo["values"]:
            self.selected_acc.set(self.acc_combo["values"][0])

        self.combined_btn.config(state="normal")
        self.single_btn.config(state="normal")

        # Update table + "most fragments"
        self.refresh_table_and_most()

        self.status.config(text="Done. Use the buttons to view gels.")

    def refresh_table_and_most(self):
        # clear table
        for row in self.tree.get_children():
            self.tree.delete(row)

        most_acc = None
        most_n = -1
        for acc in self.fragments:
            frag_list = self.fragments[acc]
            n = len(frag_list)
            if n > most_n:
                most_n = n
                most_acc = acc
            summary = ", ".join(str(x) for x in frag_list[:10])
            if len(frag_list) > 10:
                summary += ", …"
            self.tree.insert("", "end", values=(n, summary), text=acc, tags=(acc,))

        for iid in self.tree.get_children():
            acc = self.tree.item(iid, "tags")[0]
            self.tree.item(iid, text=acc)

        if most_acc is not None:
            self.most_label.config(text=f"Most fragments: {most_acc}  ({most_n} fragments)")


    def _draw_gel_axes(self, title, lanes):
        """lanes = list of (x_center, [y_band_positions])"""
        self.ax.clear()
        self.ax.set_title(title)
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1.05)
        self.ax.set_xticks([x for x, _ in lanes])
        self.ax.set_xticklabels([lbl for lbl, _ in zip([f"L{i+1}" for i in range(len(lanes))], lanes)], rotation=0)
        self.ax.set_yticks([])

        lane_width = 0.10
        for x, bands in lanes:
            # lane background
            self.ax.axvspan(x - lane_width / 2, x + lane_width / 2, alpha=0.25)
            # bands
            for y in bands:
                self.ax.hlines(y, x - lane_width / 2 + 0.01, x + lane_width / 2 - 0.01, linewidth=6)

    def show_combined_gel(self):
        if not self.migrations:
            return
        accs = list(self.migrations.keys())
        xs = np.linspace(0.1, 0.9, len(accs))
        lanes = [(x, self.migrations[acc]) for x, acc in zip(xs, accs)]
        self._draw_gel_axes("Combined gel (EcoRI digest of 10 influenza genomes)", lanes)
        # Set accession labels below ticks
        self.ax.set_xticks(xs)
        self.ax.set_xticklabels(accs, rotation=90, fontsize=8)
        self.canvas.draw()

    def show_single_gel(self):
        acc = self.selected_acc.get()
        if not acc:
            return
        bands = self.migrations.get(acc, [])
        lanes = [(0.5, bands)]
        self._draw_gel_axes(f"EcoRI digest — {acc}", lanes)
        label_bp = [3000, 1500, 500]
        a, b = 1.0, 0.32
        label_pos = a - b * np.log10(np.array(label_bp, dtype=float))
        for bp, y in zip(label_bp, label_pos):
            self.ax.text(0.02, y, f"{bp} bp", va="center", fontsize=8)
        self.ax.set_xticks([0.5])
        self.ax.set_xticklabels(["Sample"])
        self.canvas.draw()

\
if __name__ == "__main__":
    root = tk.Tk()
    app = InfluenzaEcoRIGUI(root)
    root.mainloop()
