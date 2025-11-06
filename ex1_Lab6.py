import random
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# ============================================================
#   INSERT YOUR FASTA SEQUENCE HERE (OPTIONAL)
# ============================================================

dna_sequence_user = """
TCTTTCTTTCTTTCTTTCTTTCCTTGTTTCTTTCTGTCTTTCATTACACTTTTTGCCTTTTTTTTCAGTC
TTAGAGAATTTTATATTAATATGTAGATTCTGTATATCTGCCAATATATCAGTTTTCGGCTTCCAATATT
AAAGAAATATCGGTTATCGTTATCGGCCAAAATTTACATATTGGTGCATTCCTAGATTCTATTGTGATTT
TGAACCTTCCTATAGAAAATATTATAGAATACTGTATTTACTAAACACTATATCAATTTACTATATTTTA
TTCTACAGAGTACTAGCACAACTTATAGCACTATAGTACTCTAATAAATGATATACTACAGTATGTGGGT
CAGCTTAAATGAAACCTCCAGTCATTGAGGATCGTTGTTCATCTTGAACAATGAGCTTGAAATTGGTGAA
TTTTTTGGGGCACTGTTCATCTTTTCTCAGGTGTAATATTGTCAGTCTTGGCCAAATCTCCATGGCGACA
GTGTCAGAGGTGTGATAATGTTCTAGAGGGAACCAGCAGGTCCTAATGACGATGGGACGATTGTGCCGCA
CTCGCTGTGTTCTCTCATGTTCTGTTCTTCTGGGATCTTACATGTCTCTCCTCTCTCACCCTCATGAATA
GCGTGATGTTGATTGATGCTACACTTTTACACTGGCCTGAAATCACAGCAGCTCGGTTTAAATGTGAAAG
AGCGCTTGGTGTTTGTTTATCCCTCTTCAGTGAAAGATGTCTCTGAAAGATTCTCTAATGTCTCACTTAT
TTCTCTTAGACATCAATGAGATCAAATGGAGAGCAAATGGACACTTCTGCTTAAGTCTCATGTGCGTCTC
AGATATTTCTCCTGAGAAAAAGAGTCAGAAATGTCACAAATGTGTCTGTCATTAGAGTTTGAGAAGAGAT
GTCAACAATGTGTCAAACTGAGCTCAGATTTGTCAATTCTGCAATCAAAAGTCTATTAAAGATCTCAATA
TAAACTCAAACATTTCTCATCACATTTTTCTTTCATCTACACATGTTTAATTGAAGGAGCTGAATGAACA
AAATGGAGATCAATTGTTGATTTGAATAAGACAATGTTCATATCATATTTAATGTGATGACGGGGGTGGG
GGGGTGTTCGCTCTGTAATGAGCCATAATGCCCGGCCCTCCTGGGCAGAAACTAATTGGTGCAGGATTAG
CATTAGAAAACCCAGCGACCAGAGATTTTCTTCTTTCTGCCTATCAGTGACCAACTAAATGCCCAAAACA
CACCAGATGCCCAATTTAGAGGTAACCGGCCCACATTTCTCAACATTCATTTCCGGGGAAAACTTGAGTT
TTGGAATGTGAGTCTGAGTAATAATCTCTCATGGCCAATGGCAAAAATATAACCTGCTTATTCAAGACAG
AAAATCACTGCGGCTGTATGAGACGCCGGCAACAACAGGGACAAATCTGCTGTTAAAGAACAAAAATCCT
CTTAAATTCATGATTAAAAAAACTTTTGGTTGCAAGAAAATCACCTAAAATATTTTTAACTAATATACGG
TTTATGTACCAACCTATTTTTTACTGTTTTTCAGATTTTTGTTATTATGCACAGTACATAGAATTTGAAC
AATATTTATATAGTTAATAATTGTCATGTTTATTGTATAGATTAATTAGGAAATGTATATCTTTTACTTT
ACTTTTTT
"""
dna_sequence_user = dna_sequence_user.replace("\n", "").strip()  

def generate_dna():
    """Return either user-provided FASTA or random DNA"""
    if len(dna_sequence_user) > 0:
        return dna_sequence_user

    seq_len = random.randint(1000, 3000)
    return "".join(random.choice("ACGT") for _ in range(seq_len))


def generate_fragments(dna):
    """Generate 10 random fragments from the DNA"""
    fragments = []
    seq_len = len(dna)

    for i in range(10):
        frag_len = random.randint(100, seq_len)
        start = random.randint(0, seq_len - frag_len)
        end = start + frag_len

        fragments.append({
            "id": i+1,
            "start": start,
            "end": end,
            "length": frag_len
        })

    return fragments


def simulate_gel(fragments):
    """Compute migration distance based on log10(bp)"""
    a = 1.0
    b = 0.35

    sizes = np.array([f["length"] for f in fragments], dtype=float)
    log_bp = np.log10(sizes)
    noise = np.random.normal(0, 0.015, len(sizes))

    return a - b * log_bp + noise

class GelApp:
    def __init__(self, root):
        self.root = root
        root.title("Gel Electrophoresis Simulator")
        root.geometry("800x700")

        title = tk.Label(root, text="DNA Fragment Generator & Gel Simulation", font=("Arial", 16, "bold"))
        title.pack(pady=10)

        button = tk.Button(root, text="Generate DNA & Run Simulation", font=("Arial", 12), command=self.run_simulation)
        button.pack(pady=10)

        self.tree = ttk.Treeview(root, columns=("Start", "End", "Length"), show="headings", height=10)

        self.tree.column("Start", width=150, anchor="center")
        self.tree.column("End", width=150, anchor="center")
        self.tree.column("Length", width=150, anchor="center")

        self.tree.heading("Start", text="Start (0-based)")
        self.tree.heading("End", text="End")
        self.tree.heading("Length", text="Length (bp)")

        self.tree.pack(pady=10)
        self.canvas_frame = tk.Frame(root)
        self.canvas_frame.pack(fill="both", expand=True)

        self.fig, self.ax = plt.subplots(figsize=(4, 6), dpi=120)

    def run_simulation(self):
        dna = generate_dna()
        fragments = generate_fragments(dna)
        distances = simulate_gel(fragments)

        # Clear table
        for row in self.tree.get_children():
            self.tree.delete(row)

        # Insert new rows
        for f in fragments:
            self.tree.insert("", "end", values=(f["start"], f["end"], f["length"]))

        self.plot_gel(distances)

    def plot_gel(self, distances):

        # Remove old plot
        for widget in self.canvas_frame.winfo_children():
            widget.destroy()

        self.ax.clear()

        # Lane positions
        lane_width = 0.35
        marker_x = 0.3
        sample_x = 0.7

        # Ladder
        ladder_bp = np.array([3000, 2000, 1500, 1000, 700, 500, 300, 200, 100], dtype=float)
        a, b = 1.0, 0.35
        ladder_dist = a - b * np.log10(ladder_bp)

        # Draw lanes
        self.ax.axvspan(marker_x - lane_width/2, marker_x + lane_width/2, alpha=0.25)
        self.ax.axvspan(sample_x - lane_width/2, sample_x + lane_width/2, alpha=0.25)

        # Marker bands
        for y in ladder_dist:
            self.ax.hlines(y, marker_x - lane_width/2 + 0.02, marker_x + lane_width/2 - 0.02, linewidth=6)

        # Sample bands
        for y in distances:
            self.ax.hlines(y, sample_x - lane_width/2 + 0.02, sample_x + lane_width/2 - 0.02, linewidth=6)

        # Format
        self.ax.set_ylim(0, 1.05)
        self.ax.set_xlim(0, 1)
        self.ax.set_xticks([marker_x, sample_x])
        self.ax.set_xticklabels(["Marker", "Sample"])
        self.ax.set_yticks([])
        self.ax.set_title("Simulated Agarose Gel")

        # Draw into GUI
        canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)



root = tk.Tk()
app = GelApp(root)
root.mainloop()
