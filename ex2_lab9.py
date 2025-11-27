import tkinter as tk
from tkinter import ttk, scrolledtext, filedialog, messagebox
import math
import os

ENZYMES = {
    "EcoRI":  {"site": "GAATTC", "offset": 1},
    "BamHI":  {"site": "GGATCC", "offset": 1},
    "HindIII":{"site": "AAGCTT", "offset": 1},
    "TaqI":   {"site": "TCGA",   "offset": 1},
    "HaeIII": {"site": "GGCC",   "offset": 2},
}

def clean_sequence(seq):
    """Remove FASTA headers and non-ACGT characters."""
    seq = seq.upper()
    return "".join([c for c in seq if c in "ACGT"])


def read_fasta_from_file(path):
    with open(path, "r") as f:
        seq = []
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return clean_sequence("".join(seq))


def find_cuts(seq, site, offset):
    cuts = []
    pos = seq.find(site)
    while pos != -1:
        cuts.append(pos + offset)
        pos = seq.find(site, pos + 1)
    return cuts


def fragments_from_cuts(length, cuts):
    if not cuts:
        return [length]
    positions = [0] + sorted(cuts) + [length]
    return [positions[i+1] - positions[i] for i in range(len(positions)-1)]


def digest_sequence(seq):
    """Digest with all enzymes together (combined digest)."""
    all_cuts = []
    for e in ENZYMES.values():
        all_cuts.extend(find_cuts(seq, e["site"], e["offset"]))
    all_cuts = sorted(set(all_cuts))
    return fragments_from_cuts(len(seq), all_cuts)


def simulate_gel(lanes, height=20):
    """
    lanes = [(label, [fragment sizes])]
    Returns string containing ASCII gel.
    """
    all_frag_sizes = [f for _, frags in lanes for f in frags]
    if not all_frag_sizes:
        return "No fragments to display.\n"

    max_s = max(all_frag_sizes)
    min_s = min(all_frag_sizes)

    gel = ""

    for row in range(height):
        line = ""
        for _, frags in lanes:
            line += "|"
            drawn = False
            for s in frags:
                y = (math.log10(max_s) - math.log10(s)) / (
                    math.log10(max_s) - math.log10(min_s) + 1e-9
                )
                r = int(y * (height - 1))
                if r == row:
                    drawn = True
                    break
            line += "*" if drawn else " "
            line += " "
        gel += line + "\n"

    gel += "\n"
    for i, (label, _) in enumerate(lanes):
        gel += f"Lane {i+1}: {label}\n"

    return gel



class MultiDigestGUI:
    def __init__(self, root):
        self.root = root
        root.title("Influenza Multi-Genome Restriction Digest Simulator")

        # title
        tk.Label(root, text="Load 10 Influenza Genome FASTA Files", font=("Arial", 14, "bold")).grid(row=0, column=0, sticky="w", padx=5)

        # load button
        self.btn_load = tk.Button(root, text="Load FASTA Files", command=self.load_fastas)
        self.btn_load.grid(row=1, column=0, sticky="w", padx=5, pady=5)

        # list of loaded files
        self.file_box = scrolledtext.ScrolledText(root, width=70, height=8)
        self.file_box.grid(row=2, column=0, columnspan=2, padx=5, pady=5)

        # run button
        self.btn_run = tk.Button(root, text="Run Digest for All", command=self.run_digest)
        self.btn_run.grid(row=3, column=0, sticky="w", padx=5, pady=5)

        # output box (results + merged gel)
        tk.Label(root, text="Output:", font=("Arial", 12, "bold")).grid(row=4, column=0, sticky="w", padx=5)
        self.out = scrolledtext.ScrolledText(root, width=100, height=30)
        self.out.grid(row=5, column=0, columnspan=2, padx=5, pady=5)

        self.fasta_paths = []

    def load_fastas(self):
        paths = filedialog.askopenfilenames(
            title="Select Influenza FASTA Files",
            filetypes=[("FASTA files", "*.fasta *.fa *.txt"), ("All files", "*.*")]
        )
        if not paths:
            return

        self.fasta_paths = list(paths)
        self.file_box.delete("1.0", tk.END)
        for p in self.fasta_paths:
            self.file_box.insert(tk.END, p + "\n")

    def run_digest(self):
        if not self.fasta_paths:
            messagebox.showerror("Error", "Please load FASTA files first.")
            return

        lanes = []
        output = []

        for path in self.fasta_paths:
            seq = read_fasta_from_file(path)
            frags = digest_sequence(seq)
            label = os.path.basename(path)

            output.append(f"=== {label} ===")
            output.append(f"Sequence length: {len(seq)}")
            output.append("Fragments (bp): " + ", ".join(map(str, sorted(frags, reverse=True))))
            output.append("")

            lanes.append((label, frags))

        output.append("\n=== MERGED ELECTROPHORESIS GEL (all lanes) ===\n")
        output.append(simulate_gel(lanes))

        self.out.delete("1.0", tk.END)
        self.out.insert(tk.END, "\n".join(output))


if __name__ == "__main__":
    root = tk.Tk()
    app = MultiDigestGUI(root)
    root.mainloop()
