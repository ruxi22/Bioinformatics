import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
import math

# ============================================================
#  RESTRICTION ENZYME DEFINITIONS
# ============================================================

ENZYMES = {
    "EcoRI": {"site": "GAATTC", "offset": 1},
    "BamHI": {"site": "GGATCC", "offset": 1},
    "HindIII": {"site": "AAGCTT", "offset": 1},
    "TaqI": {"site": "TCGA", "offset": 1},
    "HaeIII": {"site": "GGCC", "offset": 2},
}

# ============================================================
#  DIGEST ENGINE
# ============================================================

def clean_sequence(seq):
    """Remove FASTA header and whitespace."""
    seq = seq.upper()
    seq = "".join(c for c in seq if c in "ACGT")
    return seq


def find_cuts(seq, site, offset):
    """Find 0-based cleavage positions."""
    cuts = []
    pos = seq.find(site)
    while pos != -1:
        cuts.append(pos + offset)
        pos = seq.find(site, pos + 1)
    return cuts


def compute_fragments(length, cuts):
    """Return fragment lengths for a linear DNA."""
    if not cuts:
        return [length]
    cuts = [0] + sorted(cuts) + [length]
    return [cuts[i+1] - cuts[i] for i in range(len(cuts)-1)]


def simulate_gel(lanes, height=20):
    """
    lanes = [(name, [fragment sizes])]
    Returns ASCII gel text.
    """
    txt = "Simulated Agarose Gel:\n(top = wells, bottom = migration)\n\n"

    all_sizes = [s for _, frags in lanes for s in frags]
    if not all_sizes:
        return "No fragments."

    max_s = max(all_sizes)
    min_s = min(all_sizes)

    lane_rows = []
    for name, frags in lanes:
        rows = set()
        for s in frags:
            y = (math.log10(max_s) - math.log10(s)) / (
                math.log10(max_s) - math.log10(min_s) + 1e-6
            )
            row = int(y * (height - 1))
            rows.add(row)
        lane_rows.append((name, rows))

    for r in range(height):
        line = ""
        for _, rows in lane_rows:
            line += "|" + ("*" if r in rows else " ") + " "
        txt += line + "\n"

    txt += "\n"
    for i, (name, _) in enumerate(lane_rows):
        txt += f"Lane {i+1}: {name}\n"

    return txt

# ============================================================
#  GUI INTERFACE
# ============================================================

class RestrictionGUI:
    def __init__(self, root):
        self.root = root
        root.title("Restriction Enzyme Digest Simulator")

        # DNA input box
        tk.Label(root, text="DNA Sequence (FASTA or raw):").grid(row=0, column=0, sticky="w")
        self.seq_box = scrolledtext.ScrolledText(root, width=80, height=10)
        self.seq_box.grid(row=1, column=0, columnspan=3, padx=5, pady=5)

        # Enzyme selection
        tk.Label(root, text="Select Enzymes:").grid(row=2, column=0, sticky="w")
        self.enzyme_vars = {}

        enzyme_frame = tk.Frame(root)
        enzyme_frame.grid(row=3, column=0, sticky="w")

        for i, name in enumerate(ENZYMES.keys()):
            var = tk.BooleanVar()
            tk.Checkbutton(enzyme_frame, text=name, variable=var).grid(row=0, column=i, sticky="w")
            self.enzyme_vars[name] = var

        # Run button
        self.run_button = tk.Button(root, text="Run Digest", command=self.run_digest)
        self.run_button.grid(row=4, column=0, pady=10, sticky="w")

        # Output box
        tk.Label(root, text="Output:").grid(row=5, column=0, sticky="w")
        self.output_box = scrolledtext.ScrolledText(root, width=100, height=30)
        self.output_box.grid(row=6, column=0, columnspan=3, padx=5, pady=5)

    # =====================================================
    # Perform digest
    # =====================================================
    def run_digest(self):
        dna = clean_sequence(self.seq_box.get("1.0", tk.END))
        if len(dna) < 4:
            messagebox.showerror("Error", "DNA sequence too short or invalid.")
            return

        selected = {name: ENZYMES[name] for name, v in self.enzyme_vars.items() if v.get()}
        if not selected:
            messagebox.showerror("Error", "Select at least one enzyme.")
            return

        out = []
        lanes = []

        out.append(f"Sequence length: {len(dna)} bp\n")

        # Single enzyme digests
        for name, edef in selected.items():
            cuts = find_cuts(dna, edef["site"], edef["offset"])
            frags = compute_fragments(len(dna), cuts)

            out.append(f"=== {name} digest ===")
            out.append(f"Recognizes: {edef['site']}")
            out.append(f"Cleavages: {len(cuts)}")
            out.append("Cut positions (1-based): " + (", ".join(str(c+1) for c in cuts) if cuts else "None"))
            out.append("Fragments (bp): " + ", ".join(map(str, sorted(frags, reverse=True))))
            out.append("")
            lanes.append((name, frags))

        # Multi-enzyme digest
        if len(selected) > 1:
            out.append("=== Combined digest ===")
            all_cuts = []
            for name, edef in selected.items():
                all_cuts += find_cuts(dna, edef["site"], edef["offset"])
            all_cuts = sorted(set(all_cuts))
            frags = compute_fragments(len(dna), all_cuts)
            combo_name = "+".join(selected.keys())

            out.append(f"Digest: {combo_name}")
            out.append(f"Cleavages: {len(all_cuts)}")
            out.append("Cut positions (1-based): " + (", ".join(str(c+1) for c in all_cuts) if all_cuts else "None"))
            out.append("Fragments (bp): " + ", ".join(map(str, sorted(frags, reverse=True))))
            out.append("")
            lanes.append((combo_name, frags))

        # Show gel
        out.append(simulate_gel(lanes))

        self.output_box.delete("1.0", tk.END)
        self.output_box.insert(tk.END, "\n".join(out))


# ============================================================
#  RUN APP
# ============================================================

if __name__ == "__main__":
    root = tk.Tk()
    gui = RestrictionGUI(root)
    root.mainloop()
