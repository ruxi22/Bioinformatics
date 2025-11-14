import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import Counter



def get_kmer_counts(seq, k=5):
    """Count k-mer patterns in a DNA/RNA sequence."""
    kmers = Counter()
    for i in range(len(seq) - k + 1):
        kmers[seq[i:i+k]] += 1
    return kmers


def analyze_and_plot(filename, k=5):
    """Read FASTA file, compute top k-mers, plot chart."""
    try:
        seq_record = next(SeqIO.parse(filename, "fasta"))
    except:
        messagebox.showerror("Error", f"Cannot read FASTA file:\n{filename}")
        return

    seq = str(seq_record.seq).upper()

    counts = get_kmer_counts(seq, k)
    top20 = counts.most_common(20)

    patterns = [x[0] for x in top20]
    values = [x[1] for x in top20]

    plt.figure(figsize=(10, 5))
    plt.bar(patterns, values)
    plt.xticks(rotation=60)
    plt.title(f"Top 20 Patterns (k={k}) for: {os.path.basename(filename)}")
    plt.xlabel("Pattern")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.show()



class InfluenzaApp:

    def __init__(self, master):
        self.master = master
        master.title("FASTA Pattern Analyzer")

        ttk.Label(master, text="FASTA Pattern Analyzer",
                  font=("Arial", 16)).pack(pady=10)

        # k-mer size
        ttk.Label(master, text="k-mer size:").pack()
        self.k_entry = ttk.Entry(master)
        self.k_entry.insert(0, "5")
        self.k_entry.pack()

        # upload button
        ttk.Button(master, text="Upload FASTA Files",
                   command=self.upload_files).pack(pady=5)

        # analyze button
        ttk.Button(master, text="Analyze and Plot",
                   command=self.run_analysis).pack(pady=5)

        # status label
        self.status_label = ttk.Label(master, text="No files loaded.")
        self.status_label.pack(pady=10)

        # list of selected files
        self.files = []

    # --------------------
    # Upload files
    # --------------------
    def upload_files(self):
        selected = filedialog.askopenfilenames(
            title="Select FASTA Files",
            filetypes=[
                ("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa *.frn"),
                ("All files", "*.*")
            ]
        )

        if selected:
            self.files = selected
            self.status_label.config(text=f"{len(self.files)} file(s) loaded.")
        else:
            self.status_label.config(text="No files selected.")

    # --------------------
    # Run analysis AND plot all charts in one window
    # --------------------
    def run_analysis(self):
        if not self.files:
            messagebox.showerror("Error", "Please upload FASTA files first.")
            return

        # read k-mer size
        try:
            k = int(self.k_entry.get())
        except ValueError:
            messagebox.showerror("Invalid Input", "k must be an integer.")
            return

        num_files = len(self.files)
        cols = 3
        rows = (num_files + cols - 1) // cols

        plt.figure(figsize=(16, 9))

        for i, filename in enumerate(self.files):
            seq_record = next(SeqIO.parse(filename, "fasta"))
            seq = str(seq_record.seq).upper()

            kmer_counts = get_kmer_counts(seq, k)
            top20 = kmer_counts.most_common(20)

            patterns = [x[0] for x in top20]
            values = [x[1] for x in top20]

            ax = plt.subplot(rows, cols, i + 1)
            ax.bar(patterns, values)
            ax.set_title(os.path.basename(filename), fontsize=8)
            ax.tick_params(axis='x', rotation=60)
            ax.set_ylabel("Frequency")

        plt.tight_layout()
        plt.show()






if __name__ == "__main__":
    root = tk.Tk()
    app = InfluenzaApp(root)
    root.mainloop()
