# use  a dna sequence of 50 letters in order to calculate the transition matrix that represents this sequence and store it in a a json file
import tkinter as tk
from tkinter import messagebox, filedialog
import json

NUCLEOTIDES = ["A", "C", "G", "T"]



def compute_transition_matrix(sequence):
    sequence = sequence.upper()

    counts = {
        a: {b: 0 for b in NUCLEOTIDES}
        for a in NUCLEOTIDES
    }

    for i in range(len(sequence) - 1):
        a = sequence[i]
        b = sequence[i + 1]
        counts[a][b] += 1

    matrix = {}
    for a in NUCLEOTIDES:
        total = sum(counts[a].values())
        if total == 0:
            matrix[a] = {b: 0.0 for b in NUCLEOTIDES}
        else:
            matrix[a] = {
                b: counts[a][b] / total
                for b in NUCLEOTIDES
            }

    return matrix


def calculate():
    seq = seq_input.get().strip().upper()

    if not seq:
        messagebox.showerror("Error", "DNA sequence is empty")
        return

    if any(c not in NUCLEOTIDES for c in seq):
        messagebox.showerror(
            "Error",
            "Sequence must contain only A, C, G, T"
        )
        return

    if len(seq) < 2:
        messagebox.showerror(
            "Error",
            "Sequence must contain at least 2 letters"
        )
        return

    matrix = compute_transition_matrix(seq)

    output.delete("1.0", tk.END)
    for a in NUCLEOTIDES:
        row = "  ".join(f"{matrix[a][b]:.4f}" for b in NUCLEOTIDES)
        output.insert(tk.END, f"{a} → {row}\n")

    output.matrix = matrix  # stash for saving

def save_json():
    if not hasattr(output, "matrix"):
        messagebox.showerror("Error", "Nothing to save")
        return

    path = filedialog.asksaveasfilename(
        defaultextension=".json",
        filetypes=[("JSON files", "*.json")]
    )

    if not path:
        return

    with open(path, "w") as f:
        json.dump(output.matrix, f, indent=4)

    messagebox.showinfo("Saved", "Transition matrix saved successfully")


root = tk.Tk()
root.title("DNA Transition Matrix")
root.resizable(False, False)

main = tk.Frame(root, padx=12, pady=12)
main.pack()

tk.Label(
    main,
    text="DNA sequence (A, C, G, T):"
).grid(row=0, column=0, sticky="w")

seq_input = tk.Entry(main, width=60)
seq_input.grid(row=1, column=0, columnspan=2, pady=(0, 8))

seq_input.insert(
    0,
    "ACGTGCTAGCTAGCTAGCGTACGTAGCTAGCTAGCGTAGCTAGCTAGC"
)

tk.Button(
    main,
    text="Calculate transition matrix",
    command=calculate,
    width=30
).grid(row=2, column=0, pady=6, sticky="w")

tk.Button(
    main,
    text="Save as JSON",
    command=save_json,
    width=20
).grid(row=2, column=1, pady=6, sticky="e")

tk.Label(main, text="Transition matrix:").grid(
    row=3, column=0, columnspan=2, sticky="w", pady=(8, 0)
)

output = tk.Text(
    main,
    height=8,
    width=60,
    font=("Courier", 10)
)
output.grid(row=4, column=0, columnspan=2)

root.mainloop()
