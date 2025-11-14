import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext

def load_fasta():
    """Load a FASTA file and return the sequence as a string."""
    filepath = filedialog.askopenfilename(
        filetypes=[("FASTA files", "*.fasta *.fa *.txt"), ("All files", "*.*")]
    )
    if not filepath:
        return

    try:
        with open(filepath, "r") as f:
            lines = f.readlines()
            seq = "".join(line.strip() for line in lines if not line.startswith(">"))
        text_input.delete("1.0", tk.END)
        text_input.insert(tk.END, seq.upper())
    except Exception as e:
        messagebox.showerror("Error", f"Could not load file:\n{e}")

def find_tandem_repeats(seq, min_L=3, max_L=6):
    """Return a list of tandem repeats found in the sequence."""
    results = []
    n = len(seq)

    for L in range(min_L, max_L + 1):
        i = 0
        while i + L <= n:
            unit = seq[i:i+L]
            count = 1
            j = i + L

            while j + L <= n and seq[j:j+L] == unit:
                count += 1
                j += L

            if count >= 3:
                results.append({
                    "start": i + 1,    # convert to 1-based index
                    "unit": unit,
                    "length": L,
                    "count": count,
                    "end": j
                })
                i = j  # skip the entire repeat block
            else:
                i += 1

    return results

def run_detection():
    seq = text_input.get("1.0", tk.END).upper()
    seq = "".join(c for c in seq if c in "ACGTN")

    if len(seq) == 0:
        messagebox.showwarning("Warning", "Sequence is empty or invalid.")
        return

    if any(c not in "ACGTN" for c in seq):
        messagebox.showwarning("Warning", "Sequence contains invalid characters.")
        return

    results = find_tandem_repeats(seq)

    output.delete("1.0", tk.END)

    if not results:
        output.insert(tk.END, "No repeats (3–6 bases) found.")
        return

    for r in results:
        output.insert(
            tk.END,
            f"Start: {r['start']}, Unit: {r['unit']}, "
            f"Length: {r['length']}, Copies: {r['count']}\n"
        )


root = tk.Tk()
root.title("DNA Tandem Repeat Detector (3–6 bases)")
root.geometry("700x600")

# Input Label
tk.Label(root, text="DNA Sequence Input:").pack(anchor="w", padx=10, pady=5)

# Input Text Box
text_input = scrolledtext.ScrolledText(root, height=8, wrap=tk.WORD)
text_input.pack(fill="both", padx=10, pady=5)

# Buttons
frame_buttons = tk.Frame(root)
frame_buttons.pack(pady=10)

tk.Button(frame_buttons, text="Load FASTA", command=load_fasta).grid(row=0, column=0, padx=10)
tk.Button(frame_buttons, text="Run Detection", command=run_detection).grid(row=0, column=1, padx=10)

# Output Label
tk.Label(root, text="Detected Repeats:").pack(anchor="w", padx=10, pady=5)

# Output Window
output = scrolledtext.ScrolledText(root, height=15, wrap=tk.WORD)
output.pack(fill="both", padx=10, pady=5)

root.mainloop()
