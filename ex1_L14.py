import math
import tkinter as tk
from tkinter import ttk, messagebox

ALPHABET = ["A", "C", "G", "T"]


S1 = "ATCGATTCGATATCATACACGTAT"      
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"     


def init_matrix():
    return {a: {b: 0 for b in ALPHABET} for a in ALPHABET}

def count_transitions(seq):
    m = init_matrix()
    for x, y in zip(seq[:-1], seq[1:]):
        m[x][y] += 1
    return m

def normalize(counts):
    probs = init_matrix()
    for a in ALPHABET:
        total = sum(counts[a].values())
        for b in ALPHABET:
            probs[a][b] = counts[a][b] / total if total > 0 else 0.0
    return probs

def log_likelihood(tr_p, tr_n):
    beta = init_matrix()
    for a in ALPHABET:
        for b in ALPHABET:
            p = tr_p[a][b]
            q = tr_n[a][b]
            if p == 0 and q == 0:
                beta[a][b] = 0.0
            elif p == 0:
                beta[a][b] = float("-inf")
            else:
                beta[a][b] = math.log(p / q, 2)
    return beta

def score_sequence(seq, beta):
    score = 0.0
    for x, y in zip(seq[:-1], seq[1:]):
        score += beta[x][y]
    return score

def format_matrix(m):
    lines = []
    header = "      " + "   ".join(ALPHABET)
    lines.append(header)
    for a in ALPHABET:
        row = [a]
        for b in ALPHABET:
            v = m[a][b]
            if v == float("-inf"):
                row.append("-inf")
            else:
                row.append(f"{v:.3f}")
        lines.append("{:>3}   ".format(row[0]) + "   ".join(row[1:]))
    return "\n".join(lines)

def run_analysis():
    seq = entry_seq.get().strip().upper()

    if not seq or any(c not in ALPHABET for c in seq):
        messagebox.showerror("Error", "Sequence must contain only A, C, G, T.")
        return

    count_p = count_transitions(S1)
    count_n = count_transitions(S2)

    tr_p = normalize(count_p)
    tr_n = normalize(count_n)

    beta = log_likelihood(tr_p, tr_n)
    llr = score_sequence(seq, beta)

    output.delete("1.0", tk.END)

    output.insert(tk.END, "CpG+ COUNT MATRIX\n")
    output.insert(tk.END, format_matrix(count_p) + "\n\n")

    output.insert(tk.END, "CpG- COUNT MATRIX\n")
    output.insert(tk.END, format_matrix(count_n) + "\n\n")

    output.insert(tk.END, "CpG+ PROBABILITY MATRIX\n")
    output.insert(tk.END, format_matrix(tr_p) + "\n\n")

    output.insert(tk.END, "CpG- PROBABILITY MATRIX\n")
    output.insert(tk.END, format_matrix(tr_n) + "\n\n")

    output.insert(tk.END, "LOG-LIKELIHOOD MATRIX (β)\n")
    output.insert(tk.END, format_matrix(beta) + "\n\n")

    output.insert(tk.END, f"Test sequence: {seq}\n")
    output.insert(tk.END, f"Log-likelihood ratio: {llr}\n")

    if llr > 0:
        output.insert(tk.END, "Classification: CpG ISLAND (+)\n")
    else:
        output.insert(tk.END, "Classification: NON-CpG REGION (-)\n")

root = tk.Tk()
root.title("CpG Island Detector")
root.geometry("820x700")

frame = ttk.Frame(root, padding=10)
frame.pack(fill=tk.BOTH, expand=True)

ttk.Label(frame, text="Input DNA sequence (A, C, G, T):").pack(anchor=tk.W)
entry_seq = ttk.Entry(frame, width=40)
entry_seq.pack(anchor=tk.W, pady=5)
entry_seq.insert(0, "CAGGTTGGAAACGTAA")

ttk.Button(frame, text="Run Analysis", command=run_analysis).pack(anchor=tk.W, pady=5)

output = tk.Text(frame, wrap=tk.NONE, font=("Courier", 10))
output.pack(fill=tk.BOTH, expand=True)

root.mainloop()
