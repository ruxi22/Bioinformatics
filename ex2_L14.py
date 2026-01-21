import tkinter as tk
from tkinter import ttk, messagebox
import math
import re
import matplotlib.pyplot as plt

# ---------------- Text processing ----------------

def clean_text(text):
    text = text.lower()
    text = re.sub(r"[^a-zăâîșşțţ \n]", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()

def build_alphabet(t1, t2):
    return sorted(set(t1 + t2))

# ---------------- Markov model ----------------

def count_transitions(text, alphabet):
    idx = {c: i for i, c in enumerate(alphabet)}
    n = len(alphabet)
    M = [[0]*n for _ in range(n)]
    for a, b in zip(text[:-1], text[1:]):
        M[idx[a]][idx[b]] += 1
    return M

def normalize_matrix(M, alpha=1.0):
    n = len(M)
    P = [[0.0]*n for _ in range(n)]
    for i in range(n):
        row_sum = sum(M[i]) + alpha*n
        for j in range(n):
            P[i][j] = (M[i][j] + alpha) / row_sum
    return P

def log_likelihood_matrix(P, Q):
    n = len(P)
    B = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            B[i][j] = math.log(P[i][j] / Q[i][j], 2)
    return B

def score_text(text, B, alphabet):
    idx = {c: i for i, c in enumerate(alphabet)}
    score = 0.0
    for a, b in zip(text[:-1], text[1:]):
        score += B[idx[a]][idx[b]]
    return score

def sliding_window(text, B, alphabet, win, step):
    scores, pos = [], []
    for i in range(0, len(text) - win + 1, step):
        w = text[i:i+win]
        scores.append(score_text(w, B, alphabet))
        pos.append(i + win//2)
    return pos, scores

# ---------------- GUI ----------------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Poetry Style Scanner")
        self.geometry("1200x800")
        self.create_widgets()

    def create_widgets(self):
        top = ttk.Frame(self)
        top.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        ttk.Label(top, text="Eminescu (training)").grid(row=0, column=0)
        ttk.Label(top, text="Stănescu (training)").grid(row=0, column=1)
        ttk.Label(top, text="Test (mixture)").grid(row=0, column=2)

        self.txt_A = tk.Text(top, height=15)
        self.txt_B = tk.Text(top, height=15)
        self.txt_T = tk.Text(top, height=15)

        self.txt_A.grid(row=1, column=0, sticky="nsew", padx=5)
        self.txt_B.grid(row=1, column=1, sticky="nsew", padx=5)
        self.txt_T.grid(row=1, column=2, sticky="nsew", padx=5)

        top.columnconfigure(0, weight=1)
        top.columnconfigure(1, weight=1)
        top.columnconfigure(2, weight=1)

        controls = ttk.Frame(self)
        controls.pack(fill=tk.X, padx=10)

        ttk.Label(controls, text="Window size").pack(side=tk.LEFT)
        self.win = ttk.Entry(controls, width=6)
        self.win.insert(0, "350")
        self.win.pack(side=tk.LEFT, padx=5)

        ttk.Label(controls, text="Step").pack(side=tk.LEFT)
        self.step = ttk.Entry(controls, width=6)
        self.step.insert(0, "40")
        self.step.pack(side=tk.LEFT, padx=5)

        ttk.Button(controls, text="Build models", command=self.build_models).pack(side=tk.LEFT, padx=10)
        ttk.Button(controls, text="Scan + chart", command=self.scan).pack(side=tk.LEFT)

        self.output = tk.Text(self, font=("Courier", 10))
        self.output.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    def build_models(self):
        A = clean_text(self.txt_A.get("1.0", tk.END))
        B = clean_text(self.txt_B.get("1.0", tk.END))

        if len(A) < 200 or len(B) < 200:
            messagebox.showerror("Error", "Training texts too short")
            return

        L = min(len(A), len(B))
        A, B = A[:L], B[:L]

        self.alphabet = build_alphabet(A, B)
        self.MA = count_transitions(A, self.alphabet)
        self.MB = count_transitions(B, self.alphabet)

        self.PA = normalize_matrix(self.MA)
        self.PB = normalize_matrix(self.MB)

        self.BETA = log_likelihood_matrix(self.PA, self.PB)

        self.output.delete("1.0", tk.END)
        self.output.insert(tk.END, "Models built successfully\n")
        self.output.insert(tk.END, f"Alphabet size: {len(self.alphabet)}\n")
        self.output.insert(tk.END, f"Training length: {L}\n\n")

    def scan(self):
        T = clean_text(self.txt_T.get("1.0", tk.END))
        if len(T) < 200:
            messagebox.showerror("Error", "Test text too short")
            return

        win = int(self.win.get())
        step = int(self.step.get())

        x, y = sliding_window(T, self.BETA, self.alphabet, win, step)

        plt.figure()
        plt.plot(x, y)
        plt.axhline(0)
        plt.xlabel("Text position (chars)")
        plt.ylabel("Log-likelihood ratio")
        plt.title("Eminescu (above 0) vs Stănescu (below 0)")
        plt.show()



if __name__ == "__main__":
    App().mainloop()
