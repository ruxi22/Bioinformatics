import customtkinter as ctk
from tkinter import filedialog
from pathlib import Path
import math
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("green")

# ---------- Core ----------
def clean_dna(seq: str) -> str:
    s = "".join(seq.upper().split()).replace("U", "T")
    return "".join(ch for ch in s if ch in "ACGT")

def read_fasta_sequence(path: Path) -> str:
    parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            parts.append(line)
    return clean_dna("".join(parts))

def tm_wallace(seq: str) -> float:
    c = Counter(seq)
    return 2.0 * (c["A"] + c["T"]) + 4.0 * (c["G"] + c["C"])

def tm_salt_adjusted(seq: str, na_mM: float = 50.0) -> float:
    if not seq:
        return float("nan")
    L = len(seq)
    na_M = max(na_mM, 0.0001) / 1000.0  # mM -> M; avoid log10(0)
    gc = 100.0 * (seq.count("G") + seq.count("C")) / L
    return 81.5 + 16.6 * math.log10(na_M) + 0.41 * gc - (600.0 / L)

def sliding_tm_vectors(seq: str, k: int, na_mM: float):
    X, Pw, Ps = [], [], []
    for i in range(len(seq) - k + 1):
        win = seq[i:i+k]
        X.append(i + 1)
        Pw.append(tm_wallace(win))
        Ps.append(tm_salt_adjusted(win, na_mM))
    return np.array(X), np.array(Pw, float), np.array(Ps, float)

def mask_to_segments(X: np.ndarray, mask: np.ndarray):
    segs = []
    if len(X) == 0:
        return segs
    in_seg = False
    start = None
    for i, flag in enumerate(mask):
        if flag and not in_seg:
            in_seg = True
            start = X[i]
        elif not flag and in_seg:
            in_seg = False
            stop = X[i-1]
            segs.append((int(start), int(stop - start + 1)))
    if in_seg:
        stop = X[-1]
        segs.append((int(start), int(stop - start + 1)))
    return segs

# ---------- Plot ----------
def plot_board_style_gui(X, Pw, Ps, k):
    w_cut = float(np.mean(Pw))
    s_cut = float(np.mean(Ps))

    fig = plt.figure(figsize=(10, 6))
    gs = fig.add_gridspec(nrows=3, height_ratios=[4.0, 0.6, 0.6], hspace=0.07)
    ax = fig.add_subplot(gs[0])
    ax_w = fig.add_subplot(gs[1], sharex=ax)
    ax_s = fig.add_subplot(gs[2], sharex=ax)

    line_w, = ax.plot(X, Pw, linewidth=1.8, label="Wallace Tm")
    line_s, = ax.plot(X, Ps, linewidth=1.8, label="Salt-adjusted Tm")

    ax.axhline(w_cut, linestyle="--", color=line_w.get_color(), alpha=0.9, label=f"Wallace cutoff = {w_cut:.2f}°C")
    ax.axhline(s_cut, linestyle="--", color=line_s.get_color(), alpha=0.9, label=f"Salt cutoff = {s_cut:.2f}°C")

    ax.set_ylabel("Melting temperature (°C)")
    ax.set_title(f"Sliding-Window Tm (k={k}) with Cutoff Lines and Activity Bars")
    ax.grid(alpha=0.35)
    ax.legend(ncols=2, fontsize=9)

    mask_w = Pw >= w_cut
    segs_w = mask_to_segments(X, mask_w)
    ax_w.set_yticks([]); ax_w.set_ylim(0, 1); ax_w.set_ylabel("W", rotation=0, labelpad=10, va="center")
    for start, length in segs_w:
        ax_w.broken_barh([(start, length)], (0.25, 0.5), facecolor=line_w.get_color(), alpha=0.85)
    ax_w.grid(alpha=0.25, axis="x")

    mask_s = Ps >= s_cut
    segs_s = mask_to_segments(X, mask_s)
    ax_s.set_yticks([]); ax_s.set_ylim(0, 1); ax_s.set_ylabel("S", rotation=0, labelpad=10, va="center")
    for start, length in segs_s:
        ax_s.broken_barh([(start, length)], (0.25, 0.5), facecolor=line_s.get_color(), alpha=0.85)
    ax_s.grid(alpha=0.25, axis="x")

    ax_s.set_xlabel(f"Position in sequence (start of {k}-nt window)")
    plt.setp(ax_w.get_xticklabels(), visible=False)
    fig.tight_layout()
    return fig

# ---------- GUI ----------
def browse_file():
    file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa")])
    fasta_entry.delete(0, "end")
    fasta_entry.insert(0, file_path)

def plot_gui():
    fasta_path = fasta_entry.get().strip()
    try:
        window = int(window_entry.get())
        na_mM = float(na_entry.get())
    except:
        return
    seq = read_fasta_sequence(Path(fasta_path))
    X, Pw, Ps = sliding_tm_vectors(seq, window, na_mM)
    fig = plot_board_style_gui(X, Pw, Ps, window)
    for widget in plot_frame.winfo_children():
        widget.destroy()
    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill="both", expand=True)

# ---------- Tkinter GUI ----------
app = ctk.CTk()
app.geometry("900x700")
app.title("Sliding Window Tm Calculator with Cutoff Bars")

header = ctk.CTkLabel(app, text="DNA Sliding Window Tm Calculator", font=ctk.CTkFont(size=24, weight="bold"), text_color="#FFD700")
header.pack(pady=10)

input_frame = ctk.CTkFrame(app, corner_radius=20, fg_color="#1E1E2F")
input_frame.pack(padx=20, pady=10, fill="x")

fasta_entry = ctk.CTkEntry(input_frame, width=400, placeholder_text="FASTA file path")
fasta_entry.grid(row=0, column=0, padx=10, pady=10)
browse_btn = ctk.CTkButton(input_frame, text="Browse", command=browse_file)
browse_btn.grid(row=0, column=1, padx=10, pady=10)

window_entry = ctk.CTkEntry(input_frame, width=100, placeholder_text="Window size")
window_entry.insert(0, "9")
window_entry.grid(row=1, column=0, padx=10, pady=5)

na_entry = ctk.CTkEntry(input_frame, width=100, placeholder_text="[Na+] mM")
na_entry.insert(0, "50")
na_entry.grid(row=1, column=1, padx=10, pady=5)

plot_btn = ctk.CTkButton(input_frame, text="Plot Tm", command=plot_gui, width=200, height=40, fg_color="#00BFFF", hover_color="#1E90FF")
plot_btn.grid(row=2, column=0, columnspan=2, padx=10, pady=10)

plot_frame = ctk.CTkFrame(app, corner_radius=20, fg_color="#2E2E3F")
plot_frame.pack(padx=20, pady=10, fill="both", expand=True)

app.mainloop()
