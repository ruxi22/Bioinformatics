import customtkinter as ctk
from tkinter import filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import math
import os

ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("green")

def read_fasta(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return seq.upper()

def tm_simple(dna):
    return 4*(dna.count('G')+dna.count('C')) + 2*(dna.count('A')+dna.count('T'))

def tm_salt_adjusted(dna, na_conc=0.05):
    length = len(dna)
    gc_percent = ((dna.count('G')+dna.count('C')) / length) * 100
    return 81.5 + 16.6*math.log10(na_conc) + 0.41*gc_percent - 600/length

def sliding_window_tm(seq, window_size, step, na_conc):
    positions, tm_simple_values, tm_salt_values = [], [], []
    for i in range(0, len(seq)-window_size+1, step):
        window = seq[i:i+window_size]
        positions.append(i+1)
        tm_simple_values.append(tm_simple(window))
        tm_salt_values.append(tm_salt_adjusted(window, na_conc))
    return positions, tm_simple_values, tm_salt_values

def browse_file():
    file_path = filedialog.askopenfilename(filetypes=[("FASTA files","*.fasta *.fa")])
    fasta_entry.delete(0, "end")
    fasta_entry.insert(0, file_path)

def plot_tm():
    file_path = fasta_entry.get().strip()
    if not os.path.isfile(file_path):
        result_label.configure(text="Invalid file path!", text_color="red")
        return
    try:
        window_size = int(window_entry.get())
        step_size = int(step_entry.get())
        na_conc = float(na_entry.get())
    except:
        result_label.configure(text="Enter valid numbers!", text_color="red")
        return

    seq = read_fasta(file_path)
    positions, tm_simple_values, tm_salt_values = sliding_window_tm(seq, window_size, step_size, na_conc)

    fig = plt.Figure(figsize=(7,4))
    ax = fig.add_subplot(111)
    ax.plot(positions, tm_simple_values, label="Tm Simple", linewidth=2)
    ax.plot(positions, tm_salt_values, label="Tm Salt-adjusted", linewidth=2)
    ax.set_title("Sliding Window Tm Profile")
    ax.set_xlabel("Position in Sequence")
    ax.set_ylabel("Melting Temperature (°C)")
    ax.legend()
    ax.grid(True)

    for widget in plot_frame.winfo_children():
        widget.destroy()
    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill="both", expand=True)

app = ctk.CTk()
app.title("Sliding Window DNA Tm Calculator")
app.geometry("800x600")

header = ctk.CTkLabel(app, text="DNA Sliding Window Tm Calculator", font=ctk.CTkFont(size=24, weight="bold"), text_color="#FFD700")
header.pack(pady=10)

input_frame = ctk.CTkFrame(app, corner_radius=20, fg_color="#1E1E2F")
input_frame.pack(pady=10, padx=20, fill="x")

fasta_entry = ctk.CTkEntry(input_frame, width=400, placeholder_text="FASTA file path")
fasta_entry.grid(row=0, column=0, padx=10, pady=10)
browse_button = ctk.CTkButton(input_frame, text="Browse", command=browse_file)
browse_button.grid(row=0, column=1, padx=10, pady=10)

window_entry = ctk.CTkEntry(input_frame, width=100, placeholder_text="Window size")
window_entry.insert(0,"9")
window_entry.grid(row=1, column=0, padx=10, pady=5)

step_entry = ctk.CTkEntry(input_frame, width=100, placeholder_text="Step size")
step_entry.insert(0,"1")
step_entry.grid(row=1, column=1, padx=10, pady=5)

na_entry = ctk.CTkEntry(input_frame, width=100, placeholder_text="Na+ conc (M)")
na_entry.insert(0,"0.05")
na_entry.grid(row=2, column=0, padx=10, pady=5)

calc_button = ctk.CTkButton(input_frame, text="Plot Tm", command=plot_tm, width=200, height=40, fg_color="#00BFFF", hover_color="#1E90FF")
calc_button.grid(row=2, column=1, padx=10, pady=5)

result_label = ctk.CTkLabel(app, text="", font=ctk.CTkFont(size=14))
result_label.pack(pady=5)

plot_frame = ctk.CTkFrame(app, corner_radius=20, fg_color="#2E2E3F")
plot_frame.pack(padx=20, pady=10, fill="both", expand=True)

app.mainloop()
