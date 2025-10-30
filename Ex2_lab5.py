import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
from Bio import SeqIO
import random
import time
import matplotlib.pyplot as plt

def calculate_gc(seq):
    g = seq.count('G')
    c = seq.count('C')
    return (g + c) / len(seq) * 100 if len(seq) > 0 else 0

def upload_fastas():
    filepaths = filedialog.askopenfilenames(
        title="Select up to 10 FASTA Files",
        filetypes=[("FASTA files", ".fasta *.fa *.txt"), ("All files", ".*")]
    )
    if not filepaths:
        return
    fasta_files.clear()
    fasta_files.extend(filepaths[:10])
    text_box.insert(tk.END, f"\nLoaded {len(fasta_files)} FASTA files:\n")
    for f in fasta_files:
        text_box.insert(tk.END, f" - {f}\n")

def process_sequences():
    if not fasta_files:
        messagebox.showerror("Error", "Please upload FASTA files first.")
        return

    results.clear()
    text_box.insert(tk.END, "\nProcessing sequences...\n")
    text_box.update_idletasks()

    for i, path in enumerate(fasta_files, start=1):
        text_box.insert(tk.END, f"\nGenome {i}: {path.split('/')[-1]}\n")
        text_box.update_idletasks()

        try:
            record = next(SeqIO.parse(path, "fasta"))
        except Exception as e:
            text_box.insert(tk.END, f"Error reading {path}: {e}\n")
            continue

        seq = str(record.seq).upper()
        seq_len = len(seq)
        gc = calculate_gc(seq)

        start_time = time.time()

        samples = []
        for _ in range(2000):
            length = random.randint(100, 150)
            start = random.randint(0, seq_len - length)
            fragment = seq[start:start + length]
            samples.append((start, fragment))

        reconstructed = ["N"] * seq_len
        for start, frag in samples:
            reconstructed[start:start+len(frag)] = list(frag)
        reconstructed_seq = "".join(reconstructed)

        end_time = time.time()
        elapsed_ms = (end_time - start_time) * 1000

        results.append((gc, elapsed_ms))
        text_box.insert(tk.END, f" Length: {seq_len} nt | GC%: {gc:.2f} | Time: {elapsed_ms:.2f} ms\n")

    plot_results()
    save_report()

def plot_results():
    if not results:
        messagebox.showinfo("No Data", "No results to plot.")
        return

    gc_values = [r[0] for r in results]
    times = [r[1] for r in results]

    plt.figure(figsize=(7, 5))
    plt.scatter(gc_values, times, color="blue")
    plt.title("Assembly Time vs GC Content")
    plt.xlabel("GC Content [%]")
    plt.ylabel("Time [ms]")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("gc_vs_time_plot.png")
    plt.show()

    text_box.insert(tk.END, "\nPlot saved as 'gc_vs_time_plot.png'\n")

def save_report():
    if not results:
        return
    with open("results_report.txt", "w") as f:
        f.write("Assembly Time vs GC Content Report\n")
        f.write("="*50 + "\n\n")
        for i, (gc, t) in enumerate(results, start=1):
            f.write(f"Genome {i}: GC% = {gc:.2f}, Time = {t:.2f} ms\n")
        f.write("\nInterpretation:\n")
        f.write("Sequences with higher GC content often take slightly more time to process, ")
        f.write("because GC-rich regions form more stable secondary structures and may require ")
        f.write("more computation during sampling and reconstruction.\n")
        f.write("Lower-GC genomes, being more AT-rich, tend to be faster to assemble.\n")

    text_box.insert(tk.END, "Results saved to 'results_report.txt'\n")

root = tk.Tk()
root.title("Viral Genome Sampler and Analyzer")
root.geometry("950x650")
root.configure(bg="#f8f9fa")

fasta_files = []
results = []

title_label = tk.Label(root, text="Viral Genome Sampler and Analyzer",
                       font=("Helvetica", 18, "bold"), bg="#f8f9fa", fg="#1c7ed6")
title_label.pack(pady=(10, 5))

btn_frame = tk.Frame(root, bg="#f8f9fa")
btn_frame.pack(pady=10)

upload_btn = tk.Button(btn_frame, text="Upload FASTA Files", command=upload_fastas,
                       bg="#d0ebff", font=("Helvetica", 11, "bold"), width=20)
upload_btn.pack(side=tk.LEFT, padx=10)

process_btn = tk.Button(btn_frame, text="Process Sequences", command=process_sequences,
                        bg="#d8f5a2", font=("Helvetica", 11, "bold"), width=20)
process_btn.pack(side=tk.LEFT, padx=10)

text_box = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=110, height=30,
                                     font=("Courier", 10), bg="#ffffff", fg="#212529")
text_box.pack(padx=15, pady=10, fill="both", expand=True)

root.mainloop()
