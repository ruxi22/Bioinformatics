import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
from Bio import SeqIO
import random

def upload_fasta():
    filepath = filedialog.askopenfilename(
        title="Select FASTA file",
        filetypes=[("FASTA files", ".fasta *.fa *.txt"), ("All files", ".*")]
    )
    if not filepath:
        return
    fasta_path.set(filepath)
    status_label.config(text=f"Loaded: {filepath.split('/')[-1]}", fg="#2b8a3e")
    text_box.insert(tk.END, f"\nLoaded FASTA file: {filepath}\n")

def process_sequence():
    path = fasta_path.get()
    if not path:
        messagebox.showerror("Error", "Please upload a FASTA file first.")
        return

    status_label.config(text="Processing sequence...", fg="#f08c00")
    progress_bar.start(10)
    root.update_idletasks()

    try:
        record = next(SeqIO.parse(path, "fasta"))
    except Exception as e:
        messagebox.showerror("Error", f"Failed to parse FASTA: {e}")
        progress_bar.stop()
        status_label.config(text="Error reading FASTA", fg="red")
        return

    seq = str(record.seq).upper()
    seq_len = len(seq)
    text_box.insert(tk.END, f"\nOriginal sequence length: {seq_len}\n")
    text_box.insert(tk.END, f"Sequence (first 300 nt): {seq[:300]}\n\n")

    if seq_len < 1000 or seq_len > 3000:
        text_box.insert(tk.END, "Warning: Sequence length not in 1000–3000 range\n\n")

    samples = []
    for _ in range(2000):
        length = random.randint(100, 150)
        start = random.randint(0, seq_len - length)
        fragment = seq[start:start + length]
        samples.append((start, fragment))

    text_box.insert(tk.END, f"Generated {len(samples)} random samples (100–150 nt each)\n")
    text_box.insert(tk.END, "Example samples (with positions):\n")
    for i, (pos, frag) in enumerate(samples[:5]):
        text_box.insert(tk.END, f"  Sample {i+1}: start={pos}, seq={frag[:80]}...\n")
    text_box.insert(tk.END, "\n")

    text_box.insert(tk.END, "Rebuilding original sequence...\n")
    samples.sort(key=lambda x: x[0])

    reconstructed = ["N"] * seq_len
    for start, frag in samples:
        reconstructed[start:start+len(frag)] = list(frag)
    reconstructed_seq = "".join(reconstructed)

    text_box.insert(tk.END, "Reconstruction complete.\n")
    text_box.insert(tk.END, f"Reconstructed sequence length: {len(reconstructed_seq)}\n")
    text_box.insert(tk.END, f"Reconstructed (first 300 nt): {reconstructed_seq[:300]}\n\n")

    if reconstructed_seq == seq:
        text_box.insert(tk.END, "The reconstructed sequence is IDENTICAL to the original.\n", "success")
    else:
        match_count = sum(1 for a, b in zip(seq, reconstructed_seq) if a == b)
        similarity = match_count / len(seq) * 100
        text_box.insert(tk.END, "The sequences differ.\n", "error")
        text_box.insert(tk.END, f"Similarity: {similarity:.2f}%\n")

    text_box.insert(tk.END, "\n" + "="*70 + "\n")
    progress_bar.stop()
    status_label.config(text="Processing complete.", fg="#2b8a3e")

def clear_output():
    text_box.delete(1.0, tk.END)
    status_label.config(text="Output cleared.", fg="#495057")

root = tk.Tk()
root.title("DNA Sequence Rebuilder")
root.geometry("1000x650")
root.configure(bg="#e9ecef")

fasta_path = tk.StringVar()

title_frame = tk.Frame(root, bg="#343a40", height=60)
title_frame.pack(fill="x", side="top")

title_label = tk.Label(title_frame, text="DNA Sequence Rebuilder", font=("Helvetica", 18, "bold"), bg="#343a40", fg="white")
title_label.pack(pady=10)

main_frame = tk.Frame(root, bg="#e9ecef")
main_frame.pack(fill="both", expand=True)

side_frame = tk.Frame(main_frame, width=220, bg="#dee2e6")
side_frame.pack(side="left", fill="y")

button_style = {"font": ("Helvetica", 11, "bold"), "width": 18, "height": 2, "bd": 0}

upload_btn = tk.Button(side_frame, text="Upload FASTA", command=upload_fasta, bg="#adb5bd", **button_style)
upload_btn.pack(pady=(40, 15))

process_btn = tk.Button(side_frame, text="Process Sequence", command=process_sequence, bg="#a3cfbb", **button_style)
process_btn.pack(pady=10)

clear_btn = tk.Button(side_frame, text="Clear Output", command=clear_output, bg="#ffe066", **button_style)
clear_btn.pack(pady=10)

exit_btn = tk.Button(side_frame, text="Exit", command=root.destroy, bg="#ff8787", **button_style)
exit_btn.pack(pady=(60, 0))

output_frame = tk.Frame(main_frame, bg="#f8f9fa")
output_frame.pack(side="right", fill="both", expand=True, padx=10, pady=10)

text_box = scrolledtext.ScrolledText(output_frame, wrap=tk.WORD, width=90, height=35, font=("Courier", 10), bg="white", fg="black")
text_box.pack(fill="both", expand=True, padx=5, pady=5)

text_box.tag_config("success", foreground="green", font=("Courier", 10, "bold"))
text_box.tag_config("error", foreground="red", font=("Courier", 10, "bold"))

status_frame = tk.Frame(root, bg="#dee2e6", height=30)
status_frame.pack(fill="x", side="bottom")

status_label = tk.Label(status_frame, text="No file loaded.", bg="#dee2e6", fg="#495057", font=("Helvetica", 10))
status_label.pack(side="left", padx=10)

progress_bar = ttk.Progressbar(status_frame, mode="indeterminate", length=200)
progress_bar.pack(side="right", padx=10, pady=5)

root.mainloop()
