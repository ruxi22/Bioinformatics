import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
from collections import Counter

def open_fasta():
    filepath = filedialog.askopenfilename(filetypes=[("FASTA files","*.fasta *.fa")])
    if not filepath:
        return
    seq = ""
    with open(filepath) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip().upper()
    if not seq:
        messagebox.showerror("Error","No sequence found in the file!")
        return
    alphabet = sorted(set(seq))
    percentages = {k: round(v/len(seq)*100,2) for k,v in Counter(seq).items()}
 
    result_box.config(state='normal')
    result_box.delete('1.0', tk.END)
    result_box.insert(tk.END, f"FASTA file: {filepath}\n\n")
    result_box.insert(tk.END, f"Alphabet: {', '.join(alphabet)}\n\n")
    result_box.insert(tk.END, "Percentages:\n")
    color_map = {'A':"#d01c1c",'C':"#1a0090",'G':"#087208",'T':"#BDC916"}
    for k in alphabet:
        color = color_map.get(k,'#ffffff')
        result_box.insert(tk.END, f"{k}: {percentages[k]}%\n", k)
        result_box.tag_config(k, foreground=color, font=("Courier",11,"bold"))
    result_box.config(state='disabled')

# Main window
root = tk.Tk()
root.title("DNA Sequence Analyzer")
root.geometry("600x500")
root.configure(bg="#e0f7fa")

# Header
header = tk.Label(root, text="DNA Sequence Analyzer", font=("Helvetica", 18, "bold"), bg="#e0f7fa", fg="#00796b")
header.pack(pady=10)

# Button
btn = tk.Button(root, text="Select FASTA File", font=("Helvetica", 12, "bold"), bg="#00796b", fg="white", command=open_fasta)
btn.pack(pady=5)

# Scrollable result box
result_box = scrolledtext.ScrolledText(root, width=70, height=20, font=("Courier", 11))
result_box.pack(pady=10)
result_box.config(state='disabled')

root.mainloop()
