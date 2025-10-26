import tkinter as tk
from tkinter import messagebox

# RNA codon table
codon_table = {
    # U-start codons
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
    # C-start codons
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    # A-start codons
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",  
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    # G-start codons
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

def rna_to_protein(rna_sequence):
    """
    Converts a coding RNA sequence to an amino acid sequence.
    Stops translation at a stop codon.
    """
    protein_sequence = []
    
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        if len(codon) < 3:
            break
        amino_acid = codon_table.get(codon.upper(), "")
        if amino_acid == "Stop":
            break
        if amino_acid:  
            protein_sequence.append(amino_acid)
    
    return "-".join(protein_sequence)

# GUI setup
def translate_rna():
    rna_seq = rna_entry.get().strip().upper()
    if not rna_seq:
        messagebox.showwarning("Input Error", "Please enter an RNA sequence.")
        return
    protein = rna_to_protein(rna_seq)
    result_label.config(text=f"Amino acid sequence:\n{protein}")

root = tk.Tk()
root.title("RNA to Protein Translator")

tk.Label(root, text="Enter RNA Sequence:").pack(pady=5)
rna_entry = tk.Entry(root, width=50)
rna_entry.pack(pady=5)

translate_button = tk.Button(root, text="Translate", command=translate_rna)
translate_button.pack(pady=10)

result_label = tk.Label(root, text="Amino acid sequence will appear here.", wraplength=400, justify="left")
result_label.pack(pady=10)

root.mainloop()
