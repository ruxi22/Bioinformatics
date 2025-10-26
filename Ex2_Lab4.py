import tkinter as tk
from tkinter import filedialog, messagebox
from collections import Counter
import matplotlib.pyplot as plt

# Codon table
codon_table = {
    "UUU": "Phe","UUC": "Phe","UUA": "Leu","UUG": "Leu",
    "UCU": "Ser","UCC": "Ser","UCA": "Ser","UCG": "Ser",
    "UAU": "Tyr","UAC": "Tyr","UAA": "Stop","UAG": "Stop",
    "UGU": "Cys","UGC": "Cys","UGA": "Stop","UGG": "Trp",
    "CUU": "Leu","CUC": "Leu","CUA": "Leu","CUG": "Leu",
    "CCU": "Pro","CCC": "Pro","CCA": "Pro","CCG": "Pro",
    "CAU": "His","CAC": "His","CAA": "Gln","CAG": "Gln",
    "CGU": "Arg","CGC": "Arg","CGA": "Arg","CGG": "Arg",
    "AUU": "Ile","AUC": "Ile","AUA": "Ile","AUG": "Met",
    "ACU": "Thr","ACC": "Thr","ACA": "Thr","ACG": "Thr",
    "AAU": "Asn","AAC": "Asn","AAA": "Lys","AAG": "Lys",
    "AGU": "Ser","AGC": "Ser","AGA": "Arg","AGG": "Arg",
    "GUU": "Val","GUC": "Val","GUA": "Val","GUG": "Val",
    "GCU": "Ala","GCC": "Ala","GCA": "Ala","GCG": "Ala",
    "GAU": "Asp","GAC": "Asp","GAA": "Glu","GAG": "Glu",
    "GGU": "Gly","GGC": "Gly","GGA": "Gly","GGG": "Gly",
}

# Mapping of low-amino-acid foods
low_amino_foods = {
    "Leu": ["Apples", "Oranges", "Berries", "Lettuce", "Rice", "White bread"],
    "Ala": ["Fruits", "Leafy greens", "Rice", "Corn", "Wheat products"],
    "Gly": ["Fruits", "Vegetables", "Rice", "Refined cereals"],
    "Ser": ["Fruits", "Vegetables like lettuce, cucumber, zucchini", "Rice", "White bread"],
    "Thr": ["Fruits", "Rice", "Refined grains", "Certain vegetables like lettuce, cucumber"]
}



def read_fasta(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return seq.upper()

def dna_to_rna(dna_seq):
    return dna_seq.replace("T", "U")

def codon_count(rna_seq):
    codons = [rna_seq[i:i+3] for i in range(0,len(rna_seq)-2,3)]
    codons = [c for c in codons if c in codon_table and codon_table[c] != "Stop"]
    return Counter(codons)

def amino_acid_count(codon_counter):
    aa_counter = Counter()
    for codon, count in codon_counter.items():
        aa = codon_table[codon]
        aa_counter[aa] += count
    return aa_counter

def plot_top_codons(top_codons, title):
    codons, counts = zip(*top_codons)
    plt.figure(figsize=(10,5))
    plt.bar(codons, counts, color='skyblue')
    plt.title(title)
    plt.xlabel("Codons")
    plt.ylabel("Frequency")
    plt.show()

def process_genomes(covid_file, influenza_file):
    covid_seq = dna_to_rna(read_fasta(covid_file))
    influenza_seq = dna_to_rna(read_fasta(influenza_file))

    covid_codons = codon_count(covid_seq)
    influenza_codons = codon_count(influenza_seq)
    combined_codons = covid_codons + influenza_codons

    top10_covid = covid_codons.most_common(10)
    top10_influenza = influenza_codons.most_common(10)
    top10_combined = combined_codons.most_common(10)

    covid_aa = amino_acid_count(covid_codons)
    influenza_aa = amino_acid_count(influenza_codons)

    top3_covid_aa = covid_aa.most_common(3)
    top3_influenza_aa = influenza_aa.most_common(3)

    ai_prompt = f"The top three amino acids most frequently used in the SARS-CoV-2 genome are {', '.join([aa for aa,_ in top3_covid_aa])}. Suggest foods that are low in these amino acids."

    # Generate food suggestions
    food_suggestions = {}
    for aa, _ in top3_covid_aa:
        food_suggestions[aa] = low_amino_foods.get(aa, ["No data"])

    return (top10_covid, top10_influenza, top10_combined, top3_covid_aa, top3_influenza_aa, ai_prompt, food_suggestions)

# GUI functions
def browse_covid():
    path = filedialog.askopenfilename(filetypes=[("FASTA files","*.fasta")])
    covid_entry.delete(0, tk.END)
    covid_entry.insert(0, path)

def browse_influenza():
    path = filedialog.askopenfilename(filetypes=[("FASTA files","*.fasta")])
    influenza_entry.delete(0, tk.END)
    influenza_entry.insert(0, path)

def run_analysis():
    covid_file = covid_entry.get().strip()
    influenza_file = influenza_entry.get().strip()
    if not covid_file or not influenza_file:
        messagebox.showwarning("Input Error", "Please select both FASTA files.")
        return
    results = process_genomes(covid_file, influenza_file)
    top10_covid, top10_influenza, top10_combined, top3_covid_aa, top3_influenza_aa, ai_prompt, food_suggestions = results

    # Plot charts
    plot_top_codons(top10_covid, "Top 10 COVID-19 Codons")
    plot_top_codons(top10_influenza, "Top 10 Influenza Codons")
    plot_top_codons(top10_combined, "Top 10 Combined Codons")

    # Display results
    result_text = f"Top 3 COVID-19 Amino Acids: {top3_covid_aa}\n"
    result_text += f"Top 3 Influenza Amino Acids: {top3_influenza_aa}\n\n"
    result_text += f"AI Prompt:\n{ai_prompt}\n\n"
    result_text += "Food suggestions for COVID-19 top amino acids:\n"
    for aa, foods in food_suggestions.items():
        result_text += f"{aa}: {', '.join(foods)}\n"

    result_box.config(state="normal")
    result_box.delete("1.0", tk.END)
    result_box.insert(tk.END, result_text)
    result_box.config(state="disabled")

# Tkinter GUI layout
root = tk.Tk()
root.title("Genome Codon Frequency Analyzer")

tk.Label(root, text="COVID-19 FASTA:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
covid_entry = tk.Entry(root, width=50)
covid_entry.grid(row=0, column=1, padx=5, pady=5)
tk.Button(root, text="Browse", command=browse_covid).grid(row=0, column=2, padx=5, pady=5)

tk.Label(root, text="Influenza FASTA:").grid(row=1, column=0, sticky="e", padx=5, pady=5)
influenza_entry = tk.Entry(root, width=50)
influenza_entry.grid(row=1, column=1, padx=5, pady=5)
tk.Button(root, text="Browse", command=browse_influenza).grid(row=1, column=2, padx=5, pady=5)

tk.Button(root, text="Run Analysis", command=run_analysis, bg="lightgreen").grid(row=2, column=0, columnspan=3, pady=10)

result_box = tk.Text(root, width=80, height=15, state="disabled")
result_box.grid(row=3, column=0, columnspan=3, padx=10, pady=10)

root.mainloop()
