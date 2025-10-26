# Codon & Amino Acid Analyzer

This project analyzes the coding regions of genes from genomic FASTA files, converting DNA sequences into amino acid sequences and comparing codon frequencies between genomes. It also provides insights into the most frequently used amino acids and suggests foods low in these amino acids.

---

## Background

- **Codons** are triplets of nucleotides that code for amino acids.  
- Different organisms or viruses may have **different codon usage biases**.  
- Understanding codon frequencies helps in molecular biology, vaccine design, and synthetic gene optimization.  
- Highly used amino acids can be linked to nutritional insights for diet planning.

---

## Features

1. **Gene to Amino Acid Conversion**
   - Converts the coding DNA region of a gene into an amino acid sequence using the standard genetic code.
   
2. **Genome Analysis**
   - Downloads COVID-19 and Influenza genome sequences from NCBI (FASTA files).  
   - Computes codon frequencies for each genome.  
   - Shows **top 10 most frequent codons** in a bar chart for each genome.  
   - Combines results to display **most frequent codons across both genomes**.

3. **Amino Acid Analysis**
   - Identifies the **top 3 amino acids** from each genome.  
   - Generates an **AI prompt** to suggest foods low in the most frequently used amino acids.

4. **Visualization**
   - Charts generated using Matplotlib for codon frequency comparison.  

---

## Technologies

- Python 
- Matplotlib (for charts)  
- Tkinter (GUI interface)  
- Collections (for frequency counting)  

---

##  Usage

1. Load the **FASTA files** for COVID-19 and Influenza.  
2. Run the analysis script.  
3. Outputs:
   - Top 10 codons for each genome  
   - Combined top codons across both genomes  
   - Top 3 amino acids for each genome  
   - AI prompt for nutritional suggestions based on amino acid frequency  

---

##  Notes

- Users can customize input FASTA files.  
- The AI prompt can be used with ChatGPT or other language models for dietary analysis.  
- Codon frequencies are displayed in percentage form and as bar charts for clarity.  
