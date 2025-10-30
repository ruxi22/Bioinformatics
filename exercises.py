import tkinter as tk
from tkinter import scrolledtext
from itertools import product

nucleotides = ['A', 'T', 'G', 'C']
sequence = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
di_combos = []
tri_combos = []
di_percentages = {}
tri_percentages = {}
exercise_index = 0

def exercise1():
    global di_combos, tri_combos
    di_combos = [''.join(p) for p in product(nucleotides, repeat=2)]
    tri_combos = [''.join(p) for p in product(nucleotides, repeat=3)]
    output_text.config(state='normal')
    output_text.delete(1.0, tk.END)
    output_text.insert(tk.END, "Exercise 1 Solution:\n"
                               f"Dinucleotides: {', '.join(di_combos)}\n"
                               f"Trinucleotides: {', '.join(tri_combos)}\n")
    output_text.config(state='disabled')
    next_button.pack()

def exercise2():
    output_text.config(state='normal')
    output_text.delete(1.0, tk.END)
    output_text.insert(tk.END, "Exercise 2 Solution:\n")
    output_text.insert(tk.END, "Formula for percentage of each dinucleotide/trinucleotide:\n")
    output_text.insert(tk.END, "Percentage = (count of combination in sequence / total possible positions) * 100\n")
    global di_percentages, tri_percentages
    di_percentages = {combo: sequence.count(combo)/(len(sequence)-1)*100 for combo in di_combos}
    tri_percentages = {combo: sequence.count(combo)/(len(sequence)-2)*100 for combo in tri_combos}
    output_text.insert(tk.END, "Percentages are calculated and will be printed in Exercise 3.\n")
    output_text.config(state='disabled')
    next_button.pack()

def exercise3():
    output_text.config(state='normal')
    output_text.delete(1.0, tk.END)
    output_text.insert(tk.END, "Exercise 3 Solution:\nDinucleotide percentages:\n")
    for combo, perc in di_percentages.items():
        output_text.insert(tk.END, f"{combo}: {perc:.2f}%\n")
    output_text.insert(tk.END, "\nTrinucleotide percentages:\n")
    for combo, perc in tri_percentages.items():
        output_text.insert(tk.END, f"{combo}: {perc:.2f}%\n")
    output_text.config(state='disabled')
    next_button.pack_forget()  

def next_exercise():
    global exercise_index
    exercise_index += 1
    output_text.config(state='normal')
    output_text.delete(1.0, tk.END) 
    output_text.config(state='disabled')
    next_button.pack_forget()
    if exercise_index == 1:
        instruction_label.config(text="Exercise 2: For each combination, find out the percentage inside the sequence")
        solution_button.config(command=exercise2)
        solution_button.pack()
    elif exercise_index == 2:
        instruction_label.config(text="Exercise 3: Show the percentage for each combination in the output")
        solution_button.config(command=exercise3)
        solution_button.pack()

root = tk.Tk()
root.title("Step-by-Step Bioinformatics Homework 2")

instruction_label = tk.Label(root, text="Exercise 1: Build a brute force engine to generate all dinucleotide and trinucleotide combinations",
                             bg="yellow", font=("Arial", 16), wraplength=700, justify="left")
instruction_label.pack(pady=10)

seq_label = tk.Label(root, text=f"DNA Sequence: {sequence}", font=("Arial", 12), wraplength=700, justify="left")
seq_label.pack(pady=5)

solution_button = tk.Button(root, text="Show Solution", font=("Arial", 14), command=exercise1)
solution_button.pack(pady=10)

output_text = scrolledtext.ScrolledText(root, width=100, height=25, font=("Courier", 12))
output_text.pack(pady=10)
output_text.config(state='disabled')

next_button = tk.Button(root, text="Next Exercise", font=("Arial", 14), command=next_exercise)
next_button.pack_forget()

root.mainloop()
