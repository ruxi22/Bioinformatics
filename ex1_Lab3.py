import customtkinter as ctk
from PIL import Image, ImageTk
import os
import math


ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("green")


def tm_simple(dna):
    dna = dna.upper()
    return 4 * (dna.count('G') + dna.count('C')) + 2 * (dna.count('A') + dna.count('T'))

def tm_salt_adjusted(dna, na_conc=0.05):
    dna = dna.upper()
    length = len(dna)
    gc_percent = ((dna.count('G') + dna.count('C')) / length) * 100
    return 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - 600/length


def calculate_tm():
    dna = dna_entry.get().strip().upper()
    na_conc_text = na_entry.get().strip()

    if not dna or any(base not in "ATGC" for base in dna):
        result_label.configure(text="⚠️ Invalid DNA sequence!\nUse only A, T, G, C.", text_color="red")
        return

    try:
        na_conc = float(na_conc_text) if na_conc_text else 0.05
    except ValueError:
        result_label.configure(text="⚠️ Na+ concentration must be a number.", text_color="red")
        return

    tm1 = tm_simple(dna)
    tm2 = tm_salt_adjusted(dna, na_conc)
    gc_percent = ((dna.count('G') + dna.count('C')) / len(dna)) * 100

    result_label.configure(
        text=f"Tm (Simple): {tm1:.2f} °C\n"
             f"Tm (Salt-adjusted): {tm2:.2f} °C\n"
             f"GC Content: {gc_percent:.2f}%",
        text_color="#00FF7F"
    )


app = ctk.CTk()
app.title("DNA Melting Temperature Calculator")
app.geometry("600x500")


script_dir = r"C:\Users\ruxan\Desktop\L3"
image_path = os.path.join(script_dir, "dna.png")

try:
    image = Image.open(image_path)
    image = image.resize((500, 500), Image.Resampling.LANCZOS)
    dna_image = ImageTk.PhotoImage(image)
    image_label = ctk.CTkLabel(app, image=dna_image, text="")
    image_label.pack(pady=10)
except Exception as e:
    print("Could not load image:", e)


header = ctk.CTkLabel(app, text="DNA Tm Calculator", font=ctk.CTkFont(size=28, weight="bold"), text_color="#FFD700")
header.pack(pady=10)


input_card = ctk.CTkFrame(app, corner_radius=20, fg_color="#1E1E2F")
input_card.pack(pady=10, padx=30, fill="x")

dna_label = ctk.CTkLabel(input_card, text="DNA Sequence:", font=ctk.CTkFont(size=14))
dna_label.pack(pady=(10, 2))
dna_entry = ctk.CTkEntry(input_card, width=400, placeholder_text="e.g., ATGCGTAC")
dna_entry.pack(pady=(0, 10))

na_label = ctk.CTkLabel(input_card, text="Na+ Concentration (M, optional):", font=ctk.CTkFont(size=14))
na_label.pack(pady=(5, 2))
na_entry = ctk.CTkEntry(input_card, width=200, placeholder_text="Default 0.05")
na_entry.pack(pady=(0, 10))


calc_button = ctk.CTkButton(app, text="Calculate Tm", command=calculate_tm, width=250, height=45,
                            font=ctk.CTkFont(size=16, weight="bold"), fg_color="#00BFFF", hover_color="#1E90FF")
calc_button.pack(pady=15)


result_card = ctk.CTkFrame(app, corner_radius=20, fg_color="#2E2E3F")
result_card.pack(pady=10, padx=30, fill="both", expand=True)

result_label = ctk.CTkLabel(result_card, text="", font=ctk.CTkFont(size=16), justify="left", text_color="#00FF7F")
result_label.pack(pady=20, padx=20)


app.mainloop()
