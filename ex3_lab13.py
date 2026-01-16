import tkinter as tk
from tkinter import messagebox, filedialog
import json
import string

STOP_WORDS = {
    "the", "is", "a", "an", "to", "of", "and", "in", "on", "for", "with",
    "this", "that", "it", "as", "be", "are", "was", "were", "by", "or",
    "from", "at", "which", "its"
}

def tokenize(text):
    text = text.lower()
    for ch in string.punctuation:
        text = text.replace(ch, " ")
    return [w for w in text.split() if w and w not in STOP_WORDS]

def build_model(words):
    word_to_id = {}
    id_to_word = {}

    for w in words:
        if w not in word_to_id:
            idx = len(word_to_id)
            word_to_id[w] = idx
            id_to_word[idx] = w

    n = len(word_to_id)

    counts = {i: {} for i in range(n)}

    for i in range(len(words) - 1):
        a = word_to_id[words[i]]
        b = word_to_id[words[i + 1]]
        counts[a][b] = counts[a].get(b, 0) + 1

    matrix = {}
    for i in counts:
        total = sum(counts[i].values())
        matrix[i] = {
            j: counts[i][j] / total
            for j in counts[i]
        } if total > 0 else {}

    return word_to_id, id_to_word, matrix

def calculate():
    text = text_input.get("1.0", tk.END).strip()

    if len(text) < 300:
        messagebox.showerror(
            "Error",
            "Text must contain at least 300 characters"
        )
        return

    words = tokenize(text)

    if len(words) < 2:
        messagebox.showerror(
            "Error",
            "Not enough words after filtering"
        )
        return

    word_to_id, id_to_word, matrix = build_model(words)

    output.delete("1.0", tk.END)
    output.insert(
        tk.END,
        f"Words used: {len(words)}\n"
        f"Unique symbols: {len(word_to_id)}\n\n"
    )

    for i in list(matrix.keys())[:12]:
        transitions = ", ".join(
            f"{j}:{matrix[i][j]:.3f}"
            for j in matrix[i]
        )
        output.insert(
            tk.END,
            f"{i} ({id_to_word[i]}) → {transitions}\n"
        )

    output.model = {
        "word_to_id": word_to_id,
        "id_to_word": id_to_word,
        "transition_matrix": matrix
    }

def save_json():
    if not hasattr(output, "model"):
        messagebox.showerror("Error", "Nothing to save")
        return

    path = filedialog.asksaveasfilename(
        defaultextension=".json",
        filetypes=[("JSON files", "*.json")]
    )

    if not path:
        return

    with open(path, "w", encoding="utf-8") as f:
        json.dump(output.model, f, indent=4)

    messagebox.showinfo("Saved", "Model saved as JSON")


root = tk.Tk()
root.title("Word Transition Model (Integer Symbols)")
root.resizable(False, False)

main = tk.Frame(root, padx=12, pady=12)
main.pack()

tk.Label(main, text="English text (min. 300 characters):").grid(
    row=0, column=0, sticky="w"
)

text_input = tk.Text(main, height=16, width=90)
text_input.grid(row=1, column=0, columnspan=2)

# --- Your text ---
text_input.insert(
    "1.0",
    "Reading isn't just decoding words on a page; it's an active journey, "
    "a portal to endless worlds and experiences, often costing nothing more "
    "than time and curiosity. Whether it's a thrilling novel that transports "
    "you to distant galaxies, a historical account that illuminates the past, "
    "or a poem that captures the essence of human emotion, books offer "
    "unparalleled escapism and insight. They allow us to walk in someone "
    "else's shoes, understand different cultures, and develop empathy in ways "
    "few other activities can, expanding our perspective beyond our immediate "
    "reality. The simple act of opening a book can transform a quiet afternoon "
    "into an adventure. A well crafted story builds worlds, introduces "
    "unforgettable characters, and presents challenges that mirror our own "
    "lives, providing comfort and understanding. Beyond entertainment, "
    "reading fosters critical thinking and lifelong learning."
)

tk.Button(
    main,
    text="Calculate model",
    command=calculate,
    width=25
).grid(row=2, column=0, pady=6, sticky="w")

tk.Button(
    main,
    text="Save JSON",
    command=save_json,
    width=20
).grid(row=2, column=1, pady=6, sticky="e")

tk.Label(main, text="Preview:").grid(
    row=3, column=0, columnspan=2, sticky="w", pady=(8, 0)
)

output = tk.Text(
    main,
    height=12,
    width=90,
    font=("Courier", 10)
)
output.grid(row=4, column=0, columnspan=2)

root.mainloop()
