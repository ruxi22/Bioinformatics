import tkinter as tk
from tkinter import filedialog, messagebox
import json
import random


def weighted_choice(prob_map):
    r = random.random()
    acc = 0.0
    for k, p in prob_map.items():
        acc += p
        if r <= acc:
            return k
    return next(iter(prob_map))

class MarkovGeneratorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Markov Generator (DNA / English)")
        self.root.resizable(False, False)

        self.model = None
        self.mode = tk.StringVar(value="DNA")

        self.build_ui()

    def build_ui(self):
        frame = tk.Frame(self.root, padx=12, pady=12)
        frame.pack()

        # Mode
        mode_frame = tk.Frame(frame)
        mode_frame.pack(fill="x")

        tk.Label(mode_frame, text="Mode:").pack(side="left")
        tk.Radiobutton(
            mode_frame, text="DNA", variable=self.mode, value="DNA"
        ).pack(side="left")
        tk.Radiobutton(
            mode_frame, text="English text", variable=self.mode, value="TEXT"
        ).pack(side="left")

        tk.Button(
            mode_frame, text="Load JSON", command=self.load_model, width=14
        ).pack(side="right")

        self.status = tk.StringVar(value="No model loaded")
        tk.Label(frame, textvariable=self.status).pack(anchor="w", pady=(4, 8))

        # Options
        opt = tk.Frame(frame)
        opt.pack(fill="x")

        tk.Label(opt, text="Length:").grid(row=0, column=0, sticky="w")
        self.length_entry = tk.Entry(opt, width=10)
        self.length_entry.grid(row=0, column=1, padx=(6, 20))
        self.length_entry.insert(0, "50")

        tk.Label(opt, text="Start (optional):").grid(row=0, column=2, sticky="w")
        self.start_entry = tk.Entry(opt, width=18)
        self.start_entry.grid(row=0, column=3, padx=(6, 0))

        tk.Button(frame, text="Generate", command=self.generate, width=18).pack(pady=8)

        self.output = tk.Text(frame, height=12, width=90, font=("Courier", 10))
        self.output.pack()


    def load_model(self):
        path = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
        if not path:
            return

        try:
            with open(path, "r", encoding="utf-8") as f:
                self.model = json.load(f)
            self.status.set(f"Loaded: {path.split('/')[-1]}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def generate(self):
        if self.model is None:
            messagebox.showerror("Error", "Load a JSON model first")
            return

        try:
            length = int(self.length_entry.get())
            if length < 1:
                raise ValueError
        except Exception:
            messagebox.showerror("Error", "Length must be a positive integer")
            return

        if self.mode.get() == "DNA":
            self.generate_dna(length)
        else:
            self.generate_text(length)


    def generate_dna(self, length):
        start = self.start_entry.get().strip().upper()
        states = list(self.model.keys())

        if start == "":
            start = random.choice(states)
        if start not in self.model:
            messagebox.showerror("Error", "Invalid start symbol for DNA")
            return

        seq = [start]
        cur = start
        for _ in range(length - 1):
            cur = weighted_choice(self.model[cur])
            seq.append(cur)

        self.output.delete("1.0", tk.END)
        self.output.insert(tk.END, "".join(seq))


    def generate_text(self, length):
        word_to_id = self.model["word_to_id"]
        id_to_word = {int(k): v for k, v in self.model["id_to_word"].items()}
        transition = {
            int(k): {int(j): p for j, p in v.items()}
            for k, v in self.model["transition_matrix"].items()
        }

        start_word = self.start_entry.get().strip().lower()

        if start_word:
            if start_word not in word_to_id:
                messagebox.showerror("Error", "Start word not in model")
                return
            cur = int(word_to_id[start_word])
        else:
            cur = random.choice([i for i in transition if transition[i]])

        words = [id_to_word[cur]]

        for _ in range(length - 1):
            row = transition.get(cur, {})
            if not row:
                cur = random.choice([i for i in transition if transition[i]])
            else:
                cur = weighted_choice(row)
            words.append(id_to_word[cur])

        self.output.delete("1.0", tk.END)
        self.output.insert(tk.END, " ".join(words))


root = tk.Tk()
app = MarkovGeneratorApp(root)
root.mainloop()
