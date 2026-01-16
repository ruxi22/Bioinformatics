import tkinter as tk
from tkinter import messagebox


def parse_matrix(text):
    rows = text.strip().split("\n")
    return [list(map(float, row.split())) for row in rows]

def parse_vector(text):
    return list(map(float, text.split()))

def transpose(A):
    return list(map(list, zip(*A)))

def multiply(A, x):
    return [
        sum(A[i][j] * x[j] for j in range(len(x)))
        for i in range(len(A))
    ]

def predict():
    try:
        A = parse_matrix(matrix_input.get("1.0", tk.END))
        x = parse_vector(vector_input.get())

        n = len(A)
        if any(len(row) != n for row in A):
            raise ValueError("Matrix must be square")
        if len(x) != n:
            raise ValueError("Vector size must match matrix size")

        AT = transpose(A)

        output.delete("1.0", tk.END)

        for step in range(1, 6):
            x = multiply(AT, x)
            vector_str = "[" + ", ".join(f"{v:.6f}" for v in x) + "]"
            output.insert(tk.END, f"Step {step}: {vector_str}\n")

    except Exception as e:
        messagebox.showerror("Error", str(e))


root = tk.Tk()
root.title("Markov Prediction (5 Steps)")
root.resizable(False, False)

main = tk.Frame(root, padx=12, pady=12)
main.pack()

# Matrix input
tk.Label(main, text="Matrix A (rows sum to 1):").grid(row=0, column=0, sticky="w")
matrix_input = tk.Text(main, height=6, width=40)
matrix_input.grid(row=1, column=0, padx=(0, 20))

matrix_input.insert(
    "1.0",
    "0.3 0.35 0.35\n"
    "0 0 1\n"
    "0.9 0 0.1"
)

# Vector input
right = tk.Frame(main)
right.grid(row=1, column=1, sticky="n")

tk.Label(right, text="Initial vector x₀:").pack(anchor="w")
vector_input = tk.Entry(right, width=25)
vector_input.pack(pady=(0, 12))
vector_input.insert(0, "1 0 0")

tk.Button(
    right,
    text="Predict 5 Steps",
    command=predict,
    width=20
).pack()

# Output
tk.Label(main, text="Prediction output:").grid(
    row=2, column=0, columnspan=2, sticky="w", pady=(12, 0)
)

output = tk.Text(
    main,
    height=8,
    width=70,
    font=("Courier", 10)
)
output.grid(row=3, column=0, columnspan=2)

root.mainloop()
