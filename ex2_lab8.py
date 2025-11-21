import gradio as gr

def reverse_complement(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]

def find_inverted_repeats(seq, min_len=4, max_len=6):
    irs = []
    n = len(seq)

    for L in range(min_len, max_len + 1):
        for i in range(n - L):
            left = seq[i:i + L]
            right = reverse_complement(left)

            j = seq.find(right, i + L)
            if j != -1:
                irs.append((i, j, L))

    return irs


def detect_transposons_from_ir(genome):
    irs = find_inverted_repeats(genome, 4, 6)
    elements = []

    for left, right, L in irs:
        elements.append({
            "start": left,
            "end": right + L,
            "ir_length": L,
            "length": (right + L) - left
        })

    return elements


def format_output(elements):
    if not elements:
        return "No transposable elements detected."

    out = "Detected Transposable Elements:\n\n"
    for e in elements[:5000]:
        out += (
            f"Start: {e['start']}  |  "
            f"End: {e['end']}  |  "
            f"IR size: {e['ir_length']}  |  "
            f"Length: {e['length']}\n"
        )
    return out

def run_detection(fasta_text):
    if not fasta_text:
        return "Paste a FASTA sequence first!"

    lines = fasta_text.splitlines()
    genome = "".join(l.strip() for l in lines if not l.startswith(">"))

    if len(genome) < 50:
        return "FASTA too short or invalid."

    elements = detect_transposons_from_ir(genome)
    return format_output(elements)


with gr.Blocks(title="Transposon Detector (FASTA Paste)") as demo:
    gr.Markdown("# 🧬 Transposon Detector (Paste FASTA)")
    gr.Markdown("""
Paste a **FASTA genome sequence** directly below.

The app will:
- clean FASTA,
- detect inverted repeats (4–6bp),
- detect transposable elements,
- output start/end/length.

No download needed.
""")

    fasta_box = gr.Textbox(
        label="Paste FASTA Here:",
        placeholder=">my_genome\nACGTACGT...",
        lines=20
    )

    output_box = gr.Textbox(
        label="Results",
        lines=25
    )

    detect_button = gr.Button("Detect Transposons")

    detect_button.click(
        fn=run_detection,
        inputs=fasta_box,
        outputs=output_box
    )

demo.launch()
