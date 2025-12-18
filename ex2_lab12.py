import gradio as gr
import matplotlib.pyplot as plt
import os
import tempfile


motif = "AGGTAAAGT"
L = len(motif)
bases = ["A", "C", "G", "T"]

def score_window(window):
    return sum(1 if w == m else -1 for w, m in zip(window, motif))


def read_fasta(file):
    seq = ""
    with open(file.name, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip().upper()
    return seq


def scan_genome(seq):
    positions, scores = [], []
    for i in range(len(seq) - L + 1):
        window = seq[i:i + L]
        if any(b not in bases for b in window):
            continue
        positions.append(i)
        scores.append(score_window(window))
    return positions, scores


def analyze_genomes(files):
    if files is None:
        return [None] * 10

    images = []
    tmpdir = tempfile.mkdtemp()
    BIN_SIZE = 100  # netezire semnal

    for idx, f in enumerate(files[:10]):
        seq = read_fasta(f)
        if len(seq) < L:
            images.append(None)
            continue

        positions, scores = scan_genome(seq)

        # ---- binning (smoothing) ----
        binned_pos = []
        binned_scores = []
        for i in range(0, len(scores), BIN_SIZE):
            chunk = scores[i:i + BIN_SIZE]
            if chunk:
                binned_scores.append(sum(chunk) / len(chunk))
                binned_pos.append(positions[i])

        # ---- plot ----
        filename = os.path.basename(f.name)

        plt.figure(figsize=(12, 4), dpi=150)
        plt.plot(binned_pos, binned_scores)
        plt.axhline(0, linestyle="--", linewidth=1)

        plt.xlabel("Genome position")
        plt.ylabel("Signal score")
        plt.title(f"Motif signal – {filename}")

        y_min = min(binned_scores)
        y_max = max(binned_scores)
        plt.ylim(y_min - 0.5, y_max + 0.5)

        plt.tight_layout()

        out_path = os.path.join(tmpdir, f"genome_{idx}.png")
        plt.savefig(out_path)
        plt.close()

        images.append(out_path)

    # completăm până la 10 outputs (cerință Gradio)
    while len(images) < 10:
        images.append(None)

    return images



with gr.Blocks(title="Exercise 2 – Influenza Genome Motif Scan") as demo:
    gr.Markdown("# Exercise 2 – Influenza Genome Motif Scan")
    gr.Markdown(
        "Upload up to 10 genome FASTA files. Each genome is scanned independently. "
        "For each genome, a separate signal plot is generated showing the most likely "
        "locations of functional motifs."
    )

    genome_files = gr.File(
        label="Upload genome FASTA files (max 10)",
        file_types=[".fa", ".fasta", ".txt"],
        file_count="multiple"
    )

    run_btn = gr.Button("Run scan")

    gr.Markdown("## Signal plots (one per genome)")

    img_outputs = [
        gr.Image(label=f"Genome {i+1}", interactive=False)
        for i in range(10)
    ]

    run_btn.click(
        fn=analyze_genomes,
        inputs=genome_files,
        outputs=img_outputs
    )

demo.launch()
