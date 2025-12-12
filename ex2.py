import gradio as gr
import numpy as np
import matplotlib.pyplot as plt
import requests


def smith_waterman_score_native(seq1, seq2, match=3, mismatch=-3, gap=-2):
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    prev = np.zeros(cols, dtype=int)
    curr = np.zeros(cols, dtype=int)
    max_score = 0

    for i in range(1, rows):
        for j in range(1, cols):
            diag = prev[j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            up = prev[j] + gap
            left = curr[j - 1] + gap
            curr[j] = max(0, diag, up, left)
            max_score = max(max_score, curr[j])
        prev[:] = curr[:]

    return max_score


def fetch_ncbi_sequence(accession_id):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        f"efetch.fcgi?db=nuccore&id={accession_id}&rettype=fasta&retmode=text"
    )
    try:
        r = requests.get(url)
        r.raise_for_status()
        lines = r.text.strip().splitlines()
        seq = "".join(lines[1:]).upper()
        return seq, f"Fetched {accession_id} ({len(seq)} bp)"
    except Exception as e:
        return None, str(e)


def chunk_sequence(seq, size):
    return [seq[i:i + size] for i in range(0, len(seq), size)]


def run_windowed_alignment(seq1_id, seq2_id, window_size, progress=gr.Progress()):
    progress(0, desc="Fetching genome 1")
    s1, msg1 = fetch_ncbi_sequence(seq1_id)
    if not s1:
        return msg1, None

    progress(0.1, desc="Fetching genome 2")
    s2, msg2 = fetch_ncbi_sequence(seq2_id)
    if not s2:
        return msg2, None

    w = int(window_size)
    c1 = chunk_sequence(s1, w)
    c2 = chunk_sequence(s2, w)

    mat = np.zeros((len(c2), len(c1)))
    total = mat.size
    k = 0

    for i in range(len(c2)):
        for j in range(len(c1)):
            mat[i, j] = smith_waterman_score_native(c2[i], c1[j])
            k += 1
            if k % 50 == 0:
                progress(0.1 + 0.9 * k / total, desc="Aligning windows")

    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.imshow(mat, cmap="inferno", aspect="auto", interpolation="nearest")
    ax.set_title(f"Local similarity map (window = {w} bp)")
    ax.set_xlabel(seq1_id)
    ax.set_ylabel(seq2_id)
    plt.colorbar(im, ax=ax, label="Smith–Waterman score")
    plt.tight_layout()

    text = (
        f"Done.\n"
        f"{seq1_id}: {len(s1)} bp → {len(c1)} windows\n"
        f"{seq2_id}: {len(s2)} bp → {len(c2)} windows\n"
        f"Method: windowed Smith–Waterman (native)"
    )

    return text, fig


with gr.Blocks(title="Viral Genome Local Alignment") as app:
    gr.Markdown(
        """
        ## Viral genome local alignment
        Windowed Smith–Waterman comparison using NCBI reference genomes.
        """
    )

    with gr.Row():
        with gr.Column(scale=1):
            gr.Markdown("### Input")

            id1 = gr.Textbox(
                label="Genome 1 (NCBI accession)",
                value="NC_045512.2",
                elem_id="mono"
            )
            id2 = gr.Textbox(
                label="Genome 2 (NCBI accession)",
                value="NC_026433.1",
                elem_id="mono"
            )

            window = gr.Slider(
                50, 500, step=50, value=200,
                label="Window size (bp)"
            )

            run = gr.Button("Run alignment")

        with gr.Column(scale=2):
            out_text = gr.Textbox(
                label="Status",
                lines=5,
                elem_id="mono"
            )
            out_plot = gr.Plot(label="Similarity heatmap")

    app.css = """
    #mono textarea {
        font-family: Courier New, monospace;
        font-size: 13px;
    }
    """

    run.click(
        run_windowed_alignment,
        inputs=[id1, id2, window],
        outputs=[out_text, out_plot]
    )

if __name__ == "__main__":
    app.launch()
