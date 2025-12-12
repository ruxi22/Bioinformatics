import gradio as gr
import numpy as np
import matplotlib.pyplot as plt


def needleman_wunsch(seq1, seq2, match_score, mismatch_score, gap_penalty):
    n, m = len(seq1), len(seq2)
    score_matrix = np.zeros((m + 1, n + 1))

    for i in range(m + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(n + 1):
        score_matrix[0][j] = j * gap_penalty

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = score_matrix[i - 1][j - 1] + (
                match_score if seq1[j - 1] == seq2[i - 1] else mismatch_score
            )
            up = score_matrix[i - 1][j] + gap_penalty
            left = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(diag, up, left)

    align1, align2 = "", ""
    i, j = m, n
    path = [(i, j)]
    matches = 0

    while i > 0 or j > 0:
        current = score_matrix[i][j]

        if i > 0 and j > 0:
            diag_score = score_matrix[i - 1][j - 1] + (
                match_score if seq1[j - 1] == seq2[i - 1] else mismatch_score
            )
            if current == diag_score:
                align1 += seq1[j - 1]
                align2 += seq2[i - 1]
                if seq1[j - 1] == seq2[i - 1]:
                    matches += 1
                i -= 1
                j -= 1
                path.append((i, j))
                continue

        if i > 0 and current == score_matrix[i - 1][j] + gap_penalty:
            align1 += "-"
            align2 += seq2[i - 1]
            i -= 1
            path.append((i, j))
            continue

        align1 += seq1[j - 1]
        align2 += "-"
        j -= 1
        path.append((i, j))

    return score_matrix, align1[::-1], align2[::-1], path, matches


def process_alignment(s1, s2, gap, match, mismatch):
    matrix, a1, a2, path, matches = needleman_wunsch(s1, s2, match, mismatch, gap)

    length = len(a1)
    similarity = int((matches / length) * 100) if length else 0

    bars = "".join("|" if a1[i] == a2[i] and a1[i] != "-" else " " for i in range(length))

    text = (
        f"{a1}\n"
        f"{bars}\n"
        f"{a2}\n\n"
        f"Matches: {matches}\n"
        f"Length: {length}\n"
        f"Similarity: {similarity}%"
    )

    fig1, ax1 = plt.subplots(figsize=(7, 5))
    im = ax1.imshow(matrix, cmap="magma", aspect="auto")
    ax1.set_title("Alignment score matrix")
    fig1.colorbar(im)
    plt.tight_layout()

    fig2, ax2 = plt.subplots(figsize=(7, 5))
    grid = np.ones((matrix.shape[0], matrix.shape[1], 3))
    grid[:] = [1, 1, 0.9]

    for r, c in path:
        grid[r, c] = [0.85, 0.2, 0.2]

    ax2.imshow(grid)
    ax2.set_title("Traceback path")
    ax2.set_xticks([])
    ax2.set_yticks([])
    plt.tight_layout()

    return text, fig1, fig2


with gr.Blocks(title="DNA Alignment") as app:
    gr.Markdown(
        """
        ## Needleman–Wunsch DNA Alignment
        Global alignment with adjustable scoring parameters.
        """
    )

    with gr.Row():
        with gr.Column(scale=1):
            s1 = gr.Textbox(
                label="Sequence 1",
                value="ACCGTGAAGCCAATAC",
                elem_id="mono"
            )
            s2 = gr.Textbox(
                label="Sequence 2",
                value="AGCGTGCAGCCAATAC",
                elem_id="mono"
            )

            gr.Markdown("### Scoring")
            match = gr.Slider(-5, 5, value=1, step=1, label="Match")
            mismatch = gr.Slider(-5, 0, value=-1, step=1, label="Mismatch")
            gap = gr.Slider(-5, 0, value=0, step=1, label="Gap")

            run = gr.Button("Run alignment")

        with gr.Column(scale=2):
            result = gr.Textbox(
                label="Alignment",
                lines=7,
                show_copy_button=True,
                elem_id="mono"
            )
            heatmap = gr.Plot(label="Score matrix")
            traceback = gr.Plot(label="Traceback")

    app.css = """
    #mono textarea {
        font-family: Courier New, monospace;
        font-size: 14px;
    }
    """

    run.click(
        process_alignment,
        inputs=[s1, s2, gap, match, mismatch],
        outputs=[result, heatmap, traceback]
    )

if __name__ == "__main__":
    app.launch()
