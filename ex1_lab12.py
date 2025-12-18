import gradio as gr
from math import log
import matplotlib.pyplot as plt
import pandas as pd


motifs = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
]

bases = ["A", "C", "G", "T"]
L = len(motifs[0])
N = len(motifs)
PSEUDOCOUNT = 1


def analyze_sequence(S):
    # 1. Count matrix with pseudocounts
    count = {b: [PSEUDOCOUNT]*L for b in bases}
    for m in motifs:
        for i, c in enumerate(m):
            count[c][i] += 1

    count_df = pd.DataFrame(count, index=range(1, L+1)).T

    # 2–3. Weight / relative frequency matrix
    denom = N + 4 * PSEUDOCOUNT
    freq = {b: [count[b][i]/denom for i in range(L)] for b in bases}
    freq_df = pd.DataFrame(freq, index=range(1, L+1)).T

    # 4. Log-likelihood matrix
    loglik = {b: [log(freq[b][i]/0.25) for i in range(L)] for b in bases}
    loglik_df = pd.DataFrame(loglik, index=range(1, L+1)).T

    # 5. Sliding window scan
    scores = []
    positions = []
    windows = []

    for i in range(len(S) - L + 1):
        window = S[i:i+L]
        score = sum(loglik[window[j]][j] for j in range(L))
        scores.append(score)
        positions.append(i)
        windows.append(window)

    scan_df = pd.DataFrame({
        "Position": positions,
        "Window": windows,
        "Score": scores
    })

    best_idx = scan_df["Score"].idxmax()
    best_row = scan_df.loc[best_idx]

    # Plot
    plt.figure()
    plt.plot(scan_df["Position"], scan_df["Score"])
    plt.axhline(0)
    plt.xlabel("Sliding window position")
    plt.ylabel("Log-likelihood score")
    plt.title("Exon–Intron Boundary Scan (with pseudocounts)")

    conclusion = (
        f"Best candidate at position {int(best_row['Position'])}, "
        f"window {best_row['Window']}, score = {best_row['Score']:.3f}.\n"
        "A positive log-likelihood peak indicates a likely exon–intron boundary."
    )

    return count_df, freq_df, loglik_df, scan_df, conclusion, plt


with gr.Blocks(title="Exon–Intron Boundary Detection") as demo:
    gr.Markdown("# Exon–Intron Boundary Detection")
    gr.Markdown("PWM-based motif scanning using log-likelihood scores")

    gr.Markdown("## Input")
    seq_input = gr.Textbox(
        value="CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA",
        label="DNA sequence S",
        lines=2
    )

    run_btn = gr.Button("Run analysis")

    gr.Markdown("## 1. Count Matrix")
    count_out = gr.Dataframe(interactive=False)

    gr.Markdown("## 2. Weight / Relative Frequency Matrix")
    freq_out = gr.Dataframe(interactive=False)

    gr.Markdown("## 3. Log-Likelihood Matrix")
    loglik_out = gr.Dataframe(interactive=False)

    gr.Markdown("## 4. Sliding Window Scores")
    scan_out = gr.Dataframe(interactive=False)

    gr.Markdown("## 5. Log-Likelihood Plot")
    plot_out = gr.Plot()

    gr.Markdown("## Conclusion")
    conclusion_out = gr.Textbox(lines=4, interactive=False)

    run_btn.click(
        fn=analyze_sequence,
        inputs=seq_input,
        outputs=[count_out, freq_out, loglik_out, scan_out, conclusion_out, plot_out]
    )

demo.launch()