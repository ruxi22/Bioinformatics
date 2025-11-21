import random
import gradio as gr

BASES = "ACGT"


def random_seq(n):
    return ''.join(random.choice(BASES) for _ in range(n))


def reverse_complement(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def make_transposon(name):
    ir_left = random_seq(4)
    ir_right = reverse_complement(ir_left)
    dr = random_seq(3)
    core = random_seq(random.randint(20, 40))

    full = dr + ir_left + core + ir_right + dr
    return name, full


def insert_with_overlap(host, insert_seq, forced_overlap=None):
    L = len(host)
    ins_len = len(insert_seq)

    if forced_overlap:
        start = random.randint(forced_overlap["start"], forced_overlap["end"])
        start = min(start, L - ins_len)
    else:
        start = random.randint(0, L - ins_len)

    host = host[:start] + insert_seq + host[start + ins_len:]
    return host, start, start + ins_len - 1


def embed_transposons():
    host = random_seq(random.randint(200, 400))
    transposons = [make_transposon(f"T{i+1}") for i in range(4)]

    annotations = []
    used = transposons[:random.randint(3, 4)]

    for name, seq in used:
        forced = random.choice(annotations) if annotations and random.random() < 0.6 else None
        host, s, e = insert_with_overlap(host, seq, forced)
        annotations.append({
            "name": name,
            "start": s,
            "end": e,
            "sequence": seq
        })

    return host, annotations, dict(transposons)


def detect_transposons(dna, transposon_dict):
    hits = []
    for name, seq in transposon_dict.items():
        start = 0
        while True:
            pos = dna.find(seq, start)
            if pos == -1:
                break
            hits.append({"name": name, "start": pos, "end": pos + len(seq) - 1})
            start = pos + 1
    return hits



def generate_dna():
    dna, truth, tdict = embed_transposons()

    truth_str = ""
    for t in truth:
        truth_str += f"{t['name']} | start={t['start']} end={t['end']} | len={len(t['sequence'])}\n"

    return dna, truth_str, tdict


def run_detection(dna, tdict):
    if not dna:
        return "No DNA generated yet."

    hits = detect_transposons(dna, tdict)

    if not hits:
        return "No transposons detected."

    out = ""
    for h in hits:
        out += f"{h['name']} | start={h['start']} end={h['end']}\n"

    return out


# App layout
with gr.Blocks(title="Transposon Detector") as demo:
    gr.Markdown("# 🧬 Transposon Detector (Artificial DNA Simulator)")
    gr.Markdown("Generate DNA with biologically accurate transposons (IR + DR), then detect them.")

    dna_box = gr.Textbox(label="Generated DNA", lines=6)
    truth_box = gr.Textbox(label="Inserted Transposons (Ground Truth)", lines=6)
    detect_box = gr.Textbox(label="Detection Results", lines=6)

    tdict_state = gr.State()

    gen_btn = gr.Button("Generate DNA")
    detect_btn = gr.Button("Detect Transposons")

    gen_btn.click(fn=generate_dna, inputs=None,
                  outputs=[dna_box, truth_box, tdict_state])

    detect_btn.click(fn=run_detection,
                     inputs=[dna_box, tdict_state],
                     outputs=detect_box)

demo.launch()
