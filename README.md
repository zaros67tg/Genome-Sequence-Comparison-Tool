# 🧬 Genome Sequence Comparison Tool

A production-ready, interactive **bioinformatics web application** built with Python and Streamlit for comparing biological sequences — DNA, RNA, and protein. Perform pairwise alignments, quantify similarity, classify mutations, explore interactive visualizations, and export comprehensive PDF and CSV reports — all from your browser.

youtube link- https://youtu.be/qOKotSAtG3A
---

## Table of Contents

- [Features](#features)
- [Project Structure](#project-structure)
- [Installation & Setup](#installation--setup)
- [Usage](#usage)
- [Alignment Algorithms](#alignment-algorithms)
- [Scoring Schemes](#scoring-schemes)
- [Sample Datasets](#sample-datasets)
- [Exporting Results](#exporting-results)
- [Tech Stack](#tech-stack)
- [Known Limitations](#known-limitations)

---

## Features

### 🔬 Sequence Input & Validation
- Accept **FASTA / FA / TXT** file uploads or paste raw sequences directly into the app
- Parse multi-record FASTA files with `Bio.SeqIO`
- Strict per-mode character validation:
  - **DNA** — `A C G T N`
  - **RNA** — `A C G U N`
  - **Protein** — standard 20 single-letter amino acid codes
- Live validation badge (✅ / ❌) with detailed error messages and invalid character listing

### ⚡ Pairwise Alignment Engine
- **Global alignment** — Needleman-Wunsch end-to-end alignment
- **Local alignment** — Smith-Waterman best sub-sequence alignment
- Powered by `Bio.Align.PairwiseAligner` (modern Biopython ≥ 1.80 API)
- Fully configurable scoring parameters via the sidebar
- **BLOSUM62 substitution matrix** applied automatically for protein sequences — no flat match/mismatch penalties that would treat Leu↔Ile the same as Leu↔Pro
- Results cached with `@st.cache_data` — switching tabs never re-runs the aligner

### 📊 Similarity Metrics
- Sequence lengths (Seq 1, Seq 2, and alignment length)
- Alignment score
- **Percentage identity** (matches / alignment length × 100, BLAST-style)
- Match count, mismatch count, total gap positions, gap-open events
- **GC content (%)** for DNA / RNA sequences

### 🧪 Mutation Detection
Column-by-column classification of every non-match alignment position:

| Type | Definition |
|---|---|
| **Substitution** | Non-gap residue in Seq 1 differs from non-gap residue in Seq 2 |
| **Insertion** | Gap in Seq 1, residue in Seq 2 |
| **Deletion** | Residue in Seq 1, gap in Seq 2 |
| **Duplication** | Insertion or deletion flanked by ≥ 2 consecutive identical residues (tandem repeat heuristic) |

Output is a fully scrollable, colour-coded Pandas DataFrame.

### 📈 Interactive Visualizations (Plotly)
1. **Alignment Match / Mismatch Map** — horizontal scatter strip: 🟢 Match · 🔴 Mismatch · 🟠 Gap
2. **Dot Plot** — k-mer window comparison of Seq 1 vs Seq 2 (adjustable window size 1–10 bp/aa, capped at 2 000 bp for performance)
3. **Base / Residue Composition Chart** — grouped bar chart of nucleotide or amino acid frequencies
4. **Mutation Type Distribution** — donut pie chart of Substitution / Insertion / Deletion / Duplication counts

### 📥 Export & Reporting
- **CSV** — full mutation events table (UTF-8, Excel-compatible)
- **PDF** — multi-page A4 report (ReportLab) containing:
  - Metadata (date, sequence IDs, alignment mode, scoring parameters)
  - Similarity metrics table
  - Alignment preview (first 1 000 characters)
  - Embedded alignment map chart (requires `kaleido`)
  - Mutation summary table
  - Colour-coded detailed mutation table

---

## Project Structure

```
Genome-Sequence-Comparison-Tool/
│
├── app.py                   # Main Streamlit dashboard — 5-tab layout + sidebar
├── sequence_reader.py       # FASTA parsing (Bio.SeqIO) and sequence validation
├── alignment.py             # PairwiseAligner engine — NW global / SW local / BLOSUM62
├── similarity.py            # Identity %, GC content, match/mismatch/gap metrics
├── mutation_detection.py    # Substitution / Insertion / Deletion / Duplication classifier
├── visualization.py         # Plotly charts — mismatch map, dot plot, composition, pie
├── reports.py               # PDF (ReportLab) and CSV export utilities
├── requirements.txt         # Python dependencies
├── verify_all.py            # End-to-end test suite (edge cases + happy path)
│
└── datasets/                # Bundled sample FASTA files for quick testing
    ├── dna_sample_1.fasta      # Reference DNA sequence
    ├── dna_sample_2.fasta      # Mutant DNA (3 SNPs)
    ├── rna_sample_1.fasta      # Reference RNA (mRNA fragment)
    ├── rna_sample_2.fasta      # Isoform RNA (2 SNPs + insertion)
    ├── protein_sample_1.fasta  # Reference TP53 domain fragment
    └── protein_sample_2.fasta  # Oncogenic mutant variant
```

---

## Installation & Setup

### Prerequisites

- **Python 3.11 or higher** — [Download Python](https://www.python.org/downloads/)
- `pip` (bundled with Python)

### Step 1 — Clone the repository

```bash
git clone https://github.com/your-username/genome-sequence-comparison-tool.git
cd genome-sequence-comparison-tool/Genome-Sequence-Comparison-Tool
```

Or download and unzip the project folder directly.

### Step 2 — (Recommended) Create a virtual environment

```bash
# Create
python -m venv .venv

# Activate — Windows
.venv\Scripts\activate

# Activate — macOS / Linux
source .venv/bin/activate
```

### Step 3 — Install dependencies

All required packages are listed in [`requirements.txt`](requirements.txt):

```bash
pip install -r requirements.txt
```

This installs:

| Package | Purpose |
|---|---|
| `streamlit` | Web application framework |
| `biopython` | Sequence parsing (`Bio.SeqIO`) and alignment (`Bio.Align`) |
| `pandas` | Mutation table and metrics DataFrames |
| `numpy` | Numerical operations |
| `plotly` | Interactive charts |
| `reportlab` | PDF report generation |
| `kaleido` | Static PNG export of Plotly charts for PDF embedding |

> **Note:** `kaleido` is optional for running the app — it is only required to embed the alignment chart image inside the downloaded PDF report. The app runs and generates PDFs without it (a fallback note is shown in the PDF instead of the chart).

---

## Usage

### Running the app locally

From inside the `Genome-Sequence-Comparison-Tool/` directory, run:

```bash
python -m streamlit run app.py
```

Streamlit will print a local URL — open it in your browser:

```
  Local URL:  http://localhost:8501
  Network URL: http://<your-ip>:8501
```

---

### App Walkthrough

#### Sidebar (always visible)

| Control | Description |
|---|---|
| **Sequence Type** | Switch between DNA, RNA, and Protein modes |
| **Alignment Algorithm** | Global (Needleman-Wunsch) or Local (Smith-Waterman) |
| **Scoring Parameters** | Match, mismatch, gap-open, gap-extend penalties (hidden for Protein — BLOSUM62 is used automatically) |
| **Load Seq 1 / Load Seq 2** | Instantly load the bundled sample FASTA for the selected mode |
| **Reset All** | Clear all sequences and results |

---

#### Tab 1 — Sequence Input

1. **Upload** a `.fasta`, `.fa`, or `.txt` file — or **paste** a raw / FASTA-formatted sequence into the text area
2. Both inputs support multi-record FASTA files; the first record is used by default
3. A live validation badge confirms the sequence is valid for the selected mode
4. When both sequences are valid, the **▶ Run Alignment** button becomes active
5. Click **▶ Run Alignment** to execute the full analysis pipeline

> **Tip:** Click **Load Seq 1** and **Load Seq 2** in the sidebar to instantly load the bundled sample FASTA files and try the tool right away.

---

#### Tab 2 — Alignment View

- **KPI cards** — alignment score, alignment length, mode
- **Colour-coded alignment** — each column is rendered in HTML:
  - 🟢 **Green** — identical residues (match)
  - 🔴 **Red** — substitution (mismatch)
  - 🟠 **Orange** — gap in either sequence
- **Formatted alignment** — classic three-line representation (Seq1 / match-line / Seq2), wrapped at 60 columns

---

#### Tab 3 — Metrics & Mutations

- **4 KPI cards** — identity %, match count, mismatch count, gap positions
- **Detailed metrics table** — all 11 similarity metrics including GC content
- **Mutation summary** — counts by type (Substitution / Insertion / Deletion / Duplication)
- **Full mutation events table** — every non-match position with 1-based position, both bases, mutation type, and description; colour-coded by mutation type

---

#### Tab 4 — Visualizations

- **Alignment Match / Mismatch Map** — scrollable Plotly scatter strip
- **Dot Plot** — use the slider to adjust k-mer window size (1–10); larger windows reduce background noise for repetitive sequences
- **Base Composition Chart** — compare nucleotide or amino acid frequencies side by side
- **Mutation Distribution** — donut pie chart of mutation type counts

All charts are fully interactive: hover for tooltips, zoom, pan, and save as PNG via the Plotly toolbar (📷 icon).

---

#### Tab 5 — Downloads

| Button | Output |
|---|---|
| **⬇ Download mutations.csv** | Mutation events table as a UTF-8 CSV file |
| **⬇ Download report.pdf** | Multi-page A4 PDF with metrics, alignment preview, embedded chart, and mutation tables |

The PDF is generated on demand each time you click the button — nothing is written to disk on the server.

---

### Running the test suite

To verify all modules and edge cases:

```bash
python verify_all.py
```

Expected output ends with:

```
=== ALL TESTS PASSED ===
```

The suite covers: syntax check of all 7 modules, DNA/RNA/Protein happy paths, BLOSUM62 vs flat scoring verification, empty input, invalid characters, multi-header FASTA, cross-type validation, 5 kb local alignment performance, duplication heuristic boundary conditions, and PDF generation.

---

## Alignment Algorithms

| Mode | Algorithm | Best for |
|---|---|---|
| **Global** | Needleman-Wunsch | Sequences of similar length; full end-to-end comparison |
| **Local** | Smith-Waterman | Finding the best matching sub-region; highly divergent sequences; large inputs |

---

## Scoring Schemes

### DNA / RNA (flat scoring)

| Parameter | Default | Effect |
|---|---|---|
| Match score | `+2.0` | Reward for identical residues |
| Mismatch penalty | `-1.0` | Penalty for substitutions |
| Gap open penalty | `-2.0` | Penalty for starting a new gap |
| Gap extend penalty | `-0.5` | Per-position cost inside an open gap |

Adjust all four via the sidebar. Increasing the gap-open penalty discourages fragmented alignments; reducing it allows more indels.

### Protein (BLOSUM62)

When **Protein** mode is selected, the BLOSUM62 substitution matrix replaces flat scores automatically. BLOSUM62 encodes empirical log-odds substitution frequencies derived from conserved protein block alignments, correctly penalising chemically dissimilar swaps (e.g. Leu → Pro) far more than conservative ones (e.g. Leu → Ile). Gap penalties still apply and remain fully configurable.

---

## Sample Datasets

| File | Description |
|---|---|
| `dna_sample_1.fasta` | Synthetic reference DNA sequence (244 bp) |
| `dna_sample_2.fasta` | Mutant variant with 3 single-nucleotide substitutions |
| `rna_sample_1.fasta` | Synthetic reference mRNA fragment (RNA alphabet — U, not T) |
| `rna_sample_2.fasta` | Isoform — 2 SNPs + 1 insertion |
| `protein_sample_1.fasta` | Synthetic TP53 DNA-binding domain fragment (reference) |
| `protein_sample_2.fasta` | Oncogenic mutant variant (V → I substitution + length variation) |

All files are valid FASTA and can be loaded instantly via the **Load Seq 1 / Load Seq 2** sidebar buttons.

---

## Exporting Results

### CSV

The mutations CSV contains five columns:

```
Position, Seq1_Base, Seq2_Base, Mutation_Type, Description
```

Compatible with Excel, R (`read.csv()`), and Python (`pd.read_csv()`).

### PDF

The PDF report is structured as follows:

| Section | Content |
|---|---|
| **1. Metadata** | Date, sequence IDs, alignment mode, scoring parameters |
| **2. Similarity Metrics** | Full 11-row metrics table |
| **3. Alignment Preview** | First 1 000 characters of the formatted alignment |
| **4. Alignment Map Chart** | Embedded PNG (requires `kaleido`; fallback note shown otherwise) |
| **5. Mutation Summary** | Counts by mutation type |
| **6. Detailed Mutation Table** | Per-row colour-coded mutation events |

---

## Tech Stack

| Layer | Technology |
|---|---|
| UI / Web framework | [Streamlit](https://streamlit.io/) ≥ 1.35 |
| Sequence I/O & alignment | [Biopython](https://biopython.org/) ≥ 1.83 |
| Data manipulation | [Pandas](https://pandas.pydata.org/) ≥ 2.2 · [NumPy](https://numpy.org/) ≥ 1.26 |
| Interactive charts | [Plotly](https://plotly.com/python/) ≥ 5.22 |
| PDF generation | [ReportLab](https://www.reportlab.com/) ≥ 4.1 |
| Chart image export | [Kaleido](https://github.com/plotly/Kaleido) ≥ 0.2.1 |
| Language | Python 3.11+ |

---

## Known Limitations

| Limitation | Detail |
|---|---|
| **Scale** | `PairwiseAligner` is suitable for gene-scale sequences (up to ~10 kb). Chromosome-scale FASTA files will be very slow. For genome-scale work, consider BLAST or minimap2. |
| **Multi-record FASTA** | Only the first record per file is used for alignment. All records are validated and listed. |
| **Dot plot cap** | Capped at 2 000 bp/aa to prevent browser freezes. Increase the k-mer window to reduce point density for longer sequences. |
| **Kaleido** | If not installed, the PDF alignment chart is replaced by an explanatory note. Install with `pip install kaleido`. |
| **Duplication detection** | Uses a lightweight heuristic (≥ 2 consecutive identical flanking residues). Not a full-featured tandem repeat finder — complex structural variants are out of scope. |

---

## License

This project is released under the [MIT License](LICENSE).

---

*Built with ❤️ using Python · Streamlit · Biopython · Plotly · ReportLab*
