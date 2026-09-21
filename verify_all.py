"""Comprehensive verification script for all 6 review fixes."""
import sys, io, time, random, pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, '.')

print("=== SYNTAX CHECK ===")
import py_compile
for f in ['sequence_reader.py','alignment.py','similarity.py',
          'mutation_detection.py','visualization.py','reports.py','app.py']:
    py_compile.compile(f, doraise=True)
    print(f"  {f} OK")

print("\n=== HAPPY PATH: DNA GLOBAL ===")
from sequence_reader import read_sequence_from_text, read_sequence_from_upload, get_all_records_from_upload
from alignment import run_pairwise_alignment
from similarity import calculate_similarity_metrics
from mutation_detection import detect_mutations, mutation_summary, _flag_duplications, MUTATION_TYPES
from visualization import plot_alignment_map, plot_dot_plot, plot_base_composition, plot_mutation_distribution
from reports import export_mutations_csv, generate_pdf_report

r1 = read_sequence_from_text('ATGCATGCATGCATGC', 'DNA')
r2 = read_sequence_from_text('ATGCATGCATGCATCC', 'DNA')
result = run_pairwise_alignment(r1.sequence, r2.sequence, r1.seq_id, r2.seq_id, mode='global', seq_type='DNA')
metrics = calculate_similarity_metrics(result, 'DNA')
df = detect_mutations(result)
ms = mutation_summary(df)
print(f"  score={result.score} method={result.scoring_method} identity={metrics['percent_identity']}%")
assert result.scoring_method == 'flat', "DNA should use flat scoring"

print("\n=== PROTEIN / BLOSUM62 ===")
p1 = read_sequence_from_text('MKTIIALSYIFCLVFA', 'Protein')
p2 = read_sequence_from_text('MKTIIALSYIFCLVWA', 'Protein')
assert p1.is_valid and p2.is_valid
pres = run_pairwise_alignment(p1.sequence, p2.sequence, mode='global', seq_type='Protein')
print(f"  BLOSUM62 score={pres.score} method={pres.scoring_method}")
assert pres.scoring_method == 'BLOSUM62', "Protein must use BLOSUM62"
# Verify BLOSUM62 gives different score than flat scoring for same sequences
flat_pres = run_pairwise_alignment(p1.sequence, p2.sequence, mode='global', seq_type='DNA')
print(f"  Flat score for same seqs={flat_pres.score}  BLOSUM62={pres.score}")
print(f"  Scores differ: {pres.score != flat_pres.score} (expected True — BLOSUM62 is biologically meaningful)")

print("\n=== EDGE CASE 1: Empty sequence ===")
r_empty = read_sequence_from_text('', 'DNA')
assert not r_empty.is_valid
print(f"  is_valid={r_empty.is_valid}  msg='{r_empty.error_message}'")

print("\n=== EDGE CASE 2: Invalid characters ===")
r_bad = read_sequence_from_text('ATGCXYZ123', 'DNA')
assert not r_bad.is_valid
print(f"  is_valid={r_bad.is_valid}  invalid={r_bad.invalid_chars}")
assert set(r_bad.invalid_chars) == {'1','2','3','X','Y','Z'}

print("\n=== EDGE CASE 3: Multi-header FASTA ===")
multi_fasta = b">seq_a\nATGCATGCATGC\n>seq_b\nTTTTTTTTTTTT\n"
recs = get_all_records_from_upload(multi_fasta, 'DNA')
print(f"  {len(recs)} records: {[r.seq_id for r in recs]}")
assert len(recs) == 2
assert all(r.is_valid for r in recs)

print("\n=== EDGE CASE 4: T in RNA mode (cross-type) ===")
r_cross = read_sequence_from_text('ATGCATGCATGC', 'RNA')
assert not r_cross.is_valid
print(f"  is_valid={r_cross.is_valid}  invalid={r_cross.invalid_chars}")

print("\n=== EDGE CASE 5: Large sequence (~5kb) performance ===")
random.seed(42)
big1 = ''.join(random.choice('ACGT') for _ in range(5000))
big2 = ''.join(random.choice('ACGT') for _ in range(5000))
t0 = time.time()
rbig = run_pairwise_alignment(big1, big2, mode='local', seq_type='DNA')  # local is faster for large seqs
elapsed = time.time() - t0
print(f"  Local 5kb alignment -> score={rbig.score:.1f} time={elapsed:.2f}s")

print("\n=== EDGE CASE 6: Duplication heuristic (>=2 consecutive) ===")
# Poly-A run of 4: should be flagged as Duplication
a1_gap  = 'GCATAAA-GCAT'
a2_ins  = 'GCATAAAAGCAT'
raw_dup = pd.DataFrame([{'Position': 8, 'Seq1_Base': '-', 'Seq2_Base': 'A',
    'Mutation_Type': MUTATION_TYPES['INSERTION'],
    'Description': 'Position 8: Insertion of A in Seq2 (gap in Seq1).'}])
flagged = _flag_duplications(raw_dup, a1_gap, a2_ins, min_repeat_len=2)
assert flagged.iloc[0]['Mutation_Type'] == MUTATION_TYPES['DUPLICATION'], \
    f"Expected Duplication, got {flagged.iloc[0]['Mutation_Type']}"
print(f"  poly-A(4) -> {flagged.iloc[0]['Mutation_Type']} (correct)")
assert 'run of' in flagged.iloc[0]['Description']

# Single flanking base: should NOT be flagged
a1_s  = 'GCATAG-GCAT'
a2_s  = 'GCATAGCGCAT'
raw_s = pd.DataFrame([{'Position': 7, 'Seq1_Base': '-', 'Seq2_Base': 'C',
    'Mutation_Type': MUTATION_TYPES['INSERTION'],
    'Description': 'Position 7: Insertion of C in Seq2 (gap in Seq1).'}])
not_dup = _flag_duplications(raw_s, a1_s, a2_s, min_repeat_len=2)
assert not_dup.iloc[0]['Mutation_Type'] == MUTATION_TYPES['INSERTION'], \
    f"Single-base should stay Insertion"
print(f"  single-flanking -> {not_dup.iloc[0]['Mutation_Type']} (not Duplication — correct)")

print("\n=== EDGE CASE 7: PDF with alignment_fig (kaleido) ===")
fig_map = plot_alignment_map(result)
pdf_b = generate_pdf_report('seq1','seq2','DNA','global',metrics,
    result.formatted_alignment, df, ms, alignment_fig=fig_map)
print(f"  PDF bytes={len(pdf_b)}")
assert len(pdf_b) > 3000

print("\n=== EDGE CASE 8: FASTA file loader ===")
with open('datasets/dna_sample_1.fasta', 'rb') as fh:
    rec = read_sequence_from_upload(fh.read(), 'DNA')
print(f"  dna_sample_1: id={rec.seq_id} valid={rec.valid if hasattr(rec,'valid') else rec.is_valid} len={rec.length}")
assert rec.is_valid

print("\n=== ALL TESTS PASSED ===")
