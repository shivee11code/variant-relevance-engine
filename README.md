# Variant Relevance Engine (VRE) - prototype

## The one flow
Annotated germline WES TSVs in -> every variant gets a tier (Very High / High / Medium / Low / Very Low) plus a line-by-line reason for every point awarded -> ranked list out.

## Setup (needs only Python 3.8+, no packages, no API keys)
```
cd prototype
python -m vre sample_data/*.tsv                 # ranked, with reasons
python -m vre sample_data/*.tsv --json          # machine-readable
python -m vre my_matched_hpo.tsv --hpo HP:0001250,HP:0001263
python -m unittest discover -s tests            # 7 tests, should print OK
```
File names decide the source type: a name containing `unmatched` or `secondary` is treated as unmatched-HPO or secondary-findings; anything else is matched-HPO.
If your TSV headers differ, edit the `ALIASES` table at the top of `vre/scoring.py`.

## What works
- Deterministic scoring from ClinVar, consequence, population frequency, CADD/REVEL, phenotype match.
- Explainable output: every score lists its reasons.
- Tolerates missing columns and values.

## What does NOT work / is not built (stated deliberately)
- **No LLM layer.** Plain-language explanations are template-based only.
- **No HPO semantic similarity.** Phenotype match is a simple term-overlap count, not ontology-aware.
- **Weights and tier cut-offs are my assumptions, not validated** against a labelled truth set.
- **Only tested on synthetic demo rows** in `sample_data/` (not real patient data). Column names are assumed; adjust `ALIASES`.
- No hosted link or demo recording yet.
- Not for clinical use.
