"""Deterministic scoring. Every point awarded is recorded as a human-readable reason."""
from dataclasses import dataclass, field

# Column aliases: first match (case-insensitive) wins. Edit here if your TSV headers differ.
ALIASES = {
    "gene": ["gene", "gene_symbol", "symbol", "gene_name"],
    "variant": ["variant", "hgvs", "hgvsc", "hgvs_c", "variant_id", "id"],
    "consequence": ["consequence", "effect", "annotation", "variant_effect"],
    "clinvar": ["clinvar", "clinvar_significance", "clin_sig", "clnsig", "clinical_significance"],
    "af": ["gnomad_af", "af", "allele_frequency", "max_af", "gnomad"],
    "cadd": ["cadd", "cadd_phred", "cadd_score"],
    "revel": ["revel", "revel_score"],
    "hpo": ["matched_hpo", "hpo_terms", "hpo", "hpo_match", "hpo_ids"],
    "zygosity": ["zygosity", "genotype", "gt"],
}

CONSEQ_POINTS = {
    "stop_gained": 20, "frameshift": 20, "splice_donor": 20, "splice_acceptor": 20,
    "start_lost": 18, "stop_lost": 16, "inframe": 10, "missense": 8, "splice_region": 6,
    "synonymous": 1, "intron": 0, "intergenic": 0, "utr": 1,
}
CLINVAR_POINTS = [  # order matters: check "likely" forms before bare forms
    ("conflicting", 12), ("likely_pathogenic", 32), ("likely pathogenic", 32),
    ("pathogenic", 40), ("uncertain", 10), ("vus", 10),
    ("likely_benign", -15), ("likely benign", -15), ("benign", -25),
]
TIERS = [(75, "Very High"), (55, "High"), (35, "Medium"), (15, "Low"), (-999, "Very Low")]


@dataclass
class Scored:
    gene: str
    variant: str
    score: int
    tier: str
    reasons: list = field(default_factory=list)
    source: str = ""


def _num(x):
    try:
        return float(str(x).strip())
    except (ValueError, TypeError):
        return None


def pick(row, key):
    lowered = {k.lower().strip(): v for k, v in row.items()}
    for a in ALIASES[key]:
        if a in lowered and str(lowered[a]).strip() not in ("", ".", "NA", "nan"):
            return str(lowered[a]).strip()
    return ""


def tier_for(score):
    for cut, name in TIERS:
        if score >= cut:
            return name


def score_row(row, source="matched_hpo", patient_hpo=None):
    reasons, s = [], 0
    cv = pick(row, "clinvar").lower()
    for key, pts in CLINVAR_POINTS:
        if key in cv:
            s += pts
            reasons.append(f"ClinVar '{cv}' -> {pts:+d}")
            break
    else:
        reasons.append("No ClinVar assertion -> +0")

    conseq = pick(row, "consequence").lower()
    best = max((p for k, p in CONSEQ_POINTS.items() if k in conseq), default=0)
    s += best
    reasons.append(f"Consequence '{conseq or 'unknown'}' -> {best:+d}")

    af = _num(pick(row, "af"))
    if af is None:
        pts, why = 5, "frequency unknown (treated as possibly rare)"
    elif af < 0.0001: pts, why = 15, f"AF {af:g} very rare"
    elif af < 0.001: pts, why = 10, f"AF {af:g} rare"
    elif af < 0.01: pts, why = 4, f"AF {af:g} low-frequency"
    elif af < 0.05: pts, why = -10, f"AF {af:g} common-ish"
    else: pts, why = -25, f"AF {af:g} common"
    s += pts
    reasons.append(f"{why} -> {pts:+d}")

    cadd, revel = _num(pick(row, "cadd")), _num(pick(row, "revel"))
    ins = 0
    if cadd is not None and cadd >= 25: ins += 5
    elif cadd is not None and cadd >= 15: ins += 2
    if revel is not None and revel >= 0.75: ins += 5
    elif revel is not None and revel >= 0.5: ins += 2
    s += ins
    reasons.append(f"In-silico (CADD={cadd}, REVEL={revel}) -> {ins:+d}")

    if source == "matched_hpo":
        terms = {t.strip() for t in pick(row, "hpo").replace(";", ",").replace("|", ",").split(",") if t.strip()}
        overlap = len(terms & set(patient_hpo)) if patient_hpo else len(terms)
        pts = min(15, 5 * overlap) if overlap else 3
        s += pts
        reasons.append(f"Phenotype: {overlap} matched HPO term(s) -> {pts:+d}")
    elif source == "unmatched_hpo":
        reasons.append("Gene has no HPO match to patient phenotype -> +0")
    elif source == "secondary_findings":
        s += 10
        reasons.append("Secondary finding list (reportable regardless of phenotype) -> +10")

    return Scored(pick(row, "gene") or "?", pick(row, "variant") or "?", s, tier_for(s), reasons, source)
