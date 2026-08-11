"""One partition of the background chemistry, so that no two roles can draw the same compound.

The binder pipeline asks the background library for three different things, and until now took all
three from the same place:

  decoy       compounds trained on as presumed-negative, to teach the model what a non-binder looks
              like when measured inactives are scarce
  threshold   compounds whose score distribution sets the decision threshold, as a quantile
  evaluation  compounds used to report the false-positive rate at that threshold

Taking the threshold and the evaluation sample from one pool makes the reported rate an arithmetic
consequence of the quantile rather than a measurement: the threshold is the 95th percentile of the
sample, so the false-positive rate on that same sample is 5 per cent whatever the model does. Taking
decoys from the same pool as either is worse, because the model was explicitly trained to score
those compounds as zero, and it is then congratulated for doing so.

The split is by a stable hash of the structure, not by shuffling. Two properties follow that a
shuffle does not give: the assignment does not depend on the order the library happens to be in, and
adding compounds later leaves every existing assignment untouched, so a threshold set today remains
comparable with a false-positive rate measured next year.

`external_negatives` answers a fourth need. The applicability-domain reference is built from the
training tables, so any sample drawn from it has a nearest neighbour at Tanimoto 1.000 and cannot
support a statement about novel chemistry. The DrugBank structures that are absent from that
reference can.
"""
from __future__ import annotations

import hashlib
import pickle
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
M = ROOT / "models_rf"
AD_REFERENCE = M / "ad_reference.pkl"
EXTERNAL_DRUGS = ROOT / "data" / "external" / "processed" / "external_drugs.csv"

# Shares out of 100. Decoys need the most, being 3x the actives for every target; the other two need
# only enough for a stable quantile and a stable rate.
SHARES = {"decoy": 60, "threshold": 20, "evaluation": 20}
ROLES = tuple(SHARES)


def _band(smiles: str, salt: str = "brainsafe-background-v1") -> int:
    """Stable 0-99 band for a structure, independent of library order and of insertion time."""
    digest = hashlib.blake2b(f"{salt}:{smiles}".encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "big") % 100


def role_of(smiles: str) -> str:
    """Which pool a structure belongs to. Total function: every structure has exactly one role."""
    b = _band(smiles)
    edge = 0
    for role, share in SHARES.items():
        edge += share
        if b < edge:
            return role
    return ROLES[-1]


def background_pools(with_fingerprints: bool = False):
    """Split the applicability-domain reference into the three disjoint role pools.

    Returns {role: [smiles]}, or {role: ([smiles], [fingerprint])} when fingerprints are wanted,
    preserving the pairing the reference file stores.
    """
    if not AD_REFERENCE.exists():
        raise FileNotFoundError(
            f"{AD_REFERENCE.relative_to(ROOT).as_posix()} is missing; build it with "
            "src/brainsafe/build_ad_reference.py before training binders."
        )
    with AD_REFERENCE.open("rb") as fh:
        smiles, fps = pickle.load(fh)

    pools: dict[str, list] = {r: [] for r in ROLES}
    fp_pools: dict[str, list] = {r: [] for r in ROLES}
    for s, fp in zip(smiles, fps):
        r = role_of(str(s))
        pools[r].append(str(s))
        fp_pools[r].append(fp)
    if with_fingerprints:
        return {r: (pools[r], fp_pools[r]) for r in ROLES}
    return pools


def external_negatives() -> list[str]:
    """Approved and experimental drugs that are absent from the applicability-domain reference.

    Used where a specificity claim has to be about chemistry the model has not seen. A sample drawn
    from the reference itself cannot support one: every compound in it is its own nearest neighbour.
    """
    if not EXTERNAL_DRUGS.exists():
        raise FileNotFoundError(
            f"{EXTERNAL_DRUGS.relative_to(ROOT).as_posix()} is missing; build it with "
            "src/brainsafe/data/integrate_external.py."
        )
    with AD_REFERENCE.open("rb") as fh:
        ref = set(map(str, pickle.load(fh)[0]))
    drugs = pd.read_csv(EXTERNAL_DRUGS)
    col = "canonical_smiles" if "canonical_smiles" in drugs.columns else drugs.columns[0]
    return sorted({s for s in drugs[col].astype(str) if s not in ref})


def summary() -> pd.DataFrame:
    """Sizes of every pool, for the run log and for the record."""
    pools = background_pools()
    rows = [{"pool": r, "n": len(pools[r]), "share_pct": SHARES[r]} for r in ROLES]
    rows.append({"pool": "external_negatives", "n": len(external_negatives()), "share_pct": None})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    print(summary().to_string(index=False))
    pools = background_pools()
    overlap = set(pools["decoy"]) & set(pools["threshold"]) | \
              set(pools["decoy"]) & set(pools["evaluation"]) | \
              set(pools["threshold"]) & set(pools["evaluation"])
    print(f"\npairwise overlap between pools: {len(overlap)} (must be 0)")
