"""Build the Supplementary Information document that accompanies the condensed manuscript.

The condensed draft (NAR's page-limited Web Server Issue submission) embeds two figures and refers
the reader to "Supplementary Figure S1/S2/S3" for three more. Those three previously existed only as
a bare file path typed into the running prose (`figures/Figure9_model_atlas.png`), which is not how
a submitted manuscript cites supplementary material and is not something a reviewer can open from the
PDF. This assembles the actual file: the three figures, full-size, each with the same caption the
full manuscript uses for them there, as one standalone document submitted alongside the manuscript.

Run:  python src/brainsafe/analysis/build_supplementary_information.py
"""
from __future__ import annotations

from pathlib import Path

import pypandoc

ROOT = Path(__file__).resolve().parents[3]
MS = ROOT / "manuscript"
REFERENCE_DOCX = ROOT / "docs" / "nar_manuscript_reference.docx"
OUT_MD = MS / "Supplementary_Information_built.md"
OUT_DOCX = MS / "Supplementary_Information.docx"

TITLE = ("# Supplementary Information\n\n"
         "**BrainSafe AI: a calibrated, exposure-gated web server for multi-endpoint prediction of "
         "small-molecule action in the human brain**\n\n"
         "Krishnasalini Gunanathan, Raghunatha Sarma, Sai Shyam, Ramya E. M., "
         "Venketesh Sivaramakrishnan\n\n"
         "This document contains the three validation figures referred to in the main text as "
         "Supplementary Figures S1 to S3, and Supplementary Table S1. All three figures, and the "
         "table, are generated directly from the repository's own results tables and model registry "
         "(https://github.com/krishna-g-999/brainsafe-ai) and reproduce with the commands named "
         "beside each.\n\n---\n\n")

FIGURES = [
    ("S1", "Figure9_model_atlas.png",
     "**Supplementary Figure S1.** The panel, one mark per estimator, so that no claim rests on a "
     "mean a reader cannot check. (**A**) Every estimator, deployed or withdrawn, placed by the "
     "number of compounds it was trained on and by the performance claimed for it, coloured by "
     "model family. Marker shape carries the metric: AUROC and R² both run to 1.0 and are not "
     "the same quantity, since 0.5 is chance for one and a respectable fit for the other, so they "
     "are distinguished rather than averaged. The five estimators withdrawn after specificity "
     "testing are drawn in outline, because a panel showing only what survived is a selection "
     "rather than an inventory. Training sets span two orders of magnitude, from 37 to 15,831 "
     "rows; binder training sets include property-matched decoys while the others are measured "
     "compounds only, which the axis states. (**B**) The same population by family, with the "
     "family median marked. The spread is the point: the binder classifiers have a median of "
     "0.945 while the exposure and ADME family, which mixes AUROC and R² endpoints, has a "
     "median of 0.574, and a single panel average would describe neither. Reproduces with "
     "`python src/brainsafe/figures/fig09_model_atlas.py`.\n\n"
     "Supplementary Table S1, the complete inventory behind this figure (training-set composition, "
     "validation scheme, calibration and fitting date for every one of the 75 estimators), is "
     "`results/tables/MODEL_INVENTORY.csv`, submitted as a separate CSV file alongside this "
     "document."),
    ("S2", "Figure6_validation.png",
     "**Supplementary Figure S2.** Four validations that a cross-validated score cannot replace. "
     "(**A**) Expected calibration error before and after isotonic regression fitted on out-of-fold "
     "predictions, so no compound contributes to the calibrator that scores it. (**B**) Recall on "
     "whole scaffold classes withheld before training, with 95 per cent Wilson intervals and marker "
     "area proportional to the number of withheld actives, so an interval that is wide because the "
     "evidence is thin looks thin. (**C**) Specificity on chemistry the server should stay quiet "
     "about, and external discrimination on approved drugs absent from the training source. "
     "(**D**) The adversarial suite, in which each check was written so that it could fail. All six "
     "pass, one of them, the applicability-domain flag, only after its controls were corrected "
     "rather than its criterion retuned; every check is drawn at the same size whatever its "
     "verdict, so a future failure would be exactly this visible. Reproduces with "
     "`python src/brainsafe/figures/fig06_validation.py`."),
    ("S3", "Figure11_external_validation.png",
     "**Supplementary Figure S3.** External validation, and an apparent temporal decay that is a "
     "composition effect. (**A**) Per endpoint, AUROC under a size-matched random split against "
     "AUROC under a time split that withholds the most recent quarter of the data and freezes the "
     "decision threshold before the cutoff. The size match matters: a time split trains on less "
     "data as well as none of the future, so without a control at the same n a drop cannot be "
     "attributed to either. (**B**) The same comparison for sensitivity at the frozen operating "
     "point, where the gap is larger. (**C**) Why the gap exists. A random split of "
     "medicinal-chemistry data holds out mostly close analogues of its own training set, because "
     "the published record is series; a time split does not. The two are not testing comparable "
     "populations. (**D**) The resolution. Recall against maximum Tanimoto similarity to the "
     "training actives, for three test sets built by unrelated rules: withheld by publication "
     "date, withheld at random, and withheld by curator, the last being compounds deposited in "
     "BindingDB and absent from ChEMBL. They trace one curve, so recall is a function of chemical "
     "distance rather than of how the set was held out, and the expected sensitivity for a "
     "submitted compound is knowable at query time. Reproduces with "
     "`python src/brainsafe/figures/fig11_external_validation.py`."),
]


def main() -> None:
    parts = [TITLE]
    for label, fname, caption in FIGURES:
        parts.append(f"## Supplementary Figure {label}\n\n")
        parts.append(f"![Supplementary Figure {label}](figures/{fname})\n\n")
        parts.append(caption + "\n\n---\n\n")
    text = "".join(parts)
    OUT_MD.write_text(text, encoding="utf-8")

    pypandoc.convert_file(str(OUT_MD), "docx", outputfile=str(OUT_DOCX),
                          extra_args=[f"--resource-path={MS}",
                                      f"--reference-doc={REFERENCE_DOCX}"])
    print(f"wrote {OUT_MD.relative_to(ROOT)} ({len(text):,} chars)")
    print(f"wrote {OUT_DOCX.relative_to(ROOT)} ({OUT_DOCX.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
