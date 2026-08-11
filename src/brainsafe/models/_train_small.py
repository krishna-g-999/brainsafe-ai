"""Train the small measured-label binder models: TAAR1, GluA2 and Nav1.1.

This was a verbatim copy of train_measured_label_holdout.py differing only in its target list, and
carried its own copy of the same leak, so a correction applied to one file left the other reporting
sensitivity on training positives. It is kept as an entry point because it is the only recorded way
the TAAR1, GluA2 and Nav1_1 models were produced, but it now delegates.

Note that Nav1_1 is withdrawn from the deployed panel (see apply_specificity_decisions.py); it is
retained here so the historical run can be reproduced, not because it is served.

Equivalent to:  python src/brainsafe/models/train_measured_label_holdout.py TAAR1 GluA2 Nav1_1
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import train_measured_label_holdout  # noqa: E402

TARGETS = ["TAAR1", "GluA2", "Nav1_1"]

if __name__ == "__main__":
    train_measured_label_holdout.TARGETS = TARGETS
    train_measured_label_holdout.main()
