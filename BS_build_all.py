"""Master build: train all endpoints -> external (cluster + conformal) -> temporal + PR.
Runs unattended; each stage prints a banner. Outputs land in models_brain/ and *_report.json."""
import time, sys
t0 = time.time()
def banner(s): print("\n" + "=" * 70 + f"\n=== {s}  (+{time.time()-t0:.0f}s)\n" + "=" * 70, flush=True)

banner("STAGE 1/3  TRAIN + CALIBRATE + EVIDENCE (all endpoints)")
import BS_train_endpoints; BS_train_endpoints.main()

banner("STAGE 2/3  EXTERNAL VALIDATION (leave-cluster-out + conformal)")
import BS_external_validation; BS_external_validation.main()

banner("STAGE 3/3  TEMPORAL + PR / THRESHOLD")
import BS_temporal_pr; BS_temporal_pr.main()

banner(f"BUILD COMPLETE in {time.time()-t0:.0f}s")
