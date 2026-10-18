"""Keep historical catalog builders explicit; v2 has a separate exporter."""
import sys


def require_legacy_export():
    if "--legacy" not in sys.argv:
        raise SystemExit("This is a historical v1 exporter. Use build_candidate_catalog.py --config config/step2.yaml "
                         "for corrected v2 candidates. Pass --legacy only to regenerate historical v1 outputs.")
