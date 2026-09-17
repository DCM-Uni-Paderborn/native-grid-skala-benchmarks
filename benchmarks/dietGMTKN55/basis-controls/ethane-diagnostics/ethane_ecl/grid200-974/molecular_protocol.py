"""Validate only the explicitly declared deviations from the paired protocol."""
import re
from baseline_protocol import validate_input as baseline_validate

def validate_input(text, case):
    for key, n in zip(("RADIAL_GRID", "LEBEDEV_GRID"), case["grid"]):
        assert re.findall(r"(?m)^\s*" + key + r"\s+(\S+)", text) == [str(n)] * case["kinds"]
    normalized = text.replace("RADIAL_GRID 200", "RADIAL_GRID 100").replace("LEBEDEV_GRID 974", "LEBEDEV_GRID 434")
    if case.get("seed"):
        assert "SCF_GUESS RESTART" in normalized and "WFN_RESTART_FILE_NAME source.wfn" in normalized
        normalized = normalized.replace("SCF_GUESS RESTART", "SCF_GUESS ATOMIC")
        normalized = re.sub(r"(?m)^\s*WFN_RESTART_FILE_NAME source.wfn\n", "", normalized)
    return baseline_validate(normalized, case)
