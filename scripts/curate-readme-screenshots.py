#!/usr/bin/env python3
"""Curate real browser-test screenshots for the project README.

Input: screenshots in benchmarks/output/ui-preview from a successful UI run.
Output: cropped, optimised PNGs in docs/screenshots.

Review the resulting images before committing them after interface changes.
Prediction examples may use synthetic fixtures; captions must disclose this.
"""
from pathlib import Path
from PIL import Image

SOURCE = Path("benchmarks/output/ui-preview")
DEST = Path("docs/screenshots")

# Coordinates refer to the tested 1440x940 / 1366x900 full-page captures.
# Crop the relevant analysis panel, not the compact structure-input toolbar.
EXAMPLES = {
    "overview.png": ("desktop-loaded.png", (355, 232, 1405, 1205)),
    "residue-inspection.png": ("desktop-residue-zoom.png", (5, 165, 1420, 1120)),
    "all-chains.png": ("all-chains-expanded.png", (380, 1030, 1040, 1690)),
    "prediction-pae.png": ("prediction-af2-pae.png", (355, 1390, 1300, 2180)),
    "ensemble.png": ("phase-c-ensemble.png", (380, 812, 1340, 1570)),
}

def main():
    DEST.mkdir(parents=True, exist_ok=True)
    for name, (original, bounds) in EXAMPLES.items():
        path = SOURCE / original
        if not path.is_file():
            raise FileNotFoundError(f"Missing {path}; rerun the real browser tests.")
        with Image.open(path) as picture:
            if picture.width < bounds[2] or picture.height < bounds[3]:
                raise ValueError(f"{original}: unexpected {picture.size}; update crop after reviewing UI.")
            crop = picture.convert("RGB").crop(bounds)
            crop.thumbnail((1150, 1500), Image.Resampling.LANCZOS)
            crop.save(DEST / name, format="PNG", optimize=True)
            print(f"{name}: {crop.width}x{crop.height}")

if __name__ == "__main__":
    main()
