#!/usr/bin/env python3
"""Check every playback PNG against its exact numerical grain IDs and orientation."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from PIL import Image


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path,
                        default=Path(__file__).resolve().parents[1] / "docs/assets/demo")
    args = parser.parse_args()
    manifest = json.loads((args.data / "evolution.json").read_text())
    count = 0
    with gzip.open(args.data / manifest["numerical_data"], "rt") as stream:
        metadata = json.loads(next(stream))
        nx, nz = metadata["dimensions"]
        palette = np.array([[int(color[i:i + 2], 16) for i in (1, 3, 5)]
                            for color in metadata["palette"]], dtype=np.uint8)
        for line in stream:
            frame = json.loads(line)
            expected_metadata = manifest["frames"][count]
            assert (frame["step"], frame["time_ms"]) == (
                expected_metadata["step"], expected_metadata["time_ms"])
            ids = np.asarray(frame["grain_id"]).reshape(nz, nx)
            # x is fastest; positive z points upward in the displayed image.
            expected = palette[ids[::-1]]
            path = args.data / "evolution" / f"frame-{frame['step']:03d}.png"
            with Image.open(path) as image:
                scale = image.width // nx
                assert scale >= 1 and image.size == (nx * scale, nz * scale)
                actual = np.asarray(image.convert("RGB"))
            expected = np.repeat(np.repeat(expected, scale, axis=0), scale, axis=1)
            assert np.array_equal(actual, expected), f"Pixel mismatch: {path}"
            count += 1
    assert count == len(manifest["frames"])
    print(f"Verified all pixels and orientation of {count} recorded grain-field PNGs.")


if __name__ == "__main__":
    main()
