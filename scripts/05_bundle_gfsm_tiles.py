"""Pack gdal2tiles XYZ output into a few bundle files for the web viewer.

Static hosts (and the artifact host's file limits) cope badly with thousands of tiny files,
so tiles are re-encoded as 6-colour palette PNGs and concatenated into bundles:
  * z <= 7           -> one bundle  "lo.png"
  * z >= 8           -> one bundle per (z, parent tile at z6), e.g. "z9_45_27.png" (concatenated PNGs)
An index maps "z/x/y" -> [bundle, offset, length]; the viewer fetches a bundle once and slices.

Output: viewer/data/gfsm/{index.json, *.png}
Usage: python scripts/05_bundle_gfsm_tiles.py
"""
import io
import json
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "interim" / "gfsm_xyz"
OUT = ROOT / "viewer" / "data" / "gfsm"
PARENT_Z = 6
# index 0 = transparent, 1..5 = GFSM classes (same colours as scripts/04_build_gfsm_tiles.sh)
PALETTE = np.array([[0, 0, 0], [44, 123, 182], [171, 217, 233], [255, 255, 191],
                    [253, 174, 97], [215, 25, 28]], dtype=np.int32)


def to_palette_png(path):
    rgba = np.asarray(Image.open(path).convert("RGBA"), dtype=np.int32)
    rgb, alpha = rgba[..., :3], rgba[..., 3]
    d = ((rgb[..., None, :] - PALETTE[None, None, 1:, :]) ** 2).sum(-1)  # nearest class colour
    idx = (d.argmin(-1) + 1).astype(np.uint8)
    idx[alpha < 128] = 0
    if not idx.any():
        return None
    img = Image.fromarray(idx, "P")
    img.putpalette(PALETTE.astype(np.uint8).ravel().tolist())
    buf = io.BytesIO()
    img.save(buf, "PNG", optimize=True, transparency=0)
    return buf.getvalue()


def bundle_name(z, x, y):
    if z <= 7:
        return "lo"
    s = z - PARENT_Z
    return f"z{z}_{x >> s}_{y >> s}"


def main():
    groups = defaultdict(list)
    for png in SRC.glob("*/*/*.png"):
        z, x, y = int(png.parts[-3]), int(png.parts[-2]), int(png.stem)
        groups[bundle_name(z, x, y)].append((z, x, y, png))

    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    names, tiles, total = sorted(groups), {}, 0
    for bi, name in enumerate(names):
        blob = bytearray()
        for z, x, y, png in sorted(groups[name]):
            data = to_palette_png(png)
            if data is None:
                continue
            tiles[f"{z}/{x}/{y}"] = [bi, len(blob), len(data)]
            blob += data
        (OUT / f"{name}.png").write_bytes(blob)
        total += len(blob)
        if len(blob) > 14e6:
            print(f"WARNING: {name}.png is {len(blob) / 1e6:.1f} MB (artifact limit 15 MB)")

    zs = sorted({int(k.split("/")[0]) for k in tiles})
    index = {"minzoom": zs[0], "maxzoom": zs[-1], "bundles": [f"{n}.png" for n in names], "tiles": tiles}
    (OUT / "index.json").write_text(json.dumps(index, separators=(",", ":")))
    print(f"{len(tiles)} tiles in {len(names)} bundles, {total / 1e6:.1f} MB, "
          f"largest {max((OUT / f'{n}.png').stat().st_size for n in names) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
