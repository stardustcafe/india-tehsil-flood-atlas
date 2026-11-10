"""Population density tiles for the viewer, from the 30" grids written by 06_population_exposure.py.

For each source found in data/processed/ ({source}2020_india_pop.tif and _pop_high.tif):
  people per km2 (true cell area) -> 6 density classes -> EPSG:3857 at the z8 pixel size (nearest,
  the grid is ~1 km so z8 already shows single cells) -> colour relief -> gdal2tiles z4-8
  -> palette PNG bundles via scripts/05_bundle_gfsm_tiles.py.

Output: viewer/data/pop_{source}_{all,high}/{index.json, *.png}
Usage: python scripts/07_build_pop_tiles.py
"""
import importlib.util
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import rasterio

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
WORK = ROOT / "data" / "interim" / "pop_tiles"
VIEW = ROOT / "viewer" / "data"
CONDA = Path(os.environ.get("CONDA_ROOT", r"C:\Users\Admin\miniconda3"))
MAXZ = 8
# people per km2; class k covers [BREAKS[k-1], BREAKS[k]), class 6 is >= 5000
BREAKS = [1, 10, 50, 250, 1000, 5000]
COLORS = ["e9e3f2", "c5b3dc", "9d84c4", "7652a8", "52307f", "2a1450"]  # light -> dark purple

spec = importlib.util.spec_from_file_location("bundle", ROOT / "scripts" / "05_bundle_gfsm_tiles.py")
bundler = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bundler)


def gdal_env():
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([str(CONDA / "Library" / "bin"), str(CONDA), str(CONDA / "Scripts"), env["PATH"]])
    env["GDAL_DATA"] = str(CONDA / "Library" / "share" / "gdal")
    env["PROJ_DATA"] = str(CONDA / "Library" / "share" / "proj")
    return env


def classify(src_tif, dst_tif):
    with rasterio.open(src_tif) as src:
        pop = src.read(1)
        profile = src.profile
        t = src.transform
    lat = t.f + t.e * (np.arange(pop.shape[0]) + 0.5)
    # spherical cell area in km2 (authalic radius); within 0.3% of ellipsoidal at India's latitudes
    R, d = 6371.0072, np.radians(abs(t.a))
    area = R ** 2 * d * np.abs(np.sin(np.radians(lat + t.e / 2)) - np.sin(np.radians(lat - t.e / 2)))
    dens = pop / area[:, None]
    cls = np.digitize(dens, BREAKS).astype(np.uint8)  # 0 below 1/km2, 1..6
    cls[pop < 0] = 0
    profile.update(dtype="uint8", nodata=0, compress="deflate")
    with rasterio.open(dst_tif, "w", **profile) as dst:
        dst.write(cls, 1)
    return dens[pop >= 0]


def main():
    env, res = gdal_env(), 40075016.685578488 / (256 * 2 ** MAXZ)
    WORK.mkdir(parents=True, exist_ok=True)
    colors = WORK / "pop_colors.txt"
    colors.write_text("nv 0 0 0 0\n" + "".join(
        f"{k} {int(c[:2], 16)} {int(c[2:4], 16)} {int(c[4:], 16)} 255\n" for k, c in enumerate(COLORS, 1)))
    palette = bundler.parse_palette(",".join(COLORS))

    for pop_tif in sorted(PROC.glob("*2020_india_pop.tif")):
        source = pop_tif.name.split("2020_")[0]
        for kind, tif in [("all", pop_tif), ("high", pop_tif.with_name(f"{source}2020_india_pop_high.tif"))]:
            name = f"pop_{source}_{kind}"
            cls = WORK / f"{name}_cls.tif"
            dens = classify(tif, cls)
            share = np.bincount(np.digitize(dens, BREAKS), minlength=7)[1:] / max(1, (dens >= 1).sum())
            print(f"[{name}] cells by class: " + " ".join(f"{s:.0%}" for s in share))
            merc, rgba, xyz = WORK / f"{name}_3857.tif", WORK / f"{name}_rgba.vrt", WORK / f"{name}_xyz"
            subprocess.run(["gdalwarp", "-q", "-overwrite", "-t_srs", "EPSG:3857", "-tr", str(res), str(res), "-tap",
                            "-r", "near", "-co", "COMPRESS=DEFLATE", str(cls), str(merc)], check=True, env=env)
            subprocess.run(["gdaldem", "color-relief", "-q", "-alpha", "-exact_color_entry", "-of", "VRT",
                            str(merc), str(colors), str(rgba)], check=True, env=env)
            shutil.rmtree(xyz, ignore_errors=True)
            subprocess.run([str(CONDA / "python.exe"), str(CONDA / "Scripts" / "gdal2tiles.py"), "--xyz", "-z", f"4-{MAXZ}",
                            "-r", "near", "-w", "none", "-x", "-q", f"--processes={max(1, (os.cpu_count() or 2) // 2)}",
                            str(rgba), str(xyz)], check=True, env=env)
            bundler.bundle(xyz, VIEW / name, palette)


if __name__ == "__main__":
    main()
