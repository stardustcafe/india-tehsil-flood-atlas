#!/usr/bin/env bash
# Build XYZ PNG tiles of the GFSM v1 classes over India with gdal2tiles.
#
#   tiles_india/*.tif (EPSG:3395, 30 m) -> mosaic VRT
#   -> one gdalwarp to EPSG:3857 at the max-zoom pixel size, -r mode (categorical-safe),
#      aligned to the XYZ tile grid and clipped to the India ADM3 outline
#   -> GFSM colour palette (RGBA VRT) -> gdal2tiles --xyz
#
# Usage: bash scripts/04_build_gfsm_tiles.sh [MAXZOOM] [PROCESSES]     (default 10, nproc-1)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
command -v cygpath >/dev/null && ROOT="$(cygpath -m "$ROOT")"  # Windows GDAL needs C:/ paths
CONDA="${CONDA_ROOT:-/c/Users/Admin/miniconda3}"
export PATH="$CONDA/Library/bin:$CONDA:$CONDA/Scripts:$PATH"
export GDAL_DATA="$CONDA/Library/share/gdal" PROJ_DATA="$CONDA/Library/share/proj"
export GDAL_CACHEMAX=2048

MAXZ="${1:-10}"
PROCS="${2:-$(( $(nproc) - 1 ))}"
WORK="$ROOT/data/interim/gfsm_tiles"
OUT="$ROOT/data/interim/gfsm_xyz"
mkdir -p "$WORK"
# web-mercator pixel size at MAXZ: equatorial circumference / (256 * 2^z)
RES=$(python -c "print(40075016.685578488 / (256 * 2 ** $MAXZ))")

# India outline from the same ADM3 boundaries used for the zonal statistics
"$ROOT/.venv/Scripts/python" -I -c "
import geopandas as gpd, sys
g = gpd.read_file(sys.argv[1])
gpd.GeoDataFrame(geometry=[g.union_all().simplify(0.002)], crs=g.crs).to_file(sys.argv[2], driver='GPKG')
" "$ROOT/data/raw/boundaries/geoboundaries/geoBoundaries-IND-ADM3.geojson" "$WORK/india_mask.gpkg"

ls "$ROOT"/data/raw/gfsm/tiles_india/*.tif > "$WORK/tiles.txt"
gdalbuildvrt -q -overwrite -srcnodata 0 -vrtnodata 0 -input_file_list "$WORK/tiles.txt" "$WORK/gfsm_3395.vrt"

echo "[warp] mode-resampling to ${RES} m (z${MAXZ}) in EPSG:3857"
gdalwarp -q -overwrite -t_srs EPSG:3857 -tr "$RES" "$RES" -tap -r mode \
  -cutline "$WORK/india_mask.gpkg" -crop_to_cutline -srcnodata 0 -dstnodata 0 -ot Byte \
  -multi -wo NUM_THREADS=ALL_CPUS -wm 2048 \
  -co COMPRESS=DEFLATE -co TILED=YES -co BIGTIFF=IF_SAFER \
  "$WORK/gfsm_3395.vrt" "$WORK/gfsm_3857_z${MAXZ}.tif"

# GFSM class colours (NoData -> transparent)
cat > "$WORK/gfsm_colors.txt" <<'EOF'
nv 0 0 0 0
1 44 123 182 255
2 171 217 233 255
3 255 255 191 255
4 253 174 97 255
5 215 25 28 255
EOF
gdaldem color-relief -q -alpha -exact_color_entry -of VRT \
  "$WORK/gfsm_3857_z${MAXZ}.tif" "$WORK/gfsm_colors.txt" "$WORK/gfsm_rgba.vrt"

echo "[tiles] gdal2tiles z4-${MAXZ} on ${PROCS} processes"
rm -rf "$OUT"
python "$CONDA/Scripts/gdal2tiles.py" --xyz -z "4-${MAXZ}" -r near -w none -x -q \
  --processes="$PROCS" "$WORK/gfsm_rgba.vrt" "$OUT"

echo "tiles: $(find "$OUT" -name '*.png' | wc -l), size: $(du -sh "$OUT" | cut -f1)"
