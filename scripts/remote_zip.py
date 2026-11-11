"""Read single members out of large remote zip archives with HTTP range requests (zip and zip64).

Used by scripts/12_sea_level.py for the IPCC AR6 sea-level archive (9 GB) and DeltaDTM Asia (17 GB),
of which only a few hundred MB are needed.
"""
import struct
import zlib
from pathlib import Path

import requests


def _get(url, a, b, **kw):
    r = requests.get(url, headers={"Range": f"bytes={a}-{b}"}, timeout=900, **kw)
    r.raise_for_status()
    if r.status_code != 206:
        raise RuntimeError(f"{url}: server ignored the range request")
    return r


def index(url):
    """name -> (local header offset, compressed size, method)."""
    size = int(_get(url, 0, 0).headers["Content-Range"].split("/")[1])
    tail = _get(url, max(0, size - 65536), size - 1).content
    i = tail.rfind(b"PK\x06\x06")  # zip64 end of central directory
    if i >= 0:
        cdsize, cdoff = struct.unpack("<QQ", tail[i + 40:i + 56])
    else:
        i = tail.rfind(b"PK\x05\x06")
        cdsize, cdoff = struct.unpack("<II", tail[i + 12:i + 20])
    cd = _get(url, cdoff, cdoff + cdsize - 1).content
    out, p = {}, 0
    while p < len(cd):
        f = struct.unpack("<IHHHHHHIIIHHHHHII", cd[p:p + 46])
        name = cd[p + 46:p + 46 + f[10]].decode()
        csize, lho = f[8], f[16]
        ex, q = cd[p + 46 + f[10]:p + 46 + f[10] + f[11]], 0
        while q < len(ex):  # zip64 extra field carries the 64-bit sizes / offset
            tag, ln = struct.unpack("<HH", ex[q:q + 4])
            if tag == 1:
                vals, k = struct.unpack("<" + "Q" * (ln // 8), ex[q + 4:q + 4 + ln]), 0
                if f[9] == 0xFFFFFFFF:
                    k += 1
                if f[8] == 0xFFFFFFFF:
                    csize = vals[k]; k += 1
                if f[16] == 0xFFFFFFFF:
                    lho = vals[k]
            q += 4 + ln
        out[name] = (lho, csize, f[4])
        p += 46 + f[10] + f[11] + f[12]
    return out


def extract(url, entry, dest):
    """Write one member (stored or deflated) to dest."""
    lho, csize, method = entry
    h = _get(url, lho, lho + 29).content
    start = lho + 30 + sum(struct.unpack("<HH", h[26:30]))
    dest = Path(dest)
    tmp = dest.with_name(dest.name + ".part")
    dec = zlib.decompressobj(-15) if method == 8 else None
    step = 16 << 20  # fetch in 16 MB ranges, retrying each, so a dropped connection costs one chunk
    with open(tmp, "wb") as f:
        for a in range(start, start + csize, step):
            b = min(a + step, start + csize) - 1
            for attempt in range(6):
                try:
                    ch = _get(url, a, b).content
                    if len(ch) == b - a + 1:
                        break
                except requests.RequestException:
                    pass
                if attempt == 5:
                    raise RuntimeError(f"{url}: bytes {a}-{b} failed")
            f.write(dec.decompress(ch) if dec else ch)
        if dec:
            f.write(dec.flush())
    tmp.replace(dest)
    return dest
