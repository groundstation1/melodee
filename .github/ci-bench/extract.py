"""Pull the ROOM / REVERB / HQ REVERB kernels out of firmware/src/fx.c into a standalone include."""
import sys
src = open(sys.argv[1], encoding="utf-8").read()
def cut(a, b):
    i = src.index(a); j = src.index(b, i); return src[i:j]
out = [cut("#define FX_RATE", "static void fx_coeff"),
       cut("/* Fuse interpolation/hold", "/* Full-band stereo chorus"),
       cut("/* Eight-line orthogonal", "/* The filter inserts use")]
open(sys.argv[2], "w", encoding="utf-8").write("\n".join(out))
