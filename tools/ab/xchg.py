"""Array exchange between the Python reference and AutoRIFT.jl, with the layout in the file.

Every A/B diagnostic has to move 2-D arrays across the two languages, and the two disagree about
element order: Julia is column-major, NumPy is row-major by default. Writing the raw buffer and
re-deriving the shape at the far end is where that bites, because a square array reshaped with the
wrong convention is silently transposed rather than an error, and a transposed displacement field
still looks like a displacement field. It reads as a spatial offset in whatever is being measured.

So a file written here carries its own header: magic, element type, and both dimensions. `read` uses
the header rather than a caller-supplied shape, and refuses a file it did not write. There is no
argument a caller can get wrong.

    import xchg
    xchg.write("dxf", arr)          # 2-D, any dtype in TYPES
    arr = xchg.read("dxf")          # same array, same orientation

The Julia side is `xchg.jl`, byte-compatible. `selftest()` asserts the round trip in both directions
on a deliberately non-square, non-symmetric array — the only shape that catches a transpose.
"""

import os

import numpy as np

MAGIC = b"ABXC1\0\0\0"
# Type tags, shared with `xchg.jl`. Appending is safe; renumbering is not.
TYPES = {1: np.float32, 2: np.float64, 3: np.int32, 4: np.int64, 5: np.uint8}
TAGS = {v().dtype.str: k for k, v in TYPES.items()}

DIR = os.environ.get("AB_XCHG_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "xchg"))


def path(name):
    return name if os.sep in name else os.path.join(DIR, name + ".abx")


def write(name, A):
    """Write `A` with its shape and element type in the header. Returns the path."""
    A = np.asarray(A)
    if A.ndim != 2:
        raise ValueError("xchg handles 2-D arrays only, got {}".format(A.ndim))
    tag = TAGS.get(A.dtype.str)
    if tag is None:
        raise ValueError("unsupported dtype {}; add it to TYPES in both xchg.py and xchg.jl"
                         .format(A.dtype))
    p = path(name)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "wb") as fh:
        fh.write(MAGIC)
        # Little-endian throughout: both sides run on the same machine, and stating it keeps a
        # cross-endian read an explicit gap rather than a silent scramble.
        fh.write(np.asarray([tag, A.shape[0], A.shape[1]], dtype="<i8").tobytes())
        # Column-major, so Julia's `reshape` over the raw buffer is already the same array.
        fh.write(np.asfortranarray(A).tobytes(order="F"))
    return p


def read(name):
    """The array `write` stored, in the orientation it had. The shape comes from the file."""
    p = path(name)
    with open(p, "rb") as fh:
        if fh.read(len(MAGIC)) != MAGIC:
            raise ValueError("{} was not written by xchg — refusing to guess its layout".format(p))
        tag, nr, nc = np.frombuffer(fh.read(24), dtype="<i8")
        dt = TYPES.get(int(tag))
        if dt is None:
            raise ValueError("unknown type tag {} in {}".format(tag, p))
        buf = fh.read()
    want = nr * nc * np.dtype(dt).itemsize
    if len(buf) != want:
        raise ValueError("{}: header says {}x{} of {} ({} bytes), found {}"
                         .format(p, nr, nc, np.dtype(dt).name, want, len(buf)))
    # `order="F"` undoes the write; the result is the original array, not its transpose.
    return np.ndarray((nr, nc), dtype=dt, buffer=bytearray(buf), order="F")


def selftest(verbose=True):
    """Round trip a non-square, non-symmetric array. Any transpose changes the shape or the values."""
    A = (np.arange(3 * 5, dtype=np.float32) * 1.5).reshape(3, 5)
    A[0, 4] = -7.25                                  # breaks any accidental symmetry
    B = read(write("__selftest", A))
    assert B.shape == A.shape, "shape changed: {} -> {}".format(A.shape, B.shape)
    assert np.array_equal(A, B), "values changed under round trip"
    for dt in (np.float64, np.int32, np.uint8):
        C = A.astype(dt)
        assert np.array_equal(read(write("__selftest", C)), C), "round trip failed for {}".format(dt)
    os.remove(path("__selftest"))
    if verbose:
        print("xchg.py selftest passed")
    return True


if __name__ == "__main__":
    selftest()
