from __future__ import annotations
import os
import sys
import argparse
import re 
import glob
from typing import Tuple, Optional
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
import imageio


def jet2_colormap(m: int = 64) -> np.ndarray:
    """
    Reproduce MATLAB jet2 colormap. Returns an (m,3) array of floats in [0,1].
    """

    m_requested = m
    m = 60
    n = max(round(m / 4), 1)
    n = int(n)
    x = (np.arange(1, n + 1) / n).reshape(-1, 1)
    y = (np.arange(int(n / 2), n + 1) / n).reshape(-1, 1)
    e = np.ones((len(x), 1))

    r = np.vstack([np.zeros_like(y), np.zeros_like(e), x, e, np.flipud(y)])
    g = np.vstack([np.zeros_like(y), x, e, np.flipud(x), np.zeros_like(y)])
    b = np.vstack([y, e, np.flipud(x), np.zeros_like(e), np.zeros_like(y)])

    J2 = np.hstack([r, g, b])

    # trim to m rows
    while J2.shape[0] > m:
        J2 = J2[1:, :]
        if J2.shape[0] > m:
            J2 = J2[:-1, :]

    # pad to at least 64 rows to allow specific assignments
    if J2.shape[0] < 64:
        pad_rows = 64 - J2.shape[0]
        J2 = np.vstack([J2, np.tile(J2[-1, :], (pad_rows, 1))])

    # change the last colors per MATLAB script (0-based indices)
    J2[60, :] = [0.0, 0.0, 0.0]
    J2[61, :] = [1.0, 0.0, 0.0]
    J2[62, :] = [0.0, 1.0, 0.0]
    J2[63, :] = [1.0, 1.0, 1.0]

    # resample/interpolate if user requested different length
    if m_requested != J2.shape[0]:
        old_idx = np.linspace(0, 1, J2.shape[0])
        new_idx = np.linspace(0, 1, m_requested)
        J2_resampled = np.zeros((m_requested, 3))
        for c in range(3):
            J2_resampled[:, c] = np.interp(new_idx, old_idx, J2[:, c])
        return J2_resampled
    else:
        return J2


def infer_grid_from_ctags_len(length: int, prefer_square: bool = True) -> Tuple[int, int]:
    """
    Infer NEX, NEY from the total number of elements found in a ctags file.
    Tries to return a grid close to square; falls back to (length,1) if nothing sensible found.
    """
    sqrt_n = int(round(np.sqrt(length)))
    if sqrt_n * sqrt_n == length:
        return sqrt_n, sqrt_n

    # try common defaults, e.g., 300x300
    if length % 300 == 0 and length // 300 == 300:
        return 300, 300

    # find factor pair closest to square
    factors = []
    for a in range(1, int(np.sqrt(length)) + 1):
        if length % a == 0:
            factors.append((a, length // a))
    if factors:
        return min(factors, key=lambda t: abs(t[0] - t[1]))

    return length, 1


def load_ctags_mpi(path: str, incr: int, NVX: int):
    """
    Load and stitch ctags{incr}_r*.out files for Y-decomposed MPI.

    Returns:
        ctags2d (NEY, NEX)
    """
    pattern = os.path.join(path, f"ctags{incr}_r*.out")
    files = glob.glob(pattern)

    if not files:
        return None

    # sort by rank number
    def get_rank(f):
        m = re.search(r"_r(\d+)", f)
        return int(m.group(1)) if m else 0

    files = sorted(files, key=get_rank)

    rows = []
    for f in files:
        data = np.loadtxt(f, dtype=int)

        # ensure 2D (important if single row)
        if data.ndim == 1:
            data = data.reshape(1, -1)

        rows.append(data)

    # stack along Y direction
    full = np.hstack(rows)

    # NOTE: your file writes:
    # outer loop vx, inner vy → transposed relative to expected
    # so we must transpose
    full = full.T  # now shape (NEY, NEX)

    return full


def build_bigfield_vectorized(ctags2d: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Vectorized building of bigctags (2x upsampled) and bigfield (color indices)
    Returns: bigctags (NEY2 x NEX2), bigfield (NEY2 x NEX2)
    Color indices: background=64, cell interior=62, borders=61 (consistent with MATLAB script)
    """
    # Expand to 2x2 blocks
    bigctags = np.repeat(np.repeat(ctags2d, 2, axis=0), 2, axis=1)

    # Base color array: default white-ish (index 64)
    bigfield = np.full_like(bigctags, 64, dtype=int)

    # cell interiors -> 62
    bigfield[bigctags != 0] = 62

    # border detection: compare with up/left neighbors (avoid wrap-around)
    border = np.zeros_like(bigctags, dtype=bool)
    border[1:, :] |= (bigctags[1:, :] != bigctags[:-1, :])
    border[:, 1:] |= (bigctags[:, 1:] != bigctags[:, :-1])
    bigfield[border] = 61

    return bigctags, bigfield


def load_pstrain_mpi(path: str, incr: int):
    pattern = os.path.join(path, f"pstrain{incr}_r*.out")
    files = glob.glob(pattern)

    if not files:
        return None

    def get_rank(f):
        import re
        m = re.search(r"_r(\d+)", f)
        return int(m.group(1)) if m else 0

    files = sorted(files, key=get_rank)

    data_list = []
    for f in files:
        arr = np.loadtxt(f)
        data_list.append(arr)

    return np.vstack(data_list)


def compute_quiver_points(pstrain_arr: np.ndarray, NEX: int, NEY: int,
                          REF: float = 0.1, AMPL: float = 6.0,
                          SPACING: int = 3, SHOWTR: float = 0.03):
    """
    Compute xpos, ypos, strx, stry arrays matching MATLAB script.
    """
    magnitudes = pstrain_arr[:, 0] * 1e-6
    strsc = (magnitudes / REF) * AMPL

    xpos = []
    ypos = []
    strx = []
    stry = []
    for ey in range(1, NEY + 1, SPACING):
        for ex in range(1, NEX + 1, SPACING):
            e_idx = (ex - 1) + (ey - 1) * NEX
            if e_idx >= magnitudes.size:
                continue
            if magnitudes[e_idx] > SHOWTR:
                xpos.append(ex * 2)
                ypos.append(ey * 2)
                comp_x = pstrain_arr[e_idx, 1] / 1000.0 if pstrain_arr.shape[1] > 1 else 0.0
                comp_y = pstrain_arr[e_idx, 2] / 1000.0 if pstrain_arr.shape[1] > 2 else 0.0
                strx.append(strsc[e_idx] * comp_x)
                stry.append(strsc[e_idx] * comp_y)
    return np.array(xpos), np.array(ypos), np.array(strx), np.array(stry)


def frame_rgb_from_bigfield(bigfield: np.ndarray, cmap_array: np.ndarray) -> np.ndarray:
    """
    Convert bigfield indices (1..M) to an RGB image using the custom colormap array.
    bigfield values are treated as 1-based indices like MATLAB; we convert to 0-based.
    Returns RGB uint8 image (H, W, 3).
    """
    max_idx = cmap_array.shape[0]
    idxs = np.clip(bigfield - 1, 0, max_idx - 1).astype(int)
    # cmap_array is floats in [0,1]; map to RGB uint8
    rgb = cmap_array[idxs]  # shape (H, W, 3)
    img = (rgb * 255.0).astype(np.uint8)
    return img


def make_movie(path: str,
               start: int = 0,
               stop: int = 3000,
               step_stride: int = 5,
               out_filename: str = "simmovie.mp4",
               fps: int = 25,
               NEX: Optional[int] = None,
               NEY: Optional[int] = None,
               output_size: Tuple[int, int] = (700, 700)):
    """
    Main driver that streams frames from files in `path` to `out_filename`.
    """

    # prepare colormap
    cmap_data = jet2_colormap(64)  # (64,3)
    cmap = cmap_data  # we'll use direct indexing into this array

    # set up writer (imageio)
    writer = imageio.get_writer(out_filename, fps=fps, macro_block_size=None)

    # matplotlib figure for rendering frames
    dpi = 100
    fig_w, fig_h = output_size
    fig = plt.figure(figsize=(fig_w / dpi, fig_h / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()

    frames_written = 0

    try:
        for incr in range(start, stop, step_stride):

            # load ctags
            ctags2d = load_ctags_mpi(path, incr, NVX=NEX if NEX else 300)

            # for debugging, print binary mask
            # plt.imshow(ctags2d != 0, origin="lower")
            # plt.title(f"Binary mask incr={incr}")
            # plt.savefig(f"debug_mask_{incr}.png")
            # plt.close()

            if ctags2d is None:
                print(f"[make_movie] missing ctags for incr={incr}, skipping")
                continue
            NEY_use, NEX_use = ctags2d.shape

            if args.ney is not None and NEY_use != args.ney:
                raise ValueError(f"NEY mismatch: got {NEY_use}, expected {args.ney}")

            # Build bigfield vectorized
            bigctags, bigfield = build_bigfield_vectorized(ctags2d)

            # Load pstrain if available and build quiver
            xpos = ypos = strx = stry = np.array([])  # defaults empty
            pstrain_arr = load_pstrain_mpi(path, incr)
            if pstrain_arr is not None:
                xpos, ypos, strx, stry = compute_quiver_points(pstrain_arr, NEX_use, NEY_use)
            else:
                xpos = ypos = strx = stry = np.array([])

            # Build RGB image from indices
            rgb_img = frame_rgb_from_bigfield(bigfield, cmap)

            # Render overlay (quiver) on top of rgb_img using matplotlib canvas
            ax.clear()
            ax.set_axis_off()
            ax.imshow(rgb_img, origin="lower", interpolation="nearest")
            if xpos.size:
                ax.quiver(xpos, ypos, strx, stry,
                          angles="xy", scale_units="xy", scale=1, color="k", width=0.002)
                ax.quiver(xpos, ypos, -strx, -stry,
                          angles="xy", scale_units="xy", scale=1, color="k", width=0.002)

            # render to buffer
            fig.canvas.draw()
            w, h = fig.canvas.get_width_height()
            buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
            img = buf.reshape(h, w, 4)
            img = img[:,:,:3] # drop the alpha (opacity) channel

            # write frame
            writer.append_data(img)
            frames_written += 1

            if frames_written % 10 == 0:
                print(f"[make_movie] wrote {frames_written} frames, last incr={incr}")

    finally:
        writer.close()
        plt.close(fig)
        print(f"[make_movie] finished. Output saved to '{out_filename}' ({frames_written} frames).")


def parse_args(argv):
    p = argparse.ArgumentParser(description="Create movie from ctags/pstrain outputs.")
    p.add_argument("--path", "-p", default="../output", help="Directory with ctags#.out and pstrain#.out files.")
    p.add_argument("--start", type=int, default=0, help="Start increment (inclusive).")
    p.add_argument("--stop", type=int, default=3000, required=True, help="Stop increment (exclusive).")
    p.add_argument("--step", type=int, default=5, help="Stride between increments")
    p.add_argument("--out", "-o", default="./simmovie.mp4", help="Output movie filename.")
    p.add_argument("--fps", type=int, default=25, help="Frames per second.")
    p.add_argument("--nex", type=int, default=None, required=True, help="NEX (number of elements in x)")
    p.add_argument("--ney", type=int, default=None, help="NEY (optional).")
    return p.parse_args(argv)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    make_movie(path=args.path,
               start=args.start,
               stop=args.stop,
               step_stride=args.step,
               out_filename=args.out,
               fps=args.fps,
               NEX=args.nex,
               NEY=args.ney,
               output_size=(700, 700))