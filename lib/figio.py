"""Figure export and image placeholders.

save_tiff            600 dpi, white background, RGB (no alpha channel, so nothing renders black in Word or PowerPoint), LZW TIFF.
placeholder          grey box that keeps the exact place and size of a photograph / micrograph you insert yourself.
report_placeholders  where every placeholder ends up in the exported TIFF (pixels at 600 dpi and mm), so the image can be
                     cropped to the same aspect ratio and pasted on top of it.
"""
import io
from PIL import Image

PAD_IN = 0.02                                        # pad_inches used by save_tiff (needed to locate placeholders in the TIFF)
FILL, EDGE, TEXT = '#EEF0F2', '#AEB6BF', '#6B7785'


def save_tiff(fig, path, dpi=600):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight', pad_inches=PAD_IN, facecolor='white')
    buf.seek(0)
    Image.open(buf).convert('RGB').save(path, dpi=(dpi, dpi), compression='tiff_lzw')
    print('saved', path)


def _size_mm(ax):
    ax.apply_aspect(); fw, fh = ax.figure.get_size_inches(); p = ax.get_position()
    return p.width * fw * 25.4, p.height * fh * 25.4


def placeholder(ax, label, text=None, show_size=True, fontsize=6, aspect=None):
    """Turn `ax` into a placeholder box. `label` names it in report_placeholders(); `text` is written inside the box
    (default: the label), followed by the box size in mm when show_size is True. The box fills the axes, or, with
    `aspect` (width/height), is the largest box of that shape centred in the axes (as an image drawn with imshow would be)."""
    if aspect:
        ax.set_box_aspect(1 / aspect); ax.set_anchor('C')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor(FILL)
    for s in ax.spines.values():
        s.set_visible(True); s.set_color(EDGE); s.set_linewidth(0.5); s.set_linestyle((0, (2.5, 1.5)))
    w, h = _size_mm(ax)
    main = label if text is None else text                              # text='' leaves only the size (or nothing)
    lines = [main] if main else []
    if show_size:
        lines.append(f'{w:.1f} × {h:.1f} mm')
    if lines:
        ax.text(0.5, 0.5, '\n'.join(lines), transform=ax.transAxes, ha='center', va='center', fontsize=fontsize, color=TEXT,
                linespacing=1.3)
    ax._placeholder_label = label
    return ax


def report_placeholders(fig, dpi=600, show=True):
    """Position and size of every placeholder in the TIFF written by save_tiff(fig) at `dpi`.
    Returns a list of dicts (label, left_px, top_px, width_px, height_px, width_mm, height_mm, aspect)."""
    dpi0 = fig.dpi
    fig.set_dpi(dpi)                                                    # text extents (and so the tight box) depend on the dpi
    try:
        r = fig.canvas.get_renderer()
        tb = fig.get_tightbbox(r)                                       # inches, as used by bbox_inches='tight'
        x_org, y_top = tb.x0 - PAD_IN, tb.y1 + PAD_IN
        out = []
        for ax in fig.axes:
            lab = getattr(ax, '_placeholder_label', None)
            if lab is None:
                continue
            ax.apply_aspect()
            bb = ax.get_window_extent(r).transformed(fig.dpi_scale_trans.inverted())     # inches from the lower-left corner
            out.append(dict(label=lab, left_px=round((bb.x0 - x_org) * dpi), top_px=round((y_top - bb.y1) * dpi),
                            width_px=round(bb.width * dpi), height_px=round(bb.height * dpi),
                            width_mm=round(bb.width * 25.4, 2), height_mm=round(bb.height * 25.4, 2),
                            aspect=round(bb.width / bb.height, 4)))
    finally:
        fig.set_dpi(dpi0)
    if show and out:
        print(f'Placeholders in the exported TIFF ({dpi} dpi; origin = top-left corner of the image):')
        for d in out:
            print(f"  {d['label']:<28s} left {d['left_px']:>5d} px, top {d['top_px']:>5d} px, "
                  f"{d['width_px']} x {d['height_px']} px ({d['width_mm']} x {d['height_mm']} mm), width/height {d['aspect']}")
    return out
