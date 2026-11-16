"""Shared publication typography helpers."""
from pathlib import Path

import matplotlib
from matplotlib.text import Text


_FONT_DIR = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
_BOLD = _FONT_DIR / "DejaVuSans-Bold.ttf"
_BOLD_ITALIC = _FONT_DIR / "DejaVuSans-BoldOblique.ttf"


def bold_all_text(figure):
    """Force every text artist to an installed true bold font face."""
    for text in figure.findobj(match=Text):
        properties = text.get_fontproperties().copy()
        italic = properties.get_style() in {"italic", "oblique"}
        properties.set_file(str(_BOLD_ITALIC if italic else _BOLD))
        properties.set_weight("bold")
        text.set_fontproperties(properties)
