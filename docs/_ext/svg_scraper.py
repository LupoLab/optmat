# SPDX-FileCopyrightText: 2026 John C. Travers and contributors
#
# SPDX-License-Identifier: MIT

"""The gallery's image scraper: every figure as SVG, with a PNG for its thumbnail.

The figures are vector graphics, so they stay sharp when zoomed. The SVG files are
reproducible (no date, a fixed hash salt), and their text is drawn as paths so that
they look the same without the fonts installed. sphinx-gallery resets matplotlib's
settings before each example, so they are applied as each figure is saved.

sphinx-gallery copies an SVG thumbnail whole but scales a PNG one, and looks for the
PNG first (``sphinx_gallery.scrapers._find_image_ext``), so each figure is also saved
as a PNG of the same name, which no page references.

``docs/conf.py`` names it by its import path, so that the configuration pickles.
"""

from pathlib import PurePosixPath

import matplotlib as mpl
import matplotlib.pyplot as plt
from sphinx_gallery.scrapers import matplotlib_scraper

SVG_SETTINGS = {"svg.fonttype": "path", "svg.hashsalt": "optmat"}


def svg_scraper(block: tuple, block_vars: dict, gallery_conf: dict) -> str:
    """Save the open figures as SVG (and PNG) and return their reST.

    Parameters
    ----------
    block : tuple
        The code block executed.
    block_vars : dict
        Its variables, with the ``image_path_iterator`` of the example.
    gallery_conf : dict
        The gallery's configuration.

    Returns
    -------
    str
        The reStructuredText of the figures.
    """
    iterator = block_vars["image_path_iterator"]
    for offset, number in enumerate(plt.get_fignums(), start=1):
        path = iterator.image_path.format(len(iterator) + offset)
        plt.figure(number).savefig(PurePosixPath(path).with_suffix(".png"))
    with mpl.rc_context(SVG_SETTINGS):
        return matplotlib_scraper(
            block, block_vars, gallery_conf, format="svg", metadata={"Date": None}
        )
