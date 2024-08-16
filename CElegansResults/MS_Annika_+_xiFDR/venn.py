#!/usr/bin/env python3

# VENN
# 2023 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
from typing import List
import warnings

__version = "1.0.1"
__date = "2023-01-30"

def venn(set_1: set,
         set_2: set,
         set_3: set = None,
         labels: List[str] = ["Set 1", "Set 2", "Set 3"],
         colors: List[str] = ["#4361EE", "#4CC9F0", "#F72585"],
         alpha: float = 0.6,
         contour: bool = False,
         linewidth: float = 0.5,
         title: str = "Venn Diagram",
         filename: str = None) -> plt.figure:

    """
    DESCRIPTION
    Wrapper with some defaults that I use for creating venn diagrams with the
    matplotlib-venn package (https://github.com/konstantint/matplotlib-venn).

    PARAMETERS:
        set_1: set
            first set of the venn diagram
        set_2: set
            second set of the venn diagram
        set_3: set
            [Optional]
            if specified a three set venn diagram will be drawn, if not specified
            the two set venn diagram of set_1 and set_2 will be drawn
            DEFAULT: None
        labels: List[str]
            [Optional]
            list of labels for the sets
            DEFAULT: ["Set 1", "Set 2", "Set 3"]
        colors: List[str]
            [Optional]
            list of valid matplotlib colors to use for the venn circles
            DEFAULT: ["#4361EE", "#4CC9F0", "#F72585"]
        alpha: float
            [Optional]
            color opacity
            DEFAULT: 0.6
        contour: bool
            [Optional]
            if contour should be drawn around venn circles
            DEFAULT: False
        linewidth: float
            [Optional]
            linewidth of the contour
            DEFAULT: 0.5
        title: str
            [Optional]
            title of the figure
            DEFAULT: "Venn Diagram"
        filename: str
            [Optional]
            if given the venn diagram will be saved at the specified location in
            .png and .svg format with and without title
            DEFAULT: None

    RETURNS:
        plt.figure
            matplotlib figure object of the created venn diagram
    """

    fig = plt.figure()

    if set_3 is None:

        # checks
        if len(labels) > 2:
            warnings.warn("More than two labels supplied for two sets. Using first two...", UserWarning)

        if len(labels) < 2:
            raise IndexError("At least two labels have to be given if two sets are supplied!")

        if len(colors) > 2:
            warnings.warn("More than two colors supplied for two sets. Using first two...", UserWarning)

        if len(colors) < 2:
            raise IndexError("At least two colors have to be given if two sets are supplied!")

        # create 2 set venn diagram
        venn2(subsets = (len(set_1.difference(set_2)),
                         len(set_2.difference(set_1)),
                         len(set_1.intersection(set_2))),
              set_labels = tuple(labels),
              set_colors = tuple(colors),
              alpha = alpha)

        # draw contour
        if contour:
            venn2_circles(subsets = (len(set_1.difference(set_2)),
                                     len(set_2.difference(set_1)),
                                     len(set_1.intersection(set_2))),
                          linewidth = linewidth)

        # save file
        if filename is not None:
            plt.savefig(filename.split(".")[0] + "_notitle.png", dpi = 300, transparent = True, bbox_inches = "tight")
            plt.savefig(filename.split(".")[0] + "_notitle.svg", dpi = 300, transparent = True, bbox_inches = "tight")
            plt.title(title)
            plt.savefig(filename.split(".")[0] + ".png", dpi = 300, transparent = True, bbox_inches = "tight")
            plt.savefig(filename.split(".")[0] + ".svg", dpi = 300, transparent = True, bbox_inches = "tight")
        else:
            plt.title(title)

    else:

        # checks
        if len(labels) > 3:
            warnings.warn("More than three labels supplied for three sets. Using first three...", UserWarning)

        if len(labels) < 3:
            raise IndexError("At least three labels have to be given if three sets are supplied!")

        if len(colors) > 3:
            warnings.warn("More than three colors supplied for three sets. Using first three...", UserWarning)

        if len(colors) < 3:
            raise IndexError("At least three colors have to be given if three sets are supplied!")

        # create 3 set venn diagram
        venn3(subsets = (len(set_1.difference(set_2, set_3)),
                        len(set_2.difference(set_1, set_3)),
                        len(set_1.intersection(set_2).difference(set_3)),
                        len(set_3.difference(set_1, set_2)),
                        len(set_1.intersection(set_3).difference(set_2)),
                        len(set_2.intersection(set_3).difference(set_1)),
                        len(set_1.intersection(set_3).intersection(set_2))),
              set_labels = tuple(labels),
              set_colors = tuple(colors),
              alpha = alpha)

        # draw contour
        if contour:
            venn3_circles(subsets = (len(set_1.difference(set_2, set_3)),
                                     len(set_2.difference(set_1, set_3)),
                                     len(set_1.intersection(set_2).difference(set_3)),
                                     len(set_3.difference(set_1, set_2)),
                                     len(set_1.intersection(set_3).difference(set_2)),
                                     len(set_2.intersection(set_3).difference(set_1)),
                                     len(set_1.intersection(set_3).intersection(set_2))),
                          linewidth = linewidth)

        # save file
        if filename is not None:
            plt.savefig(filename.split(".")[0] + "_notitle.png", dpi = 300, transparent = True, bbox_inches = "tight")
            plt.savefig(filename.split(".")[0] + "_notitle.svg", dpi = 300, transparent = True, bbox_inches = "tight")
            plt.title(title)
            plt.savefig(filename.split(".")[0] + ".png", dpi = 300, transparent = True, bbox_inches = "tight")
            plt.savefig(filename.split(".")[0] + ".svg", dpi = 300, transparent = True, bbox_inches = "tight")
        else:
            plt.title(title)

    return fig
