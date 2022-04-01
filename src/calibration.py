"""
calibration.py

Produces a calibration curve for XRF setup based on
reference energies of radioactive sources.

Author: Shiqi Xu
"""

from pathlib import Path
from typing import Union, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gaussian(x, height, centre, std):
    return height * np.exp(-((x - centre) ** 2) / (2 * std ** 2))


def line(x, slope, intercept):
    return slope * x + intercept


def fit_points(x_data, y_data, func, guess):
    fitted, err_cov = curve_fit(func, x_data, y_data, p0=guess)

    return fitted, err_cov


def read_data(filename):
    with open(filename, "r") as file:
        line_no = 0
        for line in file:
            line_no += 1
            if line == "<<DATA>>\n":
                header_rows = line_no
            if line == "<<END>>\n":
                footer_rows = -line_no + 1
        footer_rows += line_no
        # print("header_rows:", header_rows)
        # print("footer_rows:", footer_rows)
    return np.genfromtxt(
        filename, dtype=int, skip_header=header_rows, skip_footer=footer_rows
    )


def fit_peak(
    channels: np.ndarray,
    counts: np.ndarray,
    first_channel: int,
    last_channel: int,
    guess: List[float],
    sample: str,
    path_save: Path,
    save_fig: bool = False,
):
    peak_fit, fit_err = fit_points(
        channels[first_channel:last_channel],
        counts[first_channel:last_channel],
        gaussian,
        guess,
    )
    scale = last_channel - first_channel
    peak_fit_x = np.arange(first_channel - scale / 10, last_channel - scale / 10)
    peak_fit_y = gaussian(peak_fit_x, peak_fit[0], peak_fit[1], peak_fit[2])
    peak_centroid = round(peak_fit[1], 2)
    plt.figure()
    plt.plot(
        channels[first_channel:last_channel],
        counts[first_channel:last_channel],
        "x",
        label="data",
    )
    plt.plot(peak_fit_x, peak_fit_y, label="fit")
    plt.title(sample + " peak centred around channel " + str(peak_centroid))
    plt.xlabel("Channel")
    plt.ylabel("Count")
    plt.text(
        175,
        385,
        "centroid:   " + str(peak_centroid),
        ha="center",
        va="center",
        transform=None,
    )
    plt.legend()
    if save_fig:
        plt.savefig(path_save)
        plt.close()
    else:
        plt.show()
