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


def gaussian(x: np.ndarray, height: float, centre: float, std: float):
    return height * np.exp(-((x - centre) ** 2) / (2 * std ** 2))


def line(x: np.ndarray, slope: float, intercept: float):
    return slope * x + intercept


def fit_points(
    x_data: np.ndarray, y_data: np.ndarray, func: function, guess: List[float]
):
    fitted, err_cov = curve_fit(func, x_data, y_data, p0=guess)

    return fitted, err_cov


def read_data(filename: Path):
    """Reads FastSDD data in csv format, and returns a numpy array.

    Args:
        filename (Path): Path to data file (CSV format).

    Returns:
        np.ndarray: 1D array containing counts in each channel
            (index corresponds to channel number).
    """
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
    counts: np.ndarray,
    first_channel: int,
    last_channel: int,
    guess: List[float],
    sample: str,
    save_fig: bool = False,
    path_save: Path = None,
):
    """Fits a Gaussian to an energy peak, and calculates the peak centre.
    Produces a plot.

    Args:
        counts (np.ndarray): Counts seen in each channel.
        first_channel (int): Lower bound on channels to include in fit.
        last_channel (int): Upper bound on channels to include in fit.
        guess (List[float]): Guesses for Gaussian parameters, [height, centre, std].
        sample (str): Name of sample. Used in plot title.
        save_fig (bool, optional): Whether to save output plot. Defaults to False.
        path_save (Path, optional): Path to save output plot. Defaults to None.
    """
    channels = np.arange(0, len(counts))
    peak_fit, fit_err = fit_points(
        channels[first_channel:last_channel],
        counts[first_channel:last_channel],
        gaussian,
        guess,
    )
    scale = last_channel - first_channel
    peak_fit_x = np.arange(first_channel - scale / 10, last_channel - scale / 10)
    peak_fit_y = gaussian(peak_fit_x, peak_fit[0], peak_fit[1], peak_fit[2])
    peak_centre = round(peak_fit[1], 2)
    plt.figure()
    plt.plot(
        channels[first_channel:last_channel],
        counts[first_channel:last_channel],
        "x",
        label="data",
    )
    plt.plot(peak_fit_x, peak_fit_y, label="fit")
    plt.title(sample + " peak centred around channel " + str(peak_centre))
    plt.xlabel("Channel")
    plt.ylabel("Count")
    plt.text(
        175,
        385,
        "centre:   " + str(peak_centre),
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
