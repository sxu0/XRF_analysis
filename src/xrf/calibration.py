"""
calibration.py

Produces a calibration curve for XRF setup based on
reference energies of radioactive sources.

Author: Shiqi Xu
"""

from pathlib import Path
from typing import Tuple, List, Callable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gaussian(x: np.ndarray, height: float, centre: float, std: float):
    return height * np.exp(-((x - centre) ** 2) / (2 * std ** 2))


def line(x: np.ndarray, slope: float, intercept: float):
    return slope * x + intercept


def fit_points(
    x_data: np.ndarray, y_data: np.ndarray, func: Callable, guess: List[float]
) -> Tuple[np.ndarray, np.ndarray]:
    fitted, err_cov = curve_fit(func, x_data, y_data, p0=guess)

    return fitted, err_cov


def read_data(filename: Path) -> np.ndarray:
    """Reads FastSDD data in csv format, and returns a numpy array.

    Args:
        filename (Path): Path to data file (CSV format).

    Returns:
        np.ndarray[int]: 1D array containing counts in each channel
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
) -> Tuple[np.ndarray, np.ndarray]:
    """Fits a Gaussian to an energy peak, and calculates the peak centre.
    Produces a plot.

    Args:
        counts (np.ndarray[int]): Counts seen in each channel.
        first_channel (int): Lower bound on channels to include in fit.
        last_channel (int): Upper bound on channels to include in fit.
        guess (List[float]): Guesses for Gaussian parameters, [height, centre, std].
        sample (str): Name of sample. Used in plot title.
        save_fig (bool, optional): Whether to save output plot. Defaults to False.
        path_save (Path, optional): Path to save output plot. Defaults to None.
    
    Returns:
        (array, 2D array): Fit parameters and error covariance of fit.
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
    
    return peak_fit, fit_err


def calib_curve(
    peak_centres: np.ndarray,
    energies: np.ndarray,
    guess: List[float],
    sample: str,
    save_fig: bool = False,
    path_save: Path = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Produce a calibration curve given peak locations and test known energies from literature.

    Args:
        peak_centres (np.ndarray[float]): Locations of peaks (centre channel).
        energies (np.ndarray[float]): Energies of each peak, from literature.
        guess (List[float]): Guesses for [slope, intercept] of calibration line.
        sample (str): Name of sample used for calibration. Used in plot title.
        save_fig (bool, optional): Whether to save output plot. Defaults to False.
        path_save (Path, optional): Path to save output plot. Defaults to None.
    
    Returns:
        (array, 2D array): Linear fit parameters and error covariance of fit.
    """
    calib_fit, calib_err = fit_points(peak_centres, energies, line, guess)
    fit_x = np.arange()
    fit_y = line(fit_x, calib_fit[0], calib_fit[1])
    plt.figure()
    plt.plot(fit_x, fit_y)
    plt.plot(peak_centres, energies, ".")
    plt.title(sample + " Calibration Curve")
    plt.xlabel("Channel $N$")
    plt.ylabel("Energy $E$ (keV)")
    plt.text(
        245,
        385,
        "calibration curve: $E="
        + str(round(calib_fit[0], 4))
        + "N"
        + str(round(calib_fit[1], 3))
        + "$",
        ha="center",
        va="center",
        transform=None,
    )
    if save_fig:
        plt.savefig(path_save)
        plt.close()
    else:
        plt.show()
    
    return calib_fit, calib_err
