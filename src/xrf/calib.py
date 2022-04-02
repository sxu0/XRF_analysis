"""
calib.py

Produces a calibration curve for XRF setup based on
reference energies of radioactive sources.

Author: Shiqi Xu
"""

from pathlib import Path
from typing import Tuple, List, Callable

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chisquare


def gaussian(x: np.ndarray, height: float, centre: float, std: float):
    return height * np.exp(-((x - centre) ** 2) / (2 * std**2))


def line(x: np.ndarray, slope: float, intercept: float):
    return slope * x + intercept


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
    channels: np.ndarray,
    counts: np.ndarray,
    first_channel: float,
    last_channel: float,
    guess: List[float],
    sample: str,
    save_fig: bool = False,
    path_save: Path = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Fits a Gaussian to an energy peak, and calculates the peak centre.
    Produces a plot.

    Args:
        channels (np.ndarray[float]): Detector channel corresponding to an energy bin.
            Alternately, energy levels.
        counts (np.ndarray[int]): Counts seen in each channel.
        first_channel (float): Lower bound on channels (energy levels) included in fit.
        last_channel (float): Upper bound on channels (energy levels) included in fit.
        guess (List[float]): Guesses for Gaussian parameters, [height, centre, std].
        sample (str): Name of sample. Used in plot title.
        save_fig (bool, optional): Whether to save output plot. Defaults to False.
        path_save (Path, optional): Path to save output plot. Defaults to None.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Fit parameters and their uncertainties.
    """
    peak_fit, fit_err_cov = curve_fit(
        gaussian,
        channels[first_channel:last_channel],
        counts[first_channel:last_channel],
        p0=guess,
    )
    fit_err = np.sqrt(np.diag(fit_err_cov))

    scale = last_channel - first_channel
    peak_fit_x = np.arange(first_channel - scale / 10, last_channel - scale / 10)
    peak_fit_y = gaussian(peak_fit_x, peak_fit[0], peak_fit[1], peak_fit[2])

    peak_centre = round(peak_fit[1], 2)
    peak_centre_err = round(fit_err[1], 2)

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
        110,
        395,
        "centre:   $" + str(peak_centre) + " \pm " + str(peak_centre_err) + "$",
        ha="left",
        va="top",
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
    peak_centre_errs: np.ndarray,
    energies: np.ndarray,
    guess: List[float],
    sample: str,
    save_fig: bool = False,
    path_save: Path = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Produce a calibration curve given peak locations and test known energies
    from literature.

    Args:
        peak_centres (np.ndarray[float]): Locations of peaks (centre channel).
        energies (np.ndarray[float]): Energies of each peak, from literature.
        guess (List[float]): Guesses for [slope, intercept] of calibration line.
        sample (str): Name of sample used for calibration. Used in plot title.
        save_fig (bool, optional): Whether to save output plot. Defaults to False.
        path_save (Path, optional): Path to save output plot. Defaults to None.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Linear fit parameters and their uncertainties.
    """
    calib_fit, calib_err_cov = curve_fit(
        line,
        peak_centres,
        energies,
        p0=guess,
        sigma=peak_centre_errs,
    )
    channel_range = int(max(peak_centres) - min(peak_centres))
    fit_x = np.arange(
        int(min(peak_centres) - channel_range / 2),
        int(max(peak_centres) + channel_range / 2),
    )
    fit_y = line(fit_x, calib_fit[0], calib_fit[1])

    expected_y = line(peak_centres, calib_fit[0], calib_fit[1])
    chisq_fit, p_fit = chisquare(energies, expected_y, ddof=len(calib_fit))
    calib_err = np.diag(calib_err_cov) * max(1, np.sqrt(chisq_fit))

    plt.figure()
    plt.plot(fit_x, fit_y)
    plt.errorbar(
        peak_centres, energies, peak_centre_errs, fmt="none", ecolor="firebrick"
    )
    plt.plot(peak_centres, energies, ".")
    plt.title(sample + " Calibration Curve")
    plt.xlabel("Channel $N$")
    plt.ylabel("Energy $E$ (keV)")
    plt.text(
        105,
        400,
        "$E = ("
        + str(round(calib_fit[0], 8))
        + " \pm "
        + str(round(calib_err[0], 8))
        + ") N + ("
        + str(round(calib_fit[1], 4))
        + " \pm "
        + str(round(calib_err[1], 4))
        + ")$",
        ha="left",
        va="top",
        transform=None,
    )
    if save_fig:
        plt.savefig(path_save)
        plt.close()
    else:
        plt.show()

    return calib_fit, calib_err
