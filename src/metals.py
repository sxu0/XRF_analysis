"""
metals.py

Central XRF analysis for metal samples.

Author: Shiqi Xu
"""

from pathlib import Path
from typing import Union, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import xrf.calibration as calib


if __name__ == "__main__":

    data_path = Path.cwd() / "data"
    fig_path = Path.cwd() / "outputs"

    # region: calibration of FastSDD Default PX5 setting
    pb210_data = calib.read_data(data_path / "20220330_pb210_run1.csv")
    pb210_peak_centres = []

    pb210_peak1_fit, pb210_peak1_err = calib.fit_peak(
        pb210_data,
        846,
        886,
        [48, 866, 5],
        "Pb-210",
        save_fig = True,
        path_save = fig_path / "20220330_pb210_peak1_fit.png",
    )
    pb210_peak_centres.append(pb210_peak1_fit[1])

    pb210_peak2_fit, pb210_peak2_err = calib.fit_peak(
        pb210_data,
        1003,
        1025,
        [12, 1013, 2],
        "Pb-210",
        save_fig = True,
        path_save = fig_path / "20220330_pb210_peak2_fit.png",
    )
    pb210_peak_centres.append(pb210_peak2_fit[1])

    pb210_peak3_fit, pb210_peak3_err = calib.fit_peak(
        pb210_data,
        1027,
        1074,
        [23, 1045, 3],
        "Pb-210",
        save_fig = True,
        path_save = fig_path / "20220330_pb210_peak3_fit.png",
    )
    pb210_peak_centres.append(pb210_peak3_fit[1])

    pb210_peak4_fit, pb210_peak4_err = calib.fit_peak(
        pb210_data,
        1229,
        1263,
        [4, 1246, 3],
        "Pb-210",
        save_fig = True,
        path_save = fig_path / "20220330_pb210_peak4_fit.png",
    )
    pb210_peak_centres.append(pb210_peak4_fit[1])

    pb210_calib_fit, pb210_calib_err = calib.calib_curve(
        pb210_peak_centres,
        [10.555, 12.305, 12.618, 15.222],
        [0, 0],
        "Pb-210",
        save_fig = True,
        path_save = fig_path / "calib_fastSDD_default.png",
    )

    channels = np.arange(0, 2048)
    energies = calib.line(channels, pb210_calib_fit[0], pb210_calib_fit[1])
    # endregion: calibration of FastSDD Default PX5 setting



