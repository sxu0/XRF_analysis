"""
metals.py

Central XRF analysis for metal samples.

Author: Shiqi Xu
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from xrf import calib


## mode: "default" or "high_rate"
# mode = "default"
mode = "high_rate"


if __name__ == "__main__":

    data_path = Path.cwd() / "data"
    fig_path = Path.cwd() / "outputs"

    if mode == "default":
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
        energies_default = calib.line(channels, pb210_calib_fit[0], pb210_calib_fit[1])
        # endregion: calibration of FastSDD Default PX5 setting

    elif mode == "high_rate":
        # region: calibration of FastSDD High Rate PX5 setting
        cs137_data = calib.read_data(data_path / "20220331_cs137_high_rate.csv")
        cs137_peak_centres = []

        cs137_peak1_fit, cs137_peak1_err = calib.fit_peak(
            cs137_data,
            1251,
            1300,
            [43, 1274, 5],
            "Cs-137",
            save_fig = True,
            path_save = fig_path / "20220331_cs137_high_rate_peak1_fit.png",
        )
        cs137_peak_centres.append(cs137_peak1_fit[1])

        cs137_peak2_fit, cs137_peak2_err = calib.fit_peak(
            cs137_data,
            1429,
            1454,
            [9, 1440, 5],
            "Cs-137",
            save_fig = True,
            path_save = fig_path / "20220331_cs137_high_rate_peak2_fit.png",
        )
        cs137_peak_centres.append(cs137_peak2_fit[1])

        cs137_calib_fit, cs137_calib_err = calib.calib_curve(
            cs137_peak_centres,
            [30.973, 34.985],
            [0, 0],
            "Cs-137",
            save_fig = True,
            path_save = fig_path / "calib_fastSDD_high_rate.png",
        )

        channels = np.arange(0, 2048)
        energies = calib.line(channels, cs137_calib_fit[0], cs137_calib_fit[1])
        # endregion: calibration of FastSDD High Rate PX5 setting
