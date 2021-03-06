"""
calibration.py

Calibration of _Default_ and _High Rate_ detector settings.

Author: Shiqi Xu
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from xrf import calib


## mode: "default" and/or "high_rate"
# mode = ["default"]
# mode = ["high_rate"]
mode = ["default", "high_rate"]


data_path = Path.cwd() / "data"
fig_path = Path.cwd() / "outputs" / "calib_radioactive"

SDD_channels = np.arange(0, 2048)


if "default" in mode:
    # region: calibration of FastSDD Default PX5 setting
    pb210_data = calib.read_data(data_path / "20220330_pb210_run1.csv")
    pb210_peak_centres = []
    pb210_peak_centre_errs = []

    pb210_peak1_fit, pb210_peak1_err = calib.fit_peak(
        SDD_channels,
        pb210_data,
        846,
        886 + 5,
        [48, 866, 5],
        "Pb-210",
        save_fig=True,
        path_save=fig_path / "20220330_pb210_peak1_fit.png",
    )
    pb210_peak_centres.append(pb210_peak1_fit[1])
    pb210_peak_centre_errs.append(pb210_peak1_err[1])

    pb210_peak2_fit, pb210_peak2_err = calib.fit_peak(
        SDD_channels,
        pb210_data,
        1003 - 4,
        1025 + 1,
        [12, 1013, 2],
        "Pb-210",
        save_fig=True,
        path_save=fig_path / "20220330_pb210_peak2_fit.png",
    )
    pb210_peak_centres.append(pb210_peak2_fit[1])
    pb210_peak_centre_errs.append(pb210_peak2_err[1])

    pb210_peak3_fit, pb210_peak3_err = calib.fit_peak(
        SDD_channels,
        pb210_data,
        1027 - 6,
        1074 + 3,
        [23, 1045, 3],
        "Pb-210",
        save_fig=True,
        path_save=fig_path / "20220330_pb210_peak3_fit.png",
    )
    pb210_peak_centres.append(pb210_peak3_fit[1])
    pb210_peak_centre_errs.append(pb210_peak3_err[1])

    pb210_peak4_fit, pb210_peak4_err = calib.fit_peak(
        SDD_channels,
        pb210_data,
        1229 - 7,
        1263 + 15,
        [4, 1246, 3],
        "Pb-210",
        save_fig=True,
        path_save=fig_path / "20220330_pb210_peak4_fit.png",
    )
    pb210_peak_centres.append(pb210_peak4_fit[1])
    pb210_peak_centre_errs.append(pb210_peak4_err[1])

    pb210_peak_centres = np.array(pb210_peak_centres)
    pb210_peak_centre_errs = np.array(pb210_peak_centre_errs)

    pb210_calib_fit, pb210_calib_err = calib.calib_curve(
        pb210_peak_centres,
        pb210_peak_centre_errs,
        [
            10.555,
            12.305,
            12.618,
            15.222,
        ],
        [0, 0],
        "Pb-210",
        save_fig=True,
        path_save=fig_path / "calib_fastSDD_default.png",
    )

    energies_default = calib.line(SDD_channels, pb210_calib_fit[0], pb210_calib_fit[1])

    plt.figure()
    plt.plot(energies_default, pb210_data, ".", markersize=4)
    plt.title("Pb-210 Natural Decay, Default FastSDD Setting")
    plt.xlabel("Energy (keV)")
    plt.ylabel("Count")
    plt.savefig(fig_path / "pb210_spectrum.png")
    plt.close()
    # endregion: calibration of FastSDD Default PX5 setting

if "high_rate" in mode:
    # region: calibration of FastSDD High Rate PX5 setting
    cs137_data = calib.read_data(data_path / "20220331_cs137_high_rate.csv")
    cs137_peak_centres = []
    cs137_peak_centre_errs = []

    cs137_peak1_fit, cs137_peak1_err = calib.fit_peak(
        SDD_channels,
        cs137_data,
        1251 - 2,
        1300 + 4,
        [43, 1274, 5],
        "Cs-137",
        save_fig=True,
        path_save=fig_path / "20220331_cs137_high_rate_peak1_fit.png",
    )
    cs137_peak_centres.append(cs137_peak1_fit[1])
    cs137_peak_centre_errs.append(cs137_peak1_err[1])

    cs137_peak2_fit, cs137_peak2_err = calib.fit_peak(
        SDD_channels,
        cs137_data,
        1429 - 1,
        1454 + 1,
        [9, 1440, 5],
        "Cs-137",
        save_fig=True,
        path_save=fig_path / "20220331_cs137_high_rate_peak2_fit.png",
    )
    cs137_peak_centres.append(cs137_peak2_fit[1])
    cs137_peak_centre_errs.append(cs137_peak2_err[1])

    cs137_peak_centres = np.array(cs137_peak_centres)
    cs137_peak_centre_errs = np.array(cs137_peak_centre_errs)

    cs137_calib_fit, cs137_calib_err = calib.calib_curve(
        cs137_peak_centres,
        cs137_peak_centre_errs,
        [30.973, 34.985],
        [0, 0],
        "Cs-137",
        save_fig=True,
        path_save=fig_path / "calib_fastSDD_high_rate.png",
    )

    energies_high_rate = calib.line(
        SDD_channels, cs137_calib_fit[0], cs137_calib_fit[1]
    )

    plt.figure()
    plt.plot(energies_high_rate, cs137_data, ".", markersize=4)
    plt.title("Cs-137 Natural Decay, High Rate FastSDD Setting")
    plt.xlabel("Energy (keV)")
    plt.ylabel("Count")
    plt.savefig(fig_path / "cs137_spectrum.png")
    plt.close()
    # endregion: calibration of FastSDD High Rate PX5 setting
