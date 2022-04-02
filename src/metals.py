"""
metals.py

Central XRF analysis for metal samples.

Author: Shiqi Xu
"""

from pathlib import Path
from re import M

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from xrf import calib
import calibration


# 20220330_au_run1.csv
# 20220330_cu_run1.csv
# 20220330_pb_run1.csv
# 20220331_cd_run1.csv
# 20220331_ni_run1.csv
# 20220331_se_run1.csv
# 20220331_ti_run1.csv
# 20220331_ag_run1.csv
# 20220331_ag_high_rate.csv

# metal = "au"
metal = "cu"


if __name__ == "__main__":

    data_path = Path.cwd() / "data"
    fig_path = Path.cwd() / "outputs" / "calib_metals"

    SDD_channels = np.arange(0, 2048)
    default_energies = calibration.energies_default
    high_rate_energies = calibration.energies_high_rate


    if metal == "au":
        # region: Au spectrum calibration
        au_counts = calib.read_data(data_path / "20220330_au_run1.csv")
        au_peak_centre_channels = []

        au_peak1_fit, au_peak1_err = calib.fit_peak(
            SDD_channels,
            au_counts,
            760,
            796,
            [118, 776, 5.5],
            "Au",
            save_fig = True,
            path_save = fig_path / "20220330_au_peak1_fit.png",
        )
        au_peak_centre_channels.append(au_peak1_fit[1])

        au_peak2_fit, au_peak2_err = calib.fit_peak(
            SDD_channels,
            au_counts,
            901,
            937,
            [34, 919, 7.4],
            "Au",
            save_fig = True,
            path_save = fig_path / "20220330_au_peak2_fit.png",
        )
        au_peak_centre_channels.append(au_peak2_fit[1])

        au_peak3_fit, au_peak3_err = calib.fit_peak(
            SDD_channels,
            au_counts,
            1052,
            1092,
            [5, 1069, 0.7],
            "Au",
            save_fig = True,
            path_save = fig_path / "20220330_au_peak3_fit.png",
        )
        au_peak_centre_channels.append(au_peak3_fit[1])

        au_peak_centre_energies = calib.line(
            np.array(au_peak_centre_channels),
            calibration.pb210_calib_fit[0],
            calibration.pb210_calib_fit[1],
        )

        plt.figure()
        for i in range(len(au_peak_centre_energies)):
            plt.axvline(
                x = au_peak_centre_energies[i],
                label = "$E=" + str(round(au_peak_centre_energies[i], 2)) + "$ keV",
                color = 'gold',
            )
        plt.plot(default_energies, au_counts, '.', markersize = 4)
        plt.title("Au XRF Spectrum, Default Setting Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        plt.savefig(fig_path / "au_spectrum.png")
        plt.close()
        # endregion: Au spectrum calibration

    if metal == "cu":
        # region: Cu spectrum calibration
        cu_counts = calib.read_data(data_path / "20220330_cu_run1.csv")
        cu_peak_centre_channels = []

        cu_peak1_fit, cu_peak1_err = calib.fit_peak(
            SDD_channels,
            cu_counts,
            629,
            659,
            [115, 643, 5.4],
            "Cu",
            save_fig = True,
            path_save = fig_path / "20220330_cu_peak1_fit.png",
        )
        cu_peak_centre_channels.append(cu_peak1_fit[1])

        cu_peak2_fit, cu_peak2_err = calib.fit_peak(
            SDD_channels,
            cu_counts,
            696,
            724,
            [17, 711, 3],
            "Cu",
            save_fig = True,
            path_save = fig_path / "20220330_cu_peak2_fit.png",
        )
        cu_peak_centre_channels.append(cu_peak2_fit[1])

        cu_peak_centre_energies = calib.line(
            np.array(cu_peak_centre_channels),
            calibration.pb210_calib_fit[0],
            calibration.pb210_calib_fit[1],
        )

        plt.figure()
        for i in range(len(cu_peak_centre_energies)):
            plt.axvline(
                x = cu_peak_centre_energies[i],
                label = "$E=" + str(round(cu_peak_centre_energies[i], 2)) + "$ keV",
                color = 'gold',
            )
        plt.plot(default_energies, cu_counts, '.', markersize = 4)
        plt.title("Cu XRF Spectrum, Default Setting Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        plt.savefig(fig_path / "cu_spectrum.png")
        plt.close()
        # endregion: Cu spectrum calibration
