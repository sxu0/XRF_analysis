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
# metal = "cu"
# metal = "pb"
# metal = "ag"
# metal = "ag_HR"
# metal = "cd"
# metal = "ni"
# metal = "se"
metal = "ti_HR"


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

    if metal == "pb":
        # region: Pb spectrum calibration
        pb_counts = calib.read_data(data_path / "20220330_pb_run1.csv")
        pb_peak_centre_channels = []

        pb_peak1_fit, pb_peak1_err = calib.fit_peak(
            SDD_channels,
            pb_counts,
            726,
            745,
            [5, 734, 2],
            "Pb",
            save_fig = True,
            path_save = fig_path / "20220330_pb_peak1_fit.png",
        )
        pb_peak_centre_channels.append(pb_peak1_fit[1])

        pb_peak2_fit, pb_peak2_err = calib.fit_peak(
            SDD_channels,
            pb_counts,
            817,
            864,
            [95, 842, 6.7],
            "Pb",
            save_fig = True,
            path_save = fig_path / "20220330_pb_peak2_fit.png",
        )
        pb_peak_centre_channels.append(pb_peak2_fit[1])

        pb_peak3_fit, pb_peak3_err = calib.fit_peak(
            SDD_channels,
            pb_counts,
            990,
            1019,
            [25, 1007, 2.3],
            "Pb",
            save_fig = True,
            path_save = fig_path / "20220330_pb_peak3_fit.png",
        )
        pb_peak_centre_channels.append(pb_peak3_fit[1])

        pb_peak_centre_energies = calib.line(
            np.array(pb_peak_centre_channels),
            calibration.pb210_calib_fit[0],
            calibration.pb210_calib_fit[1],
        )

        plt.figure()
        for i in range(len(pb_peak_centre_energies)):
            plt.axvline(
                x = pb_peak_centre_energies[i],
                label = "$E=" + str(round(pb_peak_centre_energies[i], 2)) + "$ keV",
                color = 'gold',
            )
        plt.plot(default_energies, pb_counts, '.', markersize = 4)
        plt.title("Pb XRF Spectrum, Default Setting Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        plt.savefig(fig_path / "pb_spectrum.png")
        plt.close()
        # endregion: Pb spectrum calibration

    if metal == "ag":
        # region: Ag spectrum calibration
        ag_counts = calib.read_data(data_path / "20220331_ag_run1.csv")
        ag_peak_centre_channels = []

        ag_peak1_fit, ag_peak1_err = calib.fit_peak(
            SDD_channels,
            ag_counts,
            235,
            264,
            [3, 247, 10],
            "Ag",
            save_fig = True,
            path_save = fig_path / "20220331_ag_peak1_fit.png",
        )
        ag_peak_centre_channels.append(ag_peak1_fit[1])

        ag_peak2_fit, ag_peak2_err = calib.fit_peak(
            SDD_channels,
            ag_counts,
            798,
            1205,
            [11, 1051, 200],
            "Ag",
            save_fig = True,
            path_save = fig_path / "20220331_ag_peak2_fit.png",
        )
        ag_peak_centre_channels.append(ag_peak2_fit[1])

        ag_peak_centre_energies = calib.line(
            np.array(ag_peak_centre_channels),
            calibration.pb210_calib_fit[0],
            calibration.pb210_calib_fit[1],
        )

        plt.figure()
        for i in range(len(ag_peak_centre_energies)):
            plt.axvline(
                x = ag_peak_centre_energies[i],
                label = "$E=" + str(round(ag_peak_centre_energies[i], 2)) + "$ keV",
                color = 'gold',
            )
        plt.plot(default_energies, ag_counts, '.', markersize = 4)
        plt.title("Ag XRF Spectrum, Default Setting Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        plt.savefig(fig_path / "ag_spectrum.png")
        plt.close()
        # endregion: Ag spectrum calibration

    if metal == "ag_HR":
        # region: Ag high-rate spectrum calibration
        ag_HR_counts = calib.read_data(data_path / "20220331_ag_high_rate.csv")
        ag_HR_peak_centre_channels = []

        ag_HR_peak1_fit, ag_HR_peak1_err = calib.fit_peak(
            SDD_channels,
            ag_HR_counts,
            116,
            138,
            [9, 125, 9],
            "Ag",
            save_fig = True,
            path_save = fig_path / "20220331_ag_HR_peak1_fit.png",
        )
        ag_HR_peak_centre_channels.append(ag_HR_peak1_fit[1])

        ag_HR_peak2_fit, ag_HR_peak2_err = calib.fit_peak(
            SDD_channels,
            ag_HR_counts,
            397,
            603,
            [27, 526, 50],
            "Ag",
            save_fig = True,
            path_save = fig_path / "20220331_ag_HR_peak2_fit.png",
        )
        ag_HR_peak_centre_channels.append(ag_HR_peak2_fit[1])

        ag_HR_peak_centre_energies = calib.line(
            np.array(ag_HR_peak_centre_channels),
            calibration.cs137_calib_fit[0],
            calibration.cs137_calib_fit[1],
        )

        plt.figure()
        for i in range(len(ag_HR_peak_centre_energies)):
            plt.axvline(
                x = ag_HR_peak_centre_energies[i],
                label = "$E=" + str(round(ag_HR_peak_centre_energies[i], 2)) + "$ keV",
                color = 'gold',
            )
        plt.plot(high_rate_energies, ag_HR_counts, '.', markersize = 4)
        plt.title("Ag XRF Spectrum, High Rate Setting Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        plt.savefig(fig_path / "ag_HR_spectrum.png")
        plt.close()
        # endregion: Ag high-rate spectrum calibration

    if metal == "cd":
        # region: Cd spectrum calibration
        cd_counts = calib.read_data(data_path / "20220331_cd_run1.csv")
        cd_peak_centre_channels = []

        cd_peak1_fit, cd_peak1_err = calib.fit_peak(
            SDD_channels,
            cd_counts,
            245,
            278,
            [7, 260, 18],
            "Cd",
            save_fig = True,
            path_save = fig_path / "20220331_cd_peak1_fit.png",
        )
        cd_peak_centre_channels.append(cd_peak1_fit[1])

        cd_peak2_fit, cd_peak2_err = calib.fit_peak(
            SDD_channels,
            cd_counts,
            683,
            700,
            [11, 694, 10],
            "Cd",
            save_fig = True,
            path_save = fig_path / "20220331_cd_peak2_fit.png",
        )
        cd_peak_centre_channels.append(cd_peak2_fit[1])

        cd_peak3_fit, cd_peak3_err = calib.fit_peak(
            SDD_channels,
            cd_counts,
            803,
            1202,
            [11, 1045, 100],
            "Cd",
            save_fig = True,
            path_save = fig_path / "20220331_cd_peak3_fit.png",
        )
        cd_peak_centre_channels.append(cd_peak3_fit[1])

        cd_peak_centre_energies = calib.line(
            np.array(cd_peak_centre_channels),
            calibration.pb210_calib_fit[0],
            calibration.pb210_calib_fit[1],
        )

        plt.figure()
        for i in range(len(cd_peak_centre_energies)):
            plt.axvline(
                x = cd_peak_centre_energies[i],
                label = "$E=" + str(round(cd_peak_centre_energies[i], 2)) + "$ keV",
                color = 'gold',
            )
        plt.plot(default_energies, cd_counts, '.', markersize = 4)
        plt.title("Cd XRF Spectrum, Default Setting Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        plt.savefig(fig_path / "cd_spectrum.png")
        plt.close()
        # endregion: Cd spectrum calibration

    if metal == "ni":
        # region: Ni spectrum calibration
        ni_counts = calib.read_data(data_path / "20220331_ni_run1.csv")
        ni_peak_centre_channels = []

        ni_peak1_fit, ni_peak1_err = calib.fit_peak(
            SDD_channels,
            ni_counts,
            672,
            710,
            [96, 690, 5],
            "Ni",
            save_fig = True,
            path_save = fig_path / "20220331_ni_peak1_fit.png",
        )
        ni_peak_centre_channels.append(ni_peak1_fit[1])

        ni_peak2_fit, ni_peak2_err = calib.fit_peak(
            SDD_channels,
            ni_counts,
            749,
            779,
            [15, 765, 3],
            "Ni",
            save_fig = True,
            path_save = fig_path / "20220331_ni_peak2_fit.png",
        )
        ni_peak_centre_channels.append(ni_peak2_fit[1])

        ni_peak_centre_energies = calib.line(
            np.array(ni_peak_centre_channels),
            calibration.pb210_calib_fit[0],
            calibration.pb210_calib_fit[1],
        )

        plt.figure()
        for i in range(len(ni_peak_centre_energies)):
            plt.axvline(
                x = ni_peak_centre_energies[i],
                label = "$E=" + str(round(ni_peak_centre_energies[i], 2)) + "$ keV",
                color = 'gold',
            )
        plt.plot(default_energies, ni_counts, '.', markersize = 4)
        plt.title("Ni XRF Spectrum, Default Setting Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        plt.savefig(fig_path / "ni_spectrum.png")
        plt.close()
        # endregion: Ni spectrum calibration

    if metal == "se":
        # region: Se spectrum calibration
        se_counts = calib.read_data(data_path / "20220331_se_run1.csv")
        se_peak_centre_channels = []

        se_peak1_fit, se_peak1_err = calib.fit_peak(
            SDD_channels,
            se_counts,
            876,
            913,
            [170, 895, 7],
            "Se",
            save_fig = True,
            path_save = fig_path / "20220331_se_peak1_fit.png",
        )
        se_peak_centre_channels.append(se_peak1_fit[1])

        se_peak2_fit, se_peak2_err = calib.fit_peak(
            SDD_channels,
            se_counts,
            982,
            1013,
            [28, 998, 3.6],
            "Se",
            save_fig = True,
            path_save = fig_path / "20220331_se_peak2_fit.png",
        )
        se_peak_centre_channels.append(se_peak2_fit[1])

        se_peak_centre_energies = calib.line(
            np.array(se_peak_centre_channels),
            calibration.pb210_calib_fit[0],
            calibration.pb210_calib_fit[1],
        )

        plt.figure()
        for i in range(len(se_peak_centre_energies)):
            plt.axvline(
                x = se_peak_centre_energies[i],
                label = "$E=" + str(round(se_peak_centre_energies[i], 2)) + "$ keV",
                color = 'gold',
            )
        plt.plot(default_energies, se_counts, '.', markersize = 4)
        plt.title("Se XRF Spectrum, Default Setting Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        plt.savefig(fig_path / "se_spectrum.png")
        plt.close()
        # endregion: Se spectrum calibration

    if metal == "ti_HR":
        # region: Ti high-rate spectrum calibration
        ti_HR_counts = calib.read_data(data_path / "20220331_ti_high_rate.csv")
        ti_HR_peak_centre_channels = []

        ti_HR_peak1_fit, ti_HR_peak1_err = calib.fit_peak(
            SDD_channels,
            ti_HR_counts,
            175,
            189,
            [110, 181, 2.2],
            "Ti",
            save_fig = True,
            path_save = fig_path / "20220331_ti_HR_peak1_fit.png",
        )
        ti_HR_peak_centre_channels.append(ti_HR_peak1_fit[1])

        ti_HR_peak2_fit, ti_HR_peak2_err = calib.fit_peak(
            SDD_channels,
            ti_HR_counts,
            193,
            205,
            [17, 198, 1.7],
            "Ti",
            save_fig = True,
            path_save = fig_path / "20220331_ti_HR_peak2_fit.png",
        )
        ti_HR_peak_centre_channels.append(ti_HR_peak2_fit[1])

        ti_HR_peak3_fit, ti_HR_peak3_err = calib.fit_peak(
            SDD_channels,
            ti_HR_counts,
            412,
            602,
            [2, 536, 50],
            "Ti",
            save_fig = True,
            path_save = fig_path / "20220331_ti_HR_peak3_fit.png",
        )
        ti_HR_peak_centre_channels.append(ti_HR_peak3_fit[1])

        ti_HR_peak_centre_energies = calib.line(
            np.array(ti_HR_peak_centre_channels),
            calibration.cs137_calib_fit[0],
            calibration.cs137_calib_fit[1],
        )

        plt.figure()
        for i in range(len(ti_HR_peak_centre_energies)):
            plt.axvline(
                x = ti_HR_peak_centre_energies[i],
                label = "$E=" + str(round(ti_HR_peak_centre_energies[i], 2)) + "$ keV",
                color = 'gold',
            )
        plt.plot(high_rate_energies, ti_HR_counts, '.', markersize = 4)
        plt.title("Ti XRF Spectrum, High Rate Setting Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        plt.savefig(fig_path / "ti_HR_spectrum.png")
        plt.close()
        # endregion: Ti high-rate spectrum calibration
