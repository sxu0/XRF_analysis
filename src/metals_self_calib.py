"""
metals_fit.py

XRF analysis for metal samples, fitted instead of calibrated using Pb-210 and Cs-137.

Author: Shiqi Xu
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from xrf import calib
from calibration import (energies_default, energies_high_rate)


save_plots = False
metal = ["au", "cu", "pb", "ni", "se", "ti_HR"]
# metal = ["ti_HR"]


if __name__ == "__main__":

    data_path = Path.cwd() / "data"
    fig_path = Path.cwd() / "outputs" / "calib_metals_self"

    SDD_channels = np.arange(0, 2048)

    if "au" in metal:
        # region: Au spectrum calibration
        au_counts = calib.read_data(data_path / "20220330_au_run1.csv")
        au_peak_centre_channels = []
        au_peak_centre_channel_errs = []

        au_peak1_fit, au_peak1_err = calib.fit_peak(
            SDD_channels,
            au_counts,
            760,
            796,
            [118, 776, 5.5],
            "Au",
        )
        au_peak_centre_channels.append(au_peak1_fit[1])
        au_peak_centre_channel_errs.append(au_peak1_err[1])

        au_peak2_fit, au_peak2_err = calib.fit_peak(
            SDD_channels,
            au_counts,
            901 - 5,
            937 + 10,
            [34, 919, 7.4],
            "Au",
        )
        au_peak_centre_channels.append(au_peak2_fit[1])
        au_peak_centre_channel_errs.append(au_peak2_err[1])

        au_peak3_fit, au_peak3_err = calib.fit_peak(
            SDD_channels,
            au_counts,
            1011,
            1116,
            [5, 1069, 15],
            "Au",
        )
        au_peak_centre_channels.append(au_peak3_fit[1])
        au_peak_centre_channel_errs.append(au_peak3_err[1])

        au_peak_centre_channels = np.array(au_peak_centre_channels)
        au_peak_centre_channel_errs = np.array(au_peak_centre_channel_errs)

        au_calib_fit, au_calib_err = calib.calib_curve(
            au_peak_centre_channels,
            au_peak_centre_channel_errs,
            [9.705, 11.432, 13.383],
            [0, 0],
            "Au",
            save_fig=save_plots,
            path_save=fig_path / "au_calib_curve.png",
        )

        au_peak_centre_energies = calib.line(
            np.array(au_peak_centre_channels),
            au_calib_fit[0],
            au_calib_fit[1],
        )
        au_peak_centre_energy_errs = np.sqrt(
            au_peak_centre_energies**2
            * (
                (au_calib_err[0] / au_calib_fit[0]) ** 2
                + (au_peak_centre_channel_errs / au_peak_centre_channels) ** 2
            )
            + au_calib_err[1] ** 2
        )
        au_energies = calib.line(
            np.array(SDD_channels),
            au_calib_fit[0],
            au_calib_fit[1],
        )

        plt.figure()
        for i in range(len(au_peak_centre_energies)):
            plt.axvline(
                x=au_peak_centre_energies[i],
                label="$E = "
                + str(round(au_peak_centre_energies[i], 2))
                + " \pm "
                + str(round(au_peak_centre_energy_errs[i], 2))
                + "$ keV",
                color="darkorange",
            )
        plt.plot(au_energies, au_counts, ".", markersize=4)
        plt.title("Au Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        if save_plots:
            plt.savefig(fig_path / "au_spectrum_calib.png")
        plt.close()
        # endregion: Au spectrum calibration

    if "cu" in metal:
        # region: Cu spectrum calibration
        cu_counts = calib.read_data(data_path / "20220330_cu_run1.csv")
        cu_peak_centre_channels = []
        cu_peak_centre_channel_errs = []

        cu_peak1_fit, cu_peak1_err = calib.fit_peak(
            SDD_channels,
            cu_counts,
            629,
            659,
            [115, 643, 5.4],
            "Cu",
        )
        cu_peak_centre_channels.append(cu_peak1_fit[1])
        cu_peak_centre_channel_errs.append(cu_peak1_err[1])

        cu_peak2_fit, cu_peak2_err = calib.fit_peak(
            SDD_channels,
            cu_counts,
            696,
            724 + 7,
            [17, 711, 3],
            "Cu",
        )
        cu_peak_centre_channels.append(cu_peak2_fit[1])
        cu_peak_centre_channel_errs.append(cu_peak2_err[1])

        cu_peak_centre_channels = np.array(cu_peak_centre_channels)
        cu_peak_centre_channel_errs = np.array(cu_peak_centre_channel_errs)

        cu_calib_fit, cu_calib_err = calib.calib_curve(
            cu_peak_centre_channels,
            cu_peak_centre_channel_errs,
            [8.048, 8.905],
            [0, 0],
            "Cu",
            save_fig=save_plots,
            path_save=fig_path / "cu_calib_curve.png",
        )

        cu_peak_centre_energies = calib.line(
            np.array(cu_peak_centre_channels),
            cu_calib_fit[0],
            cu_calib_fit[1],
        )
        cu_peak_centre_energy_errs = np.sqrt(
            cu_peak_centre_energies**2
            * (
                (cu_calib_err[0] / cu_calib_fit[0]) ** 2
                + (cu_peak_centre_channel_errs / cu_peak_centre_channels) ** 2
            )
            + cu_calib_err[1] ** 2
        )
        cu_energies = calib.line(
            np.array(SDD_channels),
            cu_calib_fit[0],
            cu_calib_fit[1],
        )

        plt.figure()
        for i in range(len(cu_peak_centre_energies)):
            plt.axvline(
                x=cu_peak_centre_energies[i],
                label="$E = "
                + str(round(cu_peak_centre_energies[i], 2))
                + " \pm "
                + str(round(cu_peak_centre_energy_errs[i], 2))
                + "$ keV",
                color="darkorange",
            )
        plt.plot(cu_energies, cu_counts, ".", markersize=4)
        plt.title("Cu Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        if save_plots:
            plt.savefig(fig_path / "cu_spectrum_calib.png")
        plt.close()
        # endregion: Cu spectrum calibration

    if "pb" in metal:
        # region: Pb spectrum calibration
        pb_counts = calib.read_data(data_path / "20220330_pb_run1.csv")
        pb_peak_centre_channels = []
        pb_peak_centre_channel_errs = []

        pb_peak1_fit, pb_peak1_err = calib.fit_peak(
            SDD_channels,
            pb_counts,
            709,
            755,
            [5, 734, 5],
            "Pb",
        )
        pb_peak_centre_channels.append(pb_peak1_fit[1])
        pb_peak_centre_channel_errs.append(pb_peak1_err[1])

        pb_peak2_fit, pb_peak2_err = calib.fit_peak(
            SDD_channels,
            pb_counts,
            817,
            864,
            [95, 842, 6.7],
            "Pb",
        )
        pb_peak_centre_channels.append(pb_peak2_fit[1])
        pb_peak_centre_channel_errs.append(pb_peak2_err[1])

        pb_peak3_fit, pb_peak3_err = calib.fit_peak(
            SDD_channels,
            pb_counts,
            990,
            1054,
            [25, 1007, 8],
            "Pb",
        )
        pb_peak_centre_channels.append(pb_peak3_fit[1])
        pb_peak_centre_channel_errs.append(pb_peak3_err[1])

        pb_peak_centre_channels = np.array(pb_peak_centre_channels)
        pb_peak_centre_channel_errs = np.array(pb_peak_centre_channel_errs)

        pb_calib_fit, pb_calib_err = calib.calib_curve(
            pb_peak_centre_channels,
            pb_peak_centre_channel_errs,
            [9.185, 10.555, 12.618],
            [0, 0],
            "Pb",
            save_fig=save_plots,
            path_save=fig_path / "pb_calib_curve.png",
        )

        pb_peak_centre_energies = calib.line(
            np.array(pb_peak_centre_channels),
            pb_calib_fit[0],
            pb_calib_fit[1],
        )
        pb_peak_centre_energy_errs = np.sqrt(
            pb_peak_centre_energies**2
            * (
                (pb_calib_err[0] / pb_calib_fit[0]) ** 2
                + (pb_peak_centre_channel_errs / pb_peak_centre_channels) ** 2
            )
            + pb_calib_err[1] ** 2
        )
        pb_energies = calib.line(
            np.array(SDD_channels),
            pb_calib_fit[0],
            pb_calib_fit[1],
        )

        plt.figure()
        for i in range(len(pb_peak_centre_energies)):
            plt.axvline(
                x=pb_peak_centre_energies[i],
                label="$E = "
                + str(round(pb_peak_centre_energies[i], 2))
                + " \pm "
                + str(round(pb_peak_centre_energy_errs[i], 2))
                + "$ keV",
                color="darkorange",
            )
        plt.plot(pb_energies, pb_counts, ".", markersize=4)
        plt.title("Pb Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        if save_plots:
            plt.savefig(fig_path / "pb_spectrum_calib.png")
        plt.close()
        # endregion: Pb spectrum calibration

    if "ni" in metal:
        # region: Ni spectrum calibration
        ni_counts = calib.read_data(data_path / "20220331_ni_run1.csv")
        ni_peak_centre_channels = []
        ni_peak_centre_channel_errs = []

        ni_peak1_fit, ni_peak1_err = calib.fit_peak(
            SDD_channels,
            ni_counts,
            672,
            710,
            [96, 690, 5],
            "Ni",
        )
        ni_peak_centre_channels.append(ni_peak1_fit[1])
        ni_peak_centre_channel_errs.append(ni_peak1_err[1])

        ni_peak2_fit, ni_peak2_err = calib.fit_peak(
            SDD_channels,
            ni_counts,
            749,
            779,
            [15, 765, 3],
            "Ni",
        )
        ni_peak_centre_channels.append(ni_peak2_fit[1])
        ni_peak_centre_channel_errs.append(ni_peak2_err[1])

        ni_peak_centre_channels = np.array(ni_peak_centre_channels)
        ni_peak_centre_channel_errs = np.array(ni_peak_centre_channel_errs)

        ni_calib_fit, ni_calib_err = calib.calib_curve(
            ni_peak_centre_channels,
            ni_peak_centre_channel_errs,
            [7.478, 8.265],
            [0, 0],
            "Ni",
            save_fig=save_plots,
            path_save=fig_path / "ni_calib_curve.png",
        )

        ni_peak_centre_energies = calib.line(
            np.array(ni_peak_centre_channels),
            ni_calib_fit[0],
            ni_calib_fit[1],
        )
        ni_peak_centre_energy_errs = np.sqrt(
            ni_peak_centre_energies**2
            * (
                (ni_calib_err[0] / ni_calib_fit[0]) ** 2
                + (ni_peak_centre_channel_errs / ni_peak_centre_channels) ** 2
            )
            + ni_calib_err[1] ** 2
        )
        ni_energies = calib.line(
            np.array(SDD_channels),
            ni_calib_fit[0],
            ni_calib_fit[1],
        )

        plt.figure()
        for i in range(len(ni_peak_centre_energies)):
            plt.axvline(
                x=ni_peak_centre_energies[i],
                label="$E = "
                + str(round(ni_peak_centre_energies[i], 2))
                + " \pm "
                + str(round(ni_peak_centre_energy_errs[i], 2))
                + "$ keV",
                color="darkorange",
            )
        plt.plot(ni_energies, ni_counts, ".", markersize=4)
        plt.title("Ni Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        if save_plots:
            plt.savefig(fig_path / "ni_spectrum_calib.png")
        plt.close()
        # endregion: Ni spectrum calibration

    if "se" in metal:
        # region: Se spectrum calibration
        se_counts = calib.read_data(data_path / "20220331_se_run1.csv")
        se_peak_centre_channels = []
        se_peak_centre_channel_errs = []

        se_peak1_fit, se_peak1_err = calib.fit_peak(
            SDD_channels,
            se_counts,
            876,
            913,
            [170, 895, 7],
            "Se",
        )
        se_peak_centre_channels.append(se_peak1_fit[1])
        se_peak_centre_channel_errs.append(se_peak1_err[1])

        se_peak2_fit, se_peak2_err = calib.fit_peak(
            SDD_channels,
            se_counts,
            982,
            1013,
            [28, 998, 3.6],
            "Se",
        )
        se_peak_centre_channels.append(se_peak2_fit[1])
        se_peak_centre_channel_errs.append(se_peak2_err[1])

        se_peak_centre_channels = np.array(se_peak_centre_channels)
        se_peak_centre_channel_errs = np.array(se_peak_centre_channel_errs)

        se_calib_fit, se_calib_err = calib.calib_curve(
            se_peak_centre_channels,
            se_peak_centre_channel_errs,
            [11.222, 12.496],
            [0, 0],
            "Se",
            save_fig=save_plots,
            path_save=fig_path / "se_calib_curve.png",
        )

        se_peak_centre_energies = calib.line(
            np.array(se_peak_centre_channels),
            se_calib_fit[0],
            se_calib_fit[1],
        )
        se_peak_centre_energy_errs = np.sqrt(
            se_peak_centre_energies**2
            * (
                (se_calib_err[0] / se_calib_fit[0]) ** 2
                + (se_peak_centre_channel_errs / se_peak_centre_channels) ** 2
            )
            + se_calib_err[1] ** 2
        )
        se_energies = calib.line(
            np.array(SDD_channels),
            se_calib_fit[0],
            se_calib_fit[1],
        )

        plt.figure()
        for i in range(len(se_peak_centre_energies)):
            plt.axvline(
                x=se_peak_centre_energies[i],
                label="$E = "
                + str(round(se_peak_centre_energies[i], 2))
                + " \pm "
                + str(round(se_peak_centre_energy_errs[i], 2))
                + "$ keV",
                color="darkorange",
            )
        plt.plot(se_energies, se_counts, ".", markersize=4)
        plt.title("Se Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        if save_plots:
            plt.savefig(fig_path / "se_spectrum_calib.png")
        plt.close()
        # endregion: Se spectrum calibration

    if "ti_HR" in metal:
        # region: Ti high-rate spectrum calibration
        ti_HR_counts = calib.read_data(data_path / "20220331_ti_high_rate.csv")
        ti_HR_peak_centre_channels = []
        ti_HR_peak_centre_channel_errs = []

        ti_HR_peak1_fit, ti_HR_peak1_err = calib.fit_peak(
            SDD_channels,
            ti_HR_counts,
            175 - 3,
            189 + 3,
            [110, 181, 2.2],
            "Ti",
        )
        ti_HR_peak_centre_channels.append(ti_HR_peak1_fit[1])
        ti_HR_peak_centre_channel_errs.append(ti_HR_peak1_err[1])

        ti_HR_peak2_fit, ti_HR_peak2_err = calib.fit_peak(
            SDD_channels,
            ti_HR_counts,
            193 - 4,
            205 + 4,
            [17, 198, 1.7],
            "Ti",
        )
        ti_HR_peak_centre_channels.append(ti_HR_peak2_fit[1])
        ti_HR_peak_centre_channel_errs.append(ti_HR_peak2_err[1])

        # ti_HR_peak3_fit, ti_HR_peak3_err = calib.fit_peak(
        #     SDD_channels,
        #     ti_HR_counts,
        #     412-15,
        #     602,
        #     [2, 536, 50],
        #     "Ti",
        # )
        # ti_HR_peak_centre_channels.append(ti_HR_peak3_fit[1])
        # ti_HR_peak_centre_channel_errs.append(ti_HR_peak3_err[1])

        ti_HR_peak_centre_channels = np.array(ti_HR_peak_centre_channels)
        ti_HR_peak_centre_channel_errs = np.array(ti_HR_peak_centre_channel_errs)

        ti_HR_calib_fit, ti_HR_calib_err = calib.calib_curve(
            ti_HR_peak_centre_channels,
            ti_HR_peak_centre_channel_errs,
            [4.511, 4.932],
            [0, 0],
            "Ti",
            save_fig=save_plots,
            path_save=fig_path / "ti_HR_calib_curve.png",
        )

        ti_HR_peak_centre_energies = calib.line(
            np.array(ti_HR_peak_centre_channels),
            ti_HR_calib_fit[0],
            ti_HR_calib_fit[1],
        )
        ti_HR_peak_centre_energy_errs = np.sqrt(
            ti_HR_peak_centre_energies**2
            * (
                (ti_HR_calib_err[0] / ti_HR_calib_fit[0]) ** 2
                + (ti_HR_peak_centre_channel_errs / ti_HR_peak_centre_channels) ** 2
            )
            + ti_HR_calib_err[1] ** 2
        )
        ti_HR_energies = calib.line(
            np.array(SDD_channels),
            ti_HR_calib_fit[0],
            ti_HR_calib_fit[1],
        )

        plt.figure()
        for i in range(len(ti_HR_peak_centre_energies)):
            plt.axvline(
                x=ti_HR_peak_centre_energies[i],
                label="$E = "
                + str(round(ti_HR_peak_centre_energies[i], 2))
                + " \pm "
                + str(round(ti_HR_peak_centre_energy_errs[i], 2))
                + "$ keV",
                color="darkorange",
            )
        plt.plot(ti_HR_energies, ti_HR_counts, ".", markersize=4)
        plt.title("Ti Calibrated (High Rate Setting)")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        if save_plots:
            plt.savefig(fig_path / "ti_HR_spectrum_calib.png")
        plt.close()
        # endregion: Ti high-rate spectrum calibration

    if (
        ("au" in metal) and
        ("cu" in metal) and
        ("pb" in metal) and
        ("ni" in metal) and
        ("se" in metal) and
        ("ti_HR" in metal)
    ):
        # region: overlay calibration curves
        plt.figure()
        calib_curves = {}
        for i in range(len(metal)):
            if "_HR" in metal[i]:
                k = 0.5
            else:
                k = 1
            calib_curves[metal[i]] = {}
            calib_curves[metal[i]]["fit"] = locals()[metal[i] + "_calib_fit"]
            calib_curves[metal[i]]["fit"][0] *= k  # scale high rate energy by 1/2
            calib_curves[metal[i]]["err"] = locals()[metal[i] + "_calib_err"]
            energies = calib.line(
                SDD_channels,
                calib_curves[metal[i]]["fit"][0],
                calib_curves[metal[i]]["fit"][1],
            )
            plt.plot(
                SDD_channels, energies, linewidth=0.8, label=metal[i][:2].capitalize(),
            )
        plt.plot(
            SDD_channels, energies_default, '--', linewidth=0.8, label="Pb-210",
        )
        plt.plot(
            SDD_channels, 0.5 * energies_high_rate, '--', linewidth=0.8, label="Cs-137",
        )
        plt.style.use("seaborn")
        plt.title("Metal Calibration Curves, by Element/Isotope")
        plt.xlabel("Channel $N$")
        plt.ylabel("Energy $E$ (keV)")
        plt.legend()
        if save_plots:
            plt.savefig(fig_path / "metal_calib_curves.png")
        plt.close()
        # endregion: overlay calibration curves
