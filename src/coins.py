"""
coins.py

Analysis of coin composition using XRF spectra.

Author: Shiqi Xu
"""


from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from xrf import calib
from metals_self_calib import avg_calib_curve


coins = {
    "CN_old": "1600s_chinese_coin",
    "CN_new": "2000s_chinese_dime",
    "CA_old": "1800s_canadian_coin",
    "CA_new": "1964_canadian_quarter",
}
# coin = list(coins.keys())
coin = ["CN_old", "CN_new", "CA_new"]

save_plots = True


if __name__ == "__main__":

    data_path = Path.cwd() / "data"
    fig_path = Path.cwd() / "outputs" / "coins"

    SDD_channels = np.arange(0, 2048)

    plt.style.use("default")

    if "CN_new" in coin:
        # region: identifying peak energies for 2000s_chinese_dime
        CN_new_counts = calib.read_data(data_path / "20220401_2000s_chinese_dime.csv")
        CN_new_peak_centre_channels = []
        CN_new_peak_centre_channel_errs = []

        CN_new_peak1_fit, CN_new_peak1_err = calib.fit_peak(
            SDD_channels,
            CN_new_counts,
            422 - 25,
            446 + 25,
            [110, 434, 55],
            "CN_new",
            show_fig=True,
        )
        CN_new_peak_centre_channels.append(CN_new_peak1_fit[1])
        CN_new_peak_centre_channel_errs.append(CN_new_peak1_err[1])

        CN_new_peak2_fit, CN_new_peak2_err = calib.fit_peak(
            SDD_channels,
            CN_new_counts,
            463 - 15,
            488 + 7,
            [20, 476, 50],
            "CN_new",
            show_fig=True,
        )
        CN_new_peak_centre_channels.append(CN_new_peak2_fit[1])
        CN_new_peak_centre_channel_errs.append(CN_new_peak2_err[1])

        CN_new_peak3_fit, CN_new_peak3_err = calib.fit_peak(
            SDD_channels,
            CN_new_counts,
            494 - 25,
            528 + 25,
            [680, 513, 40],
            "CN_new",
            show_fig=True,
        )
        CN_new_peak_centre_channels.append(CN_new_peak3_fit[1])
        CN_new_peak_centre_channel_errs.append(CN_new_peak3_err[1])

        CN_new_peak4_fit, CN_new_peak4_err = calib.fit_peak(
            SDD_channels,
            CN_new_counts,
            548 - 15,
            580 + 25,
            [110, 565, 60],
            "CN_new",
            show_fig=True,
        )
        CN_new_peak_centre_channels.append(CN_new_peak4_fit[1])
        CN_new_peak_centre_channel_errs.append(CN_new_peak4_err[1])

        CN_new_peak_centre_channels = np.array(CN_new_peak_centre_channels)
        CN_new_peak_centre_channel_errs = np.array(CN_new_peak_centre_channel_errs)

        CN_new_peak_centre_energies = calib.line(
            np.array(CN_new_peak_centre_channels),
            avg_calib_curve[0],
            avg_calib_curve[1],
        )
        CN_new_peak_centre_energy_errs = np.sqrt(
            CN_new_peak_centre_energies**2
            * (
                (avg_calib_curve[2] / avg_calib_curve[0]) ** 2
                + (CN_new_peak_centre_channel_errs / CN_new_peak_centre_channels) ** 2
            )
            + avg_calib_curve[3] ** 2
        )

        plt.figure()
        for i in range(len(CN_new_peak_centre_energies)):
            plt.axvline(
                x=CN_new_peak_centre_energies[i],
                label="$E = "
                + str(round(CN_new_peak_centre_energies[i], 2))
                + " \pm "
                + str(round(CN_new_peak_centre_energy_errs[i], 2))
                + "$ keV",
                color="orange",
            )
        plt.plot(avg_calib_curve[4], CN_new_counts, ".", markersize=4)
        plt.title("Modern Chinese Dime, Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        if save_plots:
            plt.savefig(fig_path / "CN_new_spectrum_calib.png")
        plt.close()
        # endregion: identifying peak energies for 2000s_chinese_dime

    if "CN_old" in coin:
        # region: identifying peak energies for 1600s_chinese_coin
        CN_old_counts = calib.read_data(data_path / "20220401_1600s_chinese_coin.csv")
        CN_old_peak_centre_channels = []
        CN_old_peak_centre_channel_errs = []

        CN_old_peak1_fit, CN_old_peak1_err = calib.fit_peak(
            SDD_channels,
            CN_old_counts,
            499 - 20,
            521 + 25,
            [10, 512, 40],
            "CN_old",
            show_fig=True,
        )
        CN_old_peak_centre_channels.append(CN_old_peak1_fit[1])
        CN_old_peak_centre_channel_errs.append(CN_old_peak1_err[1])

        CN_old_peak2_fit, CN_old_peak2_err = calib.fit_peak(
            SDD_channels,
            CN_old_counts,
            625 - 25,
            660 + 20,
            [400, 643, 60],
            "CN_old",
            show_fig=True,
        )
        CN_old_peak_centre_channels.append(CN_old_peak2_fit[1])
        CN_old_peak_centre_channel_errs.append(CN_old_peak2_err[1])

        CN_old_peak3_fit, CN_old_peak3_err = calib.fit_peak(
            SDD_channels,
            CN_old_counts,
            673 - 15,
            704,
            [312, 690, 50],
            "CN_old",
            show_fig=True,
        )
        CN_old_peak_centre_channels.append(CN_old_peak3_fit[1])
        CN_old_peak_centre_channel_errs.append(CN_old_peak3_err[1])

        CN_old_peak4_fit, CN_old_peak4_err = calib.fit_peak(
            SDD_channels,
            CN_old_counts,
            705,
            725 + 25,
            [60, 713, 50],
            "CN_old",
            show_fig=True,
        )
        CN_old_peak_centre_channels.append(CN_old_peak4_fit[1])
        CN_old_peak_centre_channel_errs.append(CN_old_peak4_err[1])

        CN_old_peak5_fit, CN_old_peak5_err = calib.fit_peak(
            SDD_channels,
            CN_old_counts,
            750 - 15,
            780 + 25,
            [30, 766, 40],
            "CN_old",
            show_fig=True,
        )
        CN_old_peak_centre_channels.append(CN_old_peak5_fit[1])
        CN_old_peak_centre_channel_errs.append(CN_old_peak5_err[1])

        CN_old_peak6_fit, CN_old_peak6_err = calib.fit_peak(
            SDD_channels,
            CN_old_counts,
            829 - 25,
            858 + 25,
            [10, 843, 30],
            "CN_old",
            show_fig=True,
        )
        CN_old_peak_centre_channels.append(CN_old_peak6_fit[1])
        CN_old_peak_centre_channel_errs.append(CN_old_peak6_err[1])

        CN_old_peak_centre_channels = np.array(CN_old_peak_centre_channels)
        CN_old_peak_centre_channel_errs = np.array(CN_old_peak_centre_channel_errs)

        CN_old_peak_centre_energies = calib.line(
            np.array(CN_old_peak_centre_channels),
            avg_calib_curve[0],
            avg_calib_curve[1],
        )
        CN_old_peak_centre_energy_errs = np.sqrt(
            CN_old_peak_centre_energies**2
            * (
                (avg_calib_curve[2] / avg_calib_curve[0]) ** 2
                + (CN_old_peak_centre_channel_errs / CN_old_peak_centre_channels) ** 2
            )
            + avg_calib_curve[3] ** 2
        )

        plt.figure()
        for i in range(len(CN_old_peak_centre_energies)):
            plt.axvline(
                x=CN_old_peak_centre_energies[i],
                label="$E = "
                + str(round(CN_old_peak_centre_energies[i], 2))
                + " \pm "
                + str(round(CN_old_peak_centre_energy_errs[i], 2))
                + "$ keV",
                color="orange",
            )
        plt.plot(avg_calib_curve[4], CN_old_counts, ".", markersize=4)
        plt.title("16th-Century Chinese Dime, Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        if save_plots:
            plt.savefig(fig_path / "CN_old_spectrum_calib.png")
        plt.close()
        # endregion: identifying peak energies for 1600s_chinese_coin

    if "CA_new" in coin:
        # region: identifying peak energies for 1964_canadian_quarter
        CA_new_counts = calib.read_data(data_path / "20220401_1964_canadian_quarter.csv")
        CA_new_peak_centre_channels = []
        CA_new_peak_centre_channel_errs = []

        CA_new_peak1_fit, CA_new_peak1_err = calib.fit_peak(
            SDD_channels,
            CA_new_counts,
            103 - 20,
            145 + 20,
            [2, 120, 10],
            "CA_new",
            show_fig=True,
        )
        CA_new_peak_centre_channels.append(CA_new_peak1_fit[1])
        CA_new_peak_centre_channel_errs.append(CA_new_peak1_err[1])

        CA_new_peak2_fit, CA_new_peak2_err = calib.fit_peak(
            SDD_channels,
            CA_new_counts,
            313 - 25,
            331 + 25,
            [60, 323, 10],
            "CA_new",
            show_fig=True,
        )
        CA_new_peak_centre_channels.append(CA_new_peak2_fit[1])
        CA_new_peak_centre_channel_errs.append(CA_new_peak2_err[1])

        CA_new_peak3_fit, CA_new_peak3_err = calib.fit_peak(
            SDD_channels,
            CA_new_counts,
            349 - 18,
            365 + 25,
            [13, 357, 10],
            "CA_new",
            show_fig=True,
        )
        CA_new_peak_centre_channels.append(CA_new_peak3_fit[1])
        CA_new_peak_centre_channel_errs.append(CA_new_peak3_err[1])

        CA_new_peak4_fit, CA_new_peak4_err = calib.fit_peak(
            SDD_channels,
            CA_new_counts,
            490,
            546,
            [5, 515, 15],
            "CA_new",
            show_fig=True,
        )
        CA_new_peak_centre_channels.append(CA_new_peak4_fit[1])
        CA_new_peak_centre_channel_errs.append(CA_new_peak4_err[1])

        CA_new_peak5_fit, CA_new_peak5_err = calib.fit_peak(
            SDD_channels,
            CA_new_counts,
            541,
            600,
            [5, 515, 15],
            "CA_new",
            show_fig=True,
        )
        CA_new_peak_centre_channels.append(CA_new_peak5_fit[1])
        CA_new_peak_centre_channel_errs.append(CA_new_peak5_err[1])

        CA_new_peak_centre_channels = np.array(CA_new_peak_centre_channels)
        CA_new_peak_centre_channel_errs = np.array(CA_new_peak_centre_channel_errs)

        CA_new_peak_centre_energies = calib.line(
            np.array(CA_new_peak_centre_channels),
            avg_calib_curve[0],
            avg_calib_curve[1],
        )
        CA_new_peak_centre_energy_errs = np.sqrt(
            CA_new_peak_centre_energies**2
            * (
                (avg_calib_curve[2] / avg_calib_curve[0]) ** 2
                + (CA_new_peak_centre_channel_errs / CA_new_peak_centre_channels) ** 2
            )
            + avg_calib_curve[3] ** 2
        )

        plt.figure()
        for i in range(len(CA_new_peak_centre_energies)):
            plt.axvline(
                x=CA_new_peak_centre_energies[i],
                label="$E = "
                + str(round(CA_new_peak_centre_energies[i], 2))
                + " \pm "
                + str(round(CA_new_peak_centre_energy_errs[i], 2))
                + "$ keV",
                color="orange",
            )
        plt.plot(avg_calib_curve[4], CA_new_counts, ".", markersize=4)
        plt.title("1964 Canadian Quarter, Calibrated")
        plt.xlabel("Energy (keV)")
        plt.ylabel("Count")
        plt.legend()
        if save_plots:
            plt.savefig(fig_path / "CA_new_spectrum_calib.png")
        plt.close()
        # endregion: identifying peak energies for 1964_canadian_quarter
