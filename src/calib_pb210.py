from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gaussian(x, height, centre, std):
    return height * np.exp(-((x - centre) ** 2) / (2 * std ** 2))


def line(x, slope, intercept):
    return slope * x + intercept


def fit_points(x_data, y_data, func, guess):
    fitted, err_cov = curve_fit(func, x_data, y_data, p0=guess)

    return fitted, err_cov


def read_data(filename):
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


if __name__ == "__main__":

    output_path = Path.cwd() / "outputs"

    pb210_data = read_data(Path.cwd() / "data" / "20220317_take01_Pb210_1523count.csv")
    # print(len(pb210_data))
    channels = np.arange(0, 2048)

    ## Pb-210 peak 1
    pb210_peak1_fit, pb210_peak1_err = fit_points(
        channels[850:882], pb210_data[850:882], gaussian, [34, 866, 20]
    )
    pb210_peak1_fit_x = np.arange(800, 932)
    pb210_peak1_fit_y = gaussian(
        pb210_peak1_fit_x, pb210_peak1_fit[0], pb210_peak1_fit[1], pb210_peak1_fit[2]
    )
    pb210_peak1_centroid = round(pb210_peak1_fit[1], 2)
    plt.figure()
    plt.plot(channels[850:882], pb210_data[850:882], "x", label="data")
    plt.plot(pb210_peak1_fit_x, pb210_peak1_fit_y, label="fit")
    plt.title("Pb-210 peak centred around channel " + str(pb210_peak1_centroid))
    plt.xlabel("Channel")
    plt.ylabel("Count")
    plt.text(
        175,
        385,
        "centroid:   " + str(pb210_peak1_centroid),
        ha="center",
        va="center",
        transform=None,
    )
    plt.legend()
    # plt.savefig(output_path / "20220317_pb210_peak1.png")
    plt.close()

    ## Pb-210 peak 2
    pb210_peak2_fit, pb210_peak2_err = fit_points(
        channels[998:1023], pb210_data[998:1023], gaussian, [10, 1010, 11]
    )
    pb210_peak2_fit_x = np.arange(998 - 50, 1023 + 50)
    pb210_peak2_fit_y = gaussian(
        pb210_peak2_fit_x, pb210_peak2_fit[0], pb210_peak2_fit[1], pb210_peak2_fit[2]
    )
    pb210_peak2_centroid = round(pb210_peak2_fit[1], 2)
    plt.figure()
    plt.plot(channels[998:1023], pb210_data[998:1023], "x", label="data")
    plt.plot(pb210_peak2_fit_x, pb210_peak2_fit_y, label="fit")
    plt.title("Pb-210 peak centred around channel " + str(pb210_peak2_centroid))
    plt.xlabel("Channel")
    plt.ylabel("Count")
    plt.text(
        175,
        385,
        "centroid:   " + str(pb210_peak2_centroid),
        ha="center",
        va="center",
        transform=None,
    )
    plt.legend()
    # plt.savefig(output_path / "20220317_pb210_peak2.png")
    plt.close()

    ## Pb-210 peak 3
    pb210_peak3_fit, pb210_peak3_err = fit_points(
        channels[1024:1068], pb210_data[1024:1068], gaussian, [18, 1046, 22]
    )
    pb210_peak3_fit_x = np.arange(1024 - 50, 1068 + 50)
    pb210_peak3_fit_y = gaussian(
        pb210_peak3_fit_x, pb210_peak3_fit[0], pb210_peak3_fit[1], pb210_peak3_fit[2]
    )
    pb210_peak3_centroid = round(pb210_peak3_fit[1], 2)
    plt.figure()
    plt.plot(channels[1024:1068], pb210_data[1024:1068], "x", label="data")
    plt.plot(pb210_peak3_fit_x, pb210_peak3_fit_y, label="fit")
    plt.title("Pb-210 peak centred around channel " + str(pb210_peak3_centroid))
    plt.xlabel("Channel")
    plt.ylabel("Count")
    plt.text(
        175,
        385,
        "centroid:   " + str(pb210_peak3_centroid),
        ha="center",
        va="center",
        transform=None,
    )
    plt.legend()
    # plt.savefig(output_path / "20220317_pb210_peak3.png")
    plt.close()

    ## Pb-210 peak 4
    pb210_peak4_fit, pb210_peak4_err = fit_points(
        channels[1236:1260], pb210_data[1236:1260], gaussian, [6, 1248, 6]
    )
    pb210_peak4_fit_x = np.arange(1236 - 25, 1260 + 25)
    pb210_peak4_fit_y = gaussian(
        pb210_peak4_fit_x, pb210_peak4_fit[0], pb210_peak4_fit[1], pb210_peak4_fit[2]
    )
    pb210_peak4_centroid = round(pb210_peak4_fit[1], 2)
    plt.figure()
    plt.plot(channels[1236:1260], pb210_data[1236:1260], "x", label="data")
    plt.plot(pb210_peak4_fit_x, pb210_peak4_fit_y, label="fit")
    plt.title("Pb-210 peak centred around channel " + str(pb210_peak4_centroid))
    plt.xlabel("Channel")
    plt.ylabel("Count")
    plt.text(
        175,
        385,
        "centroid:   " + str(pb210_peak4_centroid),
        ha="center",
        va="center",
        transform=None,
    )
    plt.legend()
    # plt.savefig(output_path / "20220317_pb210_peak4.png")
    plt.close()

    ## Pb-210 calibration curve
    peak_channels = np.array(
        [
            pb210_peak1_centroid,
            pb210_peak2_centroid,
            pb210_peak3_centroid,
            pb210_peak4_centroid,
        ]
    )
    test_energies = np.array([10.555, 12.305, 12.618, 15.222])  # from literature table
    calib_guess = [0, 0]
    calib_fit, calib_err = fit_points(peak_channels, test_energies, line, calib_guess)

    fit_x = np.arange(700, 1400)
    fit_y = line(fit_x, calib_fit[0], calib_fit[1])
    plt.figure()
    plt.plot(fit_x, fit_y)
    plt.plot(peak_channels, test_energies, ".")
    plt.title("Pb-210 Calibration Curve")
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
    # plt.savefig(output_path / "20220317_pb210_calib_curve.png")
    plt.close()

    ## test Pb-210 calibration curve on Cs-137
    cs137_data = read_data(Path.cwd() / "data" / "20220323_take03_Cs137.csv")
    # print(len(cs137_data))

    cs137_peak1_fit, cs137_peak1_err = fit_points(
        channels[11:19], cs137_data[11:19], gaussian, [425, 14, 2]
    )
    cs137_peak1_fit_x = np.arange(11 - 5, 19 + 5)
    cs137_peak1_fit_y = gaussian(
        cs137_peak1_fit_x, cs137_peak1_fit[0], cs137_peak1_fit[1], cs137_peak1_fit[2]
    )
    cs137_peak1_centroid = round(cs137_peak1_fit[1], 2)
    plt.figure()
    plt.plot(channels[11:19], cs137_data[11:19], "x", label="data")
    plt.plot(cs137_peak1_fit_x, cs137_peak1_fit_y, label="fit")
    plt.title("Cs-137 peak centred around channel " + str(cs137_peak1_centroid))
    plt.xlabel("Channel")
    plt.ylabel("Count")
    plt.text(
        175,
        385,
        "centroid:   " + str(cs137_peak1_centroid),
        ha="center",
        va="center",
        transform=None,
    )
    plt.legend()
    # plt.savefig(output_path / "20220323_cs137_peak1.png")
    plt.close()

    ## energy of Cs-137 peak calculated from Pb-210 calibration curve
    print(
        "Cs-137 channel "
        + str(cs137_peak1_centroid)
        + " energy:\t"
        + str(round(line(cs137_peak1_centroid, calib_fit[0], calib_fit[1]), 4))
        + " keV"
    )
