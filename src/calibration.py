from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gaussian(x, height, centre, std):
    return height * np.exp(-(x - centre)**2 / (2 * std ** 2))

def line(x, slope, intercept):
    return slope * x + intercept

def fit_points(x_data, y_data, func, guess):
    fitted, err_cov = curve_fit(func, x_data, y_data, p0=guess)

    return fitted, err_cov


def read_data(filename):
    return np.genfromtxt(filename, dtype=int, skip_header=12, skip_footer=73)


if __name__ == "__main__":

    output_path = Path.cwd() / "outputs"

    pb210_data = read_data("data//20220317_take01_Pb210_1523count.csv")
    channels = np.arange(0, 2048)

    pb210_peak1_fit, pb210_peak1_err = fit_points(channels[850:882], pb210_data[850:882], gaussian, [34, 866, 20])
    pb210_peak1_fit_x = np.arange(800, 933)
    pb210_peak1_fit_y = gaussian(pb210_peak1_fit_x, pb210_peak1_fit[0], pb210_peak1_fit[1], pb210_peak1_fit[2])
    pb210_peak1_centroid = round(pb210_peak1_fit[1], 2)

    plt.figure()
    plt.plot(channels[850:882], pb210_data[850:882], 'x', label="data")
    plt.plot(pb210_peak1_fit_x, pb210_peak1_fit_y, label="fit")
    plt.title("Pb-210 peak centred around channel " + str(pb210_peak1_centroid))
    plt.xlabel("Channel")
    plt.ylabel("Count")
    plt.text(175, 385, "centroid:   " + str(pb210_peak1_centroid), ha='center', va='center', transform=None)
    plt.legend()
    plt.savefig(output_path / "20220317_pb210_peak1.png")

    # peak_channels = np.array([866, 1010, 1046, 1247])
    # test_energies = np.array([10.555, 12.305, 12.618, 15.222])
    # calib_guess = [0, 0]
    # calib_fit, calib_err = fit_points(peak_channels, test_energies, line, calib_guess)

    # fit_x = np.linspace(0, 2047, 2048)
    # fit_y = line(fit_x, calib_fit[0], calib_fit[1])
    # plt.figure()
    # plt.plot(fit_x, fit_y)
    # plt.plot(peak_channels, test_energies, '.')
    # plt.show()
    # print(calib_fit)
