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


if __name__ == "__main__":

    peak_channels = np.array([866, 1010, 1046, 1247])
    test_energies = np.array([10.555, 12.305, 12.618, 15.222])
    calib_guess = [0, 0]
    calib_fit, calib_err = fit_points(peak_channels, test_energies, line, calib_guess)

    fit_x = np.linspace(0, 2047, 2048)
    fit_y = line(fit_x, calib_fit[0], calib_fit[1])
    plt.figure()
    plt.plot(fit_x, fit_y)
    plt.plot(peak_channels, test_energies, '.')
    plt.show()
    print(calib_fit)
