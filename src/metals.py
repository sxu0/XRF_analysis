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

    ## begin region: calibration of FastSDD Default PX5 setting

    pb210_data = calib.read_data(data_path / "20220330_pb_run1.csv")
    

    ## end region: calibration of FastSDD Default PX5 setting
