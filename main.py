#! /usr/bin/python3

from utils.gnc_lin_trim import Lin_Trim
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    lintrim = Lin_Trim()

    altitude_ft_array = np.array([5000])
    V_fps_array = np.array([250,300,350,400,450,500])
    alpha_rad_array = np.array([-2,0,2,4,6])*(np.pi/180)

    count_array = np.arange(len(altitude_ft_array)*len(alpha_rad_array)*len(V_fps_array))

    for altitude in altitude_ft_array:
        for V_fps in V_fps_array:
            for alpha_rad in alpha_rad_array:
