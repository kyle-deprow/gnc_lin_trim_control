#! /usr/bin/env/python3

import numpy as np

# Data structure to hold values specific to gnc problem
class gnc_data_structure():
    def __init__(self,count_array):
        self.count = count_array
        self.len_count = len(self.count)

        self.altitude = np.zeros(self.len_count)
        self.V_fps = np.zeros(self.len_count)
        self.alpha_rad = np.zeros(self.len_count)
        self.beta_rad = np.zeros(self.len_count)
        self.cx = np.zeros(self.len_count)
        self.cy = np.zeros(self.len_count)
        self.cz = np.zeros(self.len_count)
        self.Nz = np.zeros(self.len_count)
        self.alpha_dot = np.zeros(self.len_count)
        self.beta_dot = np.zeros(self.len_count)
        self.del_el_deg = np.zeros(self.len_count)
        self.del_ail_deg = np.zeros(self.len_count)
        self.del_rud_deg = np.zeros(self.len_count)
        self.Alon = np.zeros([self.len_count,4,4])
        self.Alat = np.zeros([self.len_count,4,4])
        self.long_half_amp = np.zeros([self.len_count,2])
        self.lat_half_amp = np.zeros([self.len_count,3])
