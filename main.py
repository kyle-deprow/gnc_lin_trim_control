#! /usr/bin/python3

from utils.gnc_lin_trim import Lin_Trim
from utils.f16_model import F16_Model
from utils.gnc_data_structure import gnc_data_structure
from utils.plot import plot_gnc_object
from utils.gnc_control import LQR_Controller

import numpy as np
import time

if __name__ == "__main__":
    t1 = time.time()
    lintrim = Lin_Trim(F16_Model)
    controller = LQR_Controller()

    altitude_ft_array = np.array([1000,5000])
    V_fps_array = np.array([250,300,350,400,450,500])
    alpha_rad_array = np.array([-2,0,2,4,6])*(np.pi/180)

    count_array = np.arange(1,len(altitude_ft_array)*len(alpha_rad_array)*len(V_fps_array)+1)

    gnc_struct = gnc_data_structure(count_array)

    i = 0
    for altitude in altitude_ft_array:
        for V_fps in V_fps_array:
            for alpha_rad in alpha_rad_array:
                gnc_struct.altitude[i] = altitude
                gnc_struct.V_fps[i] = V_fps
                gnc_struct.alpha_rad[i] = alpha_rad
                i += 1

    Q = 250000
    R = 0.01
    Q_end = 300
    reg_array = np.array([2])

    for i in range(len(count_array)):
        V_fps = gnc_struct.V_fps[i]
        altitude = gnc_struct.altitude[i]
        alpha_rad = gnc_struct.alpha_rad[i]

        input_conditions = np.array([V_fps,alpha_rad,altitude])
        initial_control_guess = np.array([0,-1.931,0,0,0.1485])

        xd,x,u,Az,Ay = lintrim.trim(input_conditions,initial_control_guess)
        a,b = lintrim.linearize_nonlinear_model(x,u)

        Alon,Blon = lintrim.model.build_long_state_matrices(a,b)
        Blon = Blon[:,1].reshape((Blon.shape[0],1))
        
        Alat,Blat = lintrim.model.build_lat_state_matrices(a,b)
        lon_sens, lon_metric = lintrim.perform_mode_analysis(Alon)
        lat_sens, lat_metric = lintrim.perform_mode_analysis(Alat)

        gnc_struct.beta_rad[i] = x['AOS']
        gnc_struct.Nz[i] = -Az
        gnc_struct.alpha_dot[i] = xd['alpha_dot']
        gnc_struct.beta_dot[i] = xd['beta_dot']
        gnc_struct.roll_rate_dot[i] = xd['Roll_Rate_dot']
        gnc_struct.pitch_rate_dot[i] = xd['Pitch_Rate_dot']
        gnc_struct.yaw_rate_dot[i] = xd['Yaw_Rate_dot']
        gnc_struct.del_el_deg[i] = u['Elevator']
        gnc_struct.del_ail_deg[i] = u['Aileron']
        gnc_struct.del_rud_deg[i] = u['Rudder']
        gnc_struct.Alon[i] = Alon
        gnc_struct.Alat[i] = Alat
        gnc_struct.long_half_amp[i] = np.sort(np.unique(lon_metric[0][np.where(lon_metric[0]>0)[0]]))
        dutch_roll = np.unique(lat_metric[0][np.where(np.invert(np.isnan(lat_metric[1])))[0]])
        roll_spiral = np.sort(lat_metric[0][np.where(np.isnan(lat_metric[1]))[0]])
        gnc_struct.lat_half_amp[i] = np.append(dutch_roll, roll_spiral)

        ctl = controller.create_control_object(Alon, Blon, reg_array, Q, R, Q_end)
        uoptimal,yoptimal,timeoptimal = controller.step_response(ctl)
        gnc_struct.uoptimal[i] = uoptimal
        gnc_struct.yoptimal[i] = yoptimal[:,2]*180/np.pi
        gnc_struct.timeoptimal_sec[i] = timeoptimal

    print(time.time()-t1)
    plot_gnc_object(gnc_struct)
