#! /usr/bin/env/python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

def plot_gnc_object(gnc_struct):

    pp = PdfPages('Gnc_report.pdf')

    altitudes = np.unique(gnc_struct.altitude)
    # SMALL_SIZE = 8
    # #MEDIUM_SIZE = 10
    # BIGGER_SIZE = 12

    # plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    # plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    # plt.rc('xaxis', labelsize=5)    # fontsize of the x and y labels
    # plt.rc('xtick', labelsize=5)    # fontsize of the tick labels
    # plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    # plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    # plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    for altitude in altitudes:
        plt.figure(1)
        plt.clf()
        plt.subplot(411)
        plt.title('Trim Data for Alt: '+str(altitude)+' ft')
        plt.plot(gnc_struct.count, gnc_struct.alpha_rad, 'bx', label='AOA')
        plt.ylabel(r'$\alpha$ rad')
        plt.legend()
        plt.subplot(412)
        plt.plot(gnc_struct.count, gnc_struct.V_fps, 'rx', label='Airspeed')
        plt.ylabel('V ft/s')
        plt.legend()
        plt.subplot(413)
        plt.plot(gnc_struct.count, gnc_struct.beta_rad, 'mx', label='AOS')
        plt.ylabel(r'$\beta$ rad')
        plt.legend()
        plt.subplot(414)
        plt.plot(gnc_struct.count, gnc_struct.Nz, 'gx', label='Nz')
        plt.ylabel("acceleration g's")
        plt.xlabel("Simulation Number")
        plt.legend()

        plt.figure(2)
        plt.clf()
        plt.subplot(411)
        plt.title('Trim Data for Alt: '+str(altitude)+' ft')
        plt.plot(gnc_struct.count, gnc_struct.alpha_dot, 'bx', label=r'$\dot \alpha$')
        plt.ylabel(r'rad/s')
        plt.legend()
        plt.subplot(412)
        plt.plot(gnc_struct.count, gnc_struct.beta_dot, 'rx', label=r'$\dot \beta$')
        plt.ylabel(r'rad/s')
        plt.legend()
        plt.subplot(413)
        plt.plot(gnc_struct.count, gnc_struct.roll_rate_dot, 'gx', label=r'$\dot p$')
        plt.ylabel(r'rad/s')
        plt.legend()
        plt.subplot(414)
        plt.plot(gnc_struct.count, gnc_struct.pitch_rate_dot, 'mx', label=r'$\dot q$')
        plt.ylabel(r'rad/s')
        plt.xlabel(r'Simulation Number')
        plt.legend()

        plt.figure(3)
        plt.clf()
        plt.title('Elevator, Aileron, and Rudder Positions for '+str(altitude)+' ft')
        plt.plot(gnc_struct.count, gnc_struct.del_el_deg, 'bx', label='$\delta_{Elevator}$')
        plt.plot(gnc_struct.count, gnc_struct.del_ail_deg, 'ro', label='$\delta_{Aileron}$')
        plt.plot(gnc_struct.count, gnc_struct.del_rud_deg, 'g*', label='$\delta_{Rudder}$')
        plt.legend()
        plt.ylabel(r"$^\circ's$ of Actuation")
        plt.xlabel("Simulation Number")

        plt.figure(4)
        plt.clf()
        plt.subplot(2,1,1)
        plt.title('Longitudinal Metrics, $T_{1/2}$ Amplitude for '+str(altitude)+' ft')
        plt.plot(gnc_struct.count, gnc_struct.long_half_amp[:,0], 'bx', label='Short Period')
        plt.legend()
        plt.ylabel('Time (s)')
        plt.subplot(2,1,2)
        plt.plot(gnc_struct.count, gnc_struct.long_half_amp[:,1], 'rx', label='Phugoid')
        plt.legend()
        plt.ylabel('Time (s)')
        plt.xlabel('Simulation Number')

        plt.figure(5)
        plt.clf()
        plt.subplot(3,1,1)
        plt.title('Latitudinal Metrics, $T_{1/2}$ Amplitude')
        plt.plot(gnc_struct.count, gnc_struct.lat_half_amp[:,0], 'bx', label='Dutch Roll')
        plt.legend()
        plt.ylabel('Time(s)')
        plt.subplot(3,1,2)
        plt.plot(gnc_struct.count, gnc_struct.lat_half_amp[:,1], 'rx', label='Roll')
        plt.legend()
        plt.ylabel('Time(s)')
        plt.subplot(3,1,3)
        plt.plot(gnc_struct.count, gnc_struct.lat_half_amp[:,2], 'gx', label='Spiral')
        plt.legend()
        plt.ylabel('Time(s)')

        pp.savefig(plt.figure(1))
        pp.savefig(plt.figure(2))
        pp.savefig(plt.figure(3))
        pp.savefig(plt.figure(4))
        pp.savefig(plt.figure(5))

    pp.close()

