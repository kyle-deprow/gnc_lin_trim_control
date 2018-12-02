#! /usr/bin/env/python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import sys

def plot_gnc_object(gnc_struct):

    pp = PdfPages('Gnc_report.pdf')

    altitudes = np.unique(gnc_struct.altitude)
    for altitude in altitudes:
        plt.figure(1)
        plt.clf()
        plt.subplot(411)
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
        plt.subplot(211)
        plt.plot(gnc_struct.count, gnc_struct.alpha_dot, 'bx', label=r'$\dot \alpha$')
        plt.ylabel(r'$\dot \alpha$ rad/s')
        plt.legend()
        plt.subplot(212)
        plt.plot(gnc_struct.count, gnc_struct.beta_dot, 'rx', label=r'$\dot \beta$')
        plt.xlabel(r'$\dot \beta$ rad/s')
        plt.legend()

        plt.figure(3)
        plt.clf()
        plt.title('Elevator, Aileron, and Rudder Positions')
        plt.plot(gnc_struct.count, gnc_struct.del_el_deg, 'bx', label='$\delta_{Elevator}$')
        plt.plot(gnc_struct.count, gnc_struct.del_ail_deg, 'ro', label='$\delta_{Aileron}$')
        plt.plot(gnc_struct.count, gnc_struct.del_rud_deg, 'g*', label='$\delta_{Rudder}$')
        plt.legend()
        plt.ylabel(r"$^\circ's of Actuation$")
        plt.xlabel("Simulation Number")

        plt.figure(4)
        plt.clf()
        plt.subplot(2,1,1)
        plt.title('Longitudinal Metrics, $T_{1/2}$ Amplitude')
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

if __name__ == "__main__":
    plot_gnc_object(gnc_struct)
