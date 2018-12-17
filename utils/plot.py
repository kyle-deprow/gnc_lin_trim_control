#! /usr/bin/env/python3

import matplotlib.pyplot as plt
import numpy as np

from colour import Color
from matplotlib.backends.backend_pdf import PdfPages
def plot_gnc_object(gnc_struct):

    pp = PdfPages('Gnc_report.pdf')

    altitudes = np.unique(gnc_struct.altitude)
    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    for altitude in altitudes:
        ind = np.where(gnc_struct.altitude==altitude)[0]
        count = gnc_struct.count[ind]
        number = np.arange(len(count))

        f, axarr = plt.subplots(4,1)
        axarr[0].set_title('Trim Data for Alt: '+str(altitude)+' ft')
        axarr[0].plot(number, gnc_struct.alpha_rad[ind], 'bx', label='AOA')
        axarr[0].set_ylabel(r'$\alpha$ rad')
        axarr[0].set_xticklabels('')
        axarr[0].grid()
        axarr[0].yaxis.offsetText.set_fontsize(4)
        axarr[0].legend()
        axarr[1].plot(number, gnc_struct.V_fps[ind], 'rx', label='Airspeed')
        axarr[1].set_ylabel('V ft/s')
        axarr[1].set_xticklabels('')
        axarr[1].yaxis.offsetText.set_fontsize(4)
        axarr[1].legend()
        axarr[1].grid()
        axarr[2].plot(number, gnc_struct.beta_rad[ind], 'mx', label='AOS')
        axarr[2].set_ylabel(r'$\beta$ rad')
        axarr[2].set_xticklabels('')
        axarr[2].yaxis.offsetText.set_fontsize(4)
        axarr[2].legend()
        axarr[2].grid()
        axarr[3].plot(number, gnc_struct.Nz[ind], 'gx', label='Nz')
        axarr[3].set_ylabel("acceleration g's")
        axarr[3].set_xlabel("Simulation Number")
        axarr[3].yaxis.offsetText.set_fontsize(4)
        axarr[3].legend()
        axarr[3].grid()

        f.savefig('./results/'+str(altitude)+'ft/input_conditions.png',bbox_inches='tight',dpi=300)
        plt.close()

        f, axarr = plt.subplots(4,1)
        axarr[0].set_title('Trim Data for Alt: '+str(altitude)+' ft')
        axarr[0].plot(number, gnc_struct.alpha_dot[ind], 'bx', label=r'$\dot \alpha$')
        axarr[0].set_ylabel(r'rad/s')
        axarr[0].yaxis.offsetText.set_fontsize(4)
        axarr[0].legend()
        axarr[0].grid()
        axarr[0].set_xticklabels('')
        axarr[1].plot(number, gnc_struct.beta_dot[ind], 'rx', label=r'$\dot \beta$')
        axarr[1].set_ylabel(r'rad/s')
        axarr[1].set_xticklabels('')
        axarr[1].yaxis.offsetText.set_fontsize(4)
        axarr[1].legend()
        axarr[1].grid()
        axarr[2].plot(number, gnc_struct.roll_rate_dot[ind], 'gx', label=r'$\dot p$')
        axarr[2].set_ylabel(r'rad/s')
        axarr[2].set_xticklabels('')
        axarr[2].yaxis.offsetText.set_fontsize(4)
        axarr[2].legend()
        axarr[2].grid()
        axarr[3].plot(number, gnc_struct.pitch_rate_dot[ind], 'mx', label=r'$\dot q$')
        axarr[3].set_ylabel(r'rad/s')
        axarr[3].set_xlabel(r'Simulation Number')
        axarr[3].yaxis.offsetText.set_fontsize(4)
        axarr[3].legend()
        axarr[3].grid()

        f.savefig('./results/'+str(altitude)+'ft/trim_data.png',bbox_inches='tight',dpi=300)
        plt.close()

        f, axarr = plt.subplots(1,1)
        axarr.set_title('Elevator, Aileron, and Rudder Positions for '+str(altitude)+' ft')
        axarr.plot(number, gnc_struct.del_el_deg[ind], 'bx', label='$\delta_{Elevator}$')
        axarr.plot(number, gnc_struct.del_ail_deg[ind], 'ro', label='$\delta_{Aileron}$')
        axarr.plot(number, gnc_struct.del_rud_deg[ind], 'g*', label='$\delta_{Rudder}$')
        axarr.yaxis.offsetText.set_fontsize(4)
        axarr.legend()
        axarr.grid()
        axarr.set_ylabel(r"$^\circ's$ of Actuation")
        axarr.set_xlabel("Simulation Number")

        f.savefig('./results/'+str(altitude)+'ft/trim_control_surfaces.png',bbox_inches='tight',dpi=300)
        plt.close()

        f, axarr = plt.subplots(2,1)
        axarr[0].set_title('Longitudinal Metrics, $T_{1/2}$ Amplitude for '+str(altitude)+' ft')
        axarr[0].plot(number, gnc_struct.long_half_amp[:,0][ind], 'bx', label='Short Period')
        axarr[0].yaxis.offsetText.set_fontsize(4)
        axarr[0].legend()
        axarr[0].set_ylabel('Time (s)')
        axarr[0].grid()
        axarr[1].plot(number, gnc_struct.long_half_amp[:,1][ind], 'rx', label='Phugoid')
        axarr[1].yaxis.offsetText.set_fontsize(4)
        axarr[1].legend()
        axarr[1].set_ylabel('Time (s)')
        axarr[1].set_xlabel('Simulation Number')
        axarr[1].grid()

        f.savefig('./results/'+str(altitude)+'ft/linearize_longdata.png',bbox_inches='tight',dpi=300)
        plt.close()

        f, axarr = plt.subplots(3,1)
        axarr[0].set_title('Latitudinal Metrics, $T_{1/2}$ Amplitude')
        axarr[0].plot(number, gnc_struct.lat_half_amp[:,0][ind], 'bx', label='Dutch Roll')
        axarr[0].yaxis.offsetText.set_fontsize(4)
        axarr[0].legend()
        axarr[0].grid()
        axarr[0].set_xticklabels('')
        axarr[0].set_ylabel('Time(s)')
        axarr[1].plot(number, gnc_struct.lat_half_amp[:,1][ind], 'rx', label='Roll')
        axarr[1].yaxis.offsetText.set_fontsize(4)
        axarr[1].legend()
        axarr[1].grid()
        axarr[1].set_xticklabels('')
        axarr[1].set_ylabel('Time(s)')
        axarr[2].plot(number, gnc_struct.lat_half_amp[:,2][ind], 'gx', label='Spiral')
        axarr[2].legend()
        axarr[2].set_ylabel('Time(s)')
        axarr[2].set_xlabel('Simulation Number')
        axarr[2].grid()

        f.savefig('./results/'+str(altitude)+'ft/linearize_latdata.png',bbox_inches='tight',dpi=300)
        plt.close()

        f, axarr = plt.subplots(2,1)
        axarr[0].set_title('Step and Control Responses for '+str(altitude)+' ft')
        axarr[0].set_ylabel(r'$\theta$ deg')
        axarr[0].set_xticklabels('')
        axarr[0].grid()
        axarr[1].set_ylabel(r'$\delta_{Elevator}$ deg')
        axarr[1].set_xlabel('Time(s)')
        axarr[1].grid()
        blue = Color("Blue")
        colors = list(blue.range_to(Color("red"),len(ind)))
        for ind_i,c,num in zip(ind,colors,number):
            axarr[0].plot(gnc_struct.timeoptimal_sec[ind_i],gnc_struct.yoptimal[ind_i],color=(c.rgb),label='Sim '+str(num))
            axarr[1].plot(gnc_struct.timeoptimal_sec[ind_i],gnc_struct.uoptimal[ind_i],color=(c.rgb),label='Sim '+str(num))
        axarr[0].legend(fontsize=2.5)
        axarr[1].legend(fontsize=2.5)

        f.savefig('./results/'+str(altitude)+'ft/control_data.png',bbox_inches='tight',dpi=300)
        plt.close()


