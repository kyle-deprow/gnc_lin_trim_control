#! /usr/bin/python3

import numpy as np
from scipy import interpolate
from scipy import signal

class F16_Model():
    def __init__(self):
    #########################################################################################
        # Version history:
            # Nov 9, 2018    initial release, inspired by .m code written by Dr. Buckholtz
        # Data taken from:
            # Aircraft Control and Simulation, 2nd Edition
            # by Stevens and Lewis
    #########################################################################################
        # Mass Parameters of f16
        self.axx = 9496.0
        self.ayy = 55814.0
        self.azz = 63100.0
        self.axz = 982.0
        self.weight_lbs = 20490.466
        self.s_ft2 = 300
        self.b_ft = 30
        self.cbar_ft = 11.32
        self.xcg = 0.3
        self.xcgr = 0.35
        self.hx_slgft2 = 160.0

        # Setup tables and interp objects
        self._f16_tables_and_interp()

        # Precomputation parameters
        self.g_fps2 = 32.174
        self.mass_slg = self.weight_lbs/self.g_fps2
        self.r2d = 180/np.pi
        self.d2r = 1/self.r2d
        self.axzs = self.axz**2
        self.xpq = self.axz*(self.axx - self.ayy + self.azz)
        self.gam = self.axx*self.azz - self.axz**2
        self.xqr = self.azz*(self.azz - self.ayy) + self.axzs
        self.zpq = self.axx*(self.axx - self.ayy) + self.axzs
        self.ypr = self.azz - self.axx

        # Presetting Dictionary names
        self.x_names = ['Airspeed',
                        'AOA',
                        'AOS',
                        'Roll',
                        'Pitch',
                        'Yaw',
                        'Roll_Rate',
                        'Pitch_Rate',
                        'Yaw_Rate',
                        'North_Position',
                        'East_Position',
                        'Altitude',
                        'Power']

        self.xd_names = ['Airspeed_dot',
                         'alpha_dot',
                         'beta_dot',
                         'Roll_Rate',
                         'Pitch_Rate',
                         'Yaw_Rate',
                         'Roll_Rate_dot',
                         'Pitch_Rate_dot',
                         'Yaw_Rate_dot',
                         'North_Position_dot',
                         'East_Position_dot',
                         'Altitude_dot',
                         'Power_dot']

        self.u_names = ['Throttle',
                        'Elevator',
                        'Aileron',
                        'Rudder']

        self.long_ind = np.array([0,1,4,7])
        self.lat_ind = np.array([2,3,6,8])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def nonlinear_model(self, time, x, u):
    #########################################################################################
        # Nonlinear 6DOF Aircraft model

        # Inputs:
            # time            Scalar instance of time
            # x               n element dictionary OR vector, state vector
                # x(0)... True Airspeed, Vt (ft/sec)
                # x(1)... Angle of Attack, alpha (rad)
                # x(2)... Angle of Sideslip, beta (rad)
                # x(3)... Roll attitude, phi (rad)
                # x(4)... Pitch attitude, theta (rad)
                # x(5)... Yaw attitude, psi (rad)
                # x(6)... Roll rate, P (rad/sec)
                # x(7)... Pitch rate, Q (rad/sec)
                # x(8)... Yaw rate, R (rad/sec)
                # x(9)... North Position, N (ft)
                # x(10).. East Position, E (ft)
                # x(11).. Altitude, h (ft)
                # x(12).. power, (W)
            # u               m element dictionary OR vector, control vector
                # u(0)... throttle, throttle (0-1)
                # u(1)... Elevator, del_el (deg)
                # u(2)... Aileron, del_ail (deg)
                # u(3)... Rudder, del_rud (deg)

        # Outputs:
            # xd              n element vector, state derivative vector
            # Az              Scalar, vertical acceleration (g)
            # Ay              Scalar, lateral acceleration (g)
    #########################################################################################

        if type(x) == np.ndarray:
            x,u = self._x_u_vector_to_dict(x,u)

        # Assign states to local variables. Convert from radians to degrees for AOA and AOS
        Vt_fps = x['Airspeed']
        alpha_rad = x['AOA']
        alpha_deg = x['AOA']*self.r2d
        beta_rad = x['AOS']
        beta_deg = x['AOS']*self.r2d
        phi_rad = x['Roll']
        theta_rad = x['Pitch']
        psi_rad = x['Yaw']
        P_rps = x['Roll_Rate']
        Q_rps = x['Pitch_Rate']
        R_rps = x['Yaw_Rate']
        alt_ft = x['Altitude']
        power = x['Power']

        # Assign controls to local variables. Leave surfaces in terms of degrees
        throttle = u['Throttle']
        del_el_deg = u['Elevator']
        del_ail_deg = u['Aileron']
        del_rud_deg = u['Rudder']


        # Preallocate output memory
        xd = [None]*len(x)

        # Compute air data parameters
        mach, qbar_psf = self.adc(Vt_fps, alt_ft)

        # Extract data from the engine model
        cpow = self.f16_tgear(throttle)
        xd[12] = self.f16_del_power(power,cpow)
        T_lbs = self.f16_thrust(power, alt_ft, mach)

        # Aerodynamic coefficient lookups
        cxt = self.f16_cx(alpha_deg, del_el_deg)
        cyt = self.f16_cy(beta_deg, del_ail_deg, del_rud_deg)
        czt = self.f16_cz(alpha_deg, beta_deg, del_el_deg)

        clt = self.f16_cl(alpha_deg, beta_deg) + \
              self.f16_cl_dail(alpha_deg, beta_deg)*(del_ail_deg/20.) + \
              self.f16_cl_drud(alpha_deg, beta_deg)*(del_rud_deg/30.)

        cmt = self.f16_cm(alpha_deg, del_el_deg)

        cnt = self.f16_cn(alpha_deg, beta_deg) + \
              self.f16_cn_dail(alpha_deg, beta_deg)*(del_ail_deg/20.) + \
              self.f16_cn_drud(alpha_deg, beta_deg)*(del_rud_deg/30.)

        # Damping Derivatives
        tvt_s = 0.5/Vt_fps
        b2v = self.b_ft*tvt_s
        cq = self.cbar_ft*Q_rps*tvt_s

        CXq, CYr, CYp, CZq, Clr, Clp, Cmq, Cnr, Cnp  = self.f16_damp(alpha_deg)
        damping  = self.f16_damp(alpha_deg)
        cxt += cq*CXq
        cyt += b2v*(CYr*R_rps + CYp*P_rps)
        czt += cq*CZq
        clt += b2v*(Clr*R_rps + Clp*P_rps)
        cmt += cq*Cmq + czt * (self.xcgr - self.xcg)
        cnt += b2v*(Cnr*R_rps + Cnp*P_rps) - \
                cyt*(self.xcgr - self.xcg)*(self.cbar_ft/self.b_ft)

        # Pre-compute some variables for the state space model.  Computations for
        # AOA and AOS must be in radians.
        salp = np.sin(alpha_rad)
        calp = np.cos(alpha_rad)
        sbta = np.sin(beta_rad)
        cbta = np.cos(beta_rad)
        sth = np.sin(theta_rad)
        cth = np.cos(theta_rad)
        sph = np.sin(phi_rad)
        cph = np.cos(phi_rad)
        spsi= np.sin(psi_rad)
        cpsi= np.cos(psi_rad)

        qs = qbar_psf*self.s_ft2
        qsb = qs*self.b_ft
        rmqs = qs/self.mass_slg
        gcth = self.g_fps2*cth
        qsph = Q_rps*sph
        Ay_fps2 = rmqs*cyt
        Az_fps2 = rmqs*czt

        # Velocities along wind axes.
        U_fps = Vt_fps*calp*cbta
        V_fps = Vt_fps*sbta
        W_fps = Vt_fps*salp*cbta

        # Force equations.
        udot_fps2 = R_rps*V_fps - Q_rps*W_fps - self.g_fps2*sth + \
                    (qs*cxt + T_lbs) / self.mass_slg
        vdot_fps2 = P_rps*W_fps - R_rps*U_fps + gcth*sph + Ay_fps2
        wdot_fps2 = Q_rps*U_fps - P_rps*V_fps + gcth*cph + Az_fps2
        dum = (np.power(U_fps,2) + np.power(W_fps,2))

        xd[0] = (U_fps*udot_fps2 + V_fps*vdot_fps2 + W_fps*wdot_fps2) / Vt_fps
        xd[1] = (U_fps*wdot_fps2 - W_fps*udot_fps2) / dum
        xd[2] = (Vt_fps*vdot_fps2 - V_fps*xd[0])*cbta / dum

        # Kinematics.
        xd[3] = P_rps + (sth/cth)*(qsph + R_rps*cph)
        xd[4] = Q_rps*cph - R_rps*sph
        xd[5] = (qsph + R_rps*cph) / cth

        # Moments.
        roll = qsb*clt
        pitch = qs*self.cbar_ft*cmt
        yaw = qsb*cnt
        PQ = P_rps*Q_rps
        QR = Q_rps*R_rps
        PR = P_rps*R_rps
        QHX = Q_rps*self.hx_slgft2

        xd[6] = (self.xpq*PQ - self.xqr*QR + self.azz*roll + self.axz*(yaw + QHX)) / self.gam
        xd[7] = (self.ypr*PR - self.axz*(np.power(P_rps,2) - np.power(R_rps,2)) + pitch - R_rps*self.hx_slgft2) / self.ayy
        xd[8] = (self.zpq*PQ - self.xpq*QR + self.axz*roll + self.axx*(yaw + QHX)) / self.gam

        # Navigation.
        t1 = sph*cpsi
        t2 = cph*sth
        t3 = sph*spsi
        s1 = cth*cpsi
        s2 = cth*spsi
        s3 = t1*sth - cph*spsi
        s4 = t3*sth + cph*cpsi
        s5 = sph*cth
        s6 = t2*cpsi + t3
        s7 = t2*spsi - t1
        s8 = cph*cth

        xd[9] = U_fps*s1 + V_fps*s3 + W_fps*s6
        xd[10] = U_fps*s2 + V_fps*s4 + W_fps*s7
        xd[11] = U_fps*sth - V_fps*s5 - W_fps*s8

        xd_dict = dict(zip(self.xd_names,xd))

        # Set the outputs of the vertical and lateral accelerations, in g's
        Az = Az_fps2/self.g_fps2
        Ay = Ay_fps2/self.g_fps2

        return xd_dict, Az, Ay

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def build_x_u(self, ux0, input_conditions):
    #########################################################################################
        # x               n element vector, state vector
            # x(0)... True Airspeed, Vt (ft/sec)
            # x(1)... Angle of Attack, alpha (rad)
            # x(2)... Angle of Sideslip, beta (rad)
            # x(3)... Roll attitude, phi (rad)
            # x(4)... Pitch attitude, theta (rad)
            # x(5)... Yaw attitude, psi (rad)
            # x(6)... Roll rate, P (rad/sec)
            # x(7)... Pitch rate, Q (rad/sec)
            # x(8)... Yaw rate, R (rad/sec)
            # x(9)... North Position, N (ft)
            # x(10).. East Position, E (ft)
            # x(11).. Altitude, h (ft)
            # x(12).. power, (W)
        # u               m element vector, control vector
            # u(0)... throttle, throttle (0-1)
            # u(1)... Elevator, del_el (deg)
            # u(2)... Aileron, del_ail (deg)
            # u(3)... Rudder, del_rud (deg)

        # IMPORTANT -- FUNCTION RETURNS A DICTIONARY RESPONDING TO STATE, CONTROL VECTOR
    #########################################################################################
        # compute power
        power = self.f16_tgear(ux0[4])

        # build the state vector:
        x = [input_conditions[0],
             input_conditions[1],
             ux0[0],
             0,
             0,
             0,
             0,
             0,
             0,
             0,
             0,
             input_conditions[2],
             power]

        x_dict = dict(zip(self.x_names,x))

        # build the control vector:
        u = [ux0[4],                # throttle
             ux0[1],                # del_el
             ux0[3],                # del_aileron
             ux0[2]]                # del_rudder

        u_dict = dict(zip(self.u_names,u))
        return x_dict, u_dict

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def build_state_matrices(self, a, b):
        Alon = a[self.long_ind,:][:,self.long_ind]
        Blon = b[self.long_ind,:][0]
        Alat = a[self.lat_ind,:][:,self.lat_ind]
        Blat = b[self.long_ind,:][1]
        return Alon, Blon, Alat, Blat

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def adc(self, v_fps, h_ft):
    #########################################################################################
        # This function estimates the mach number and dynamic pressure for a
        # particular flight condition.

        # Inputs:
            # v_fps           scalar, velocity (fps)
            # h_ft            scalar, altitude (ft)

        # Outputs:
            # mach            scalar, mach number
            # qbar_psf        scalar, dynamic pressure (psf)
    #########################################################################################
        # Parameters
        # Sea level density, R0
        R0_cf = 2.377e-3

        # Temperature estimation
        tfac = 1.0 - 0.703e-5 * h_ft
        if h_ft >= 35000.:
            T = 290
        else:
            T = 519*tfac

        # Estimate the air density at h
        rho_cf = R0_cf * (np.power(tfac,4.14))

        # Estimate Mach number
        mach = v_fps/np.sqrt(1.4*1716.3*T)

        # Estimate dynamic pressure
        qbar_psf = 0.5 * rho_cf * np.power(v_fps,2)

        return mach, qbar_psf

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_thrust(self, power, alt_ft, mach):
    #########################################################################################
        # F-16 data and interpolation

        # Inputs:
            # power         percentage of power
            # alt_ft        aircraft altitude (ft)
            # mach          Mach number

        # Outputs:
            # ti_lbs        thrust (lbs)
    #########################################################################################
        # Compute thrust from the military table
        milthr = self.mil_interp(alt_ft, mach)[0]
        if power < 50:
            # Compute thrust between idle and military tables
            idlethr = self.idle_interp(alt_ft, mach)[0]
            ti_lbs = idlethr + (milthr - idlethr)*power*0.02
        else:
            # Compute thrust between military and max tables
            maxthr = self.maxtab_interp(alt_ft, mach)[0]
            ti_lbs = milthr + (maxthr - milthr)*(power-50)*0.02
        return ti_lbs

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_cx(self, alpha_deg, elev_deg):
    #########################################################################################
        # F-16 axial force coefficient

        # Inputs:
            # alpha_deg           angle of attack (deg)
            # elev_deg            elevator deflection (deg)

        # Outputs:
            # cxi                 axial force coefficient
    #########################################################################################
        # Bound inputs to saturation ranges
        alpha_ltd_deg = max(-10, min(45, alpha_deg))
        elev_ltd_deg = max(-24, min(24, elev_deg))
        return self.cx_interp(alpha_ltd_deg, elev_ltd_deg)[0]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_cy(self, beta_deg, ail_deg, rud_deg):
    #########################################################################################
        # F-16 y-force coefficient

        # Inputs:
            # beta_deg            sideslip angle (deg)
            # ail_deg             aileron deflection (deg)
            # rud_deg             rudder deflection (deg)

        # Outputs:
            # cyi                 lateral force coefficient
    #########################################################################################
        # Limit inputs to saturation range
        ail_ltd_deg = min(20, max(-20, ail_deg))
        rud_ltd_deg = min(30, max(-30, rud_deg))

        # Compute the lateral force coefficient
        return -0.02*beta_deg + 0.021*(ail_ltd_deg/20.) + 0.086*(rud_ltd_deg/30.)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_cz(self, alpha_deg, beta_deg, elev_deg):
    #########################################################################################
        # F-16 z-force coefficient

        # Inputs:
            # alpha_deg           angle of attack (deg)
            # beta_deg            sideslip angle (deg)
            # elev_deg            elevator deflection (deg)

        # Outputs:
            # czi                 Z force coefficient
    #########################################################################################
        # Limit inputs to saturation range
        alpha_ltd_deg = min(45, max(-10, alpha_deg))
        elev_ltd_deg = min(25, max(-25, elev_deg))
        c = self.cz_interp(alpha_ltd_deg)

        # Compute the Z force coefficient
        return c*(1 - (beta_deg*self.d2r)**2) - 0.19*(elev_ltd_deg/25)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_cl(self, alpha_deg, beta_deg):
    #########################################################################################
        # F-16 Rolling Moment Coefficient Interpolation

        # Inputs:
            # alpha_deg           angle of attack (deg)
            # beta_deg            sideslip angle (deg)

        # Outputs:
            # cli                 Rolling Moment coefficient
    #########################################################################################
        # Limit inputs to saturation range
        alpha_ltd_deg = max(-10, min(45, alpha_deg))
        beta_ltd_deg = max(-30, min(30, beta_deg))

        # Interpretation
        cli = self.cl_interp(alpha_ltd_deg, abs(beta_ltd_deg))[0]

        # Account for the sign of the sideslip
        if beta_deg < 0.0:
            cli = -cli
        return cli

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_cl_dail(self, alpha_deg, beta_deg):
    #########################################################################################
        # F-16 Rolling Moment Coefficient due to Aileron delfection Interpolation

        # Inputs:
            # alpha_deg           angle of attack (deg)
            # beta_deg            sideslip angle (deg)

        # Outputs:
            # clda                Rolling Moment due to aileron coefficient
    #########################################################################################
        # Limit inputs to saturation range
        alpha_ltd_deg = max(-10, min(45, alpha_deg))
        beta_ltd_deg = max(-15, min(15, beta_deg))

        # Interpretation
        return self.cldail_interp(alpha_ltd_deg, beta_ltd_deg)[0]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_cl_drud(self, alpha_deg, beta_deg):
    #########################################################################################
        # F-16 Rolling Moment Coefficient due to rudder delfection Interpolation

        # Inputs:
            # alpha_deg           angle of attack (deg)
            # beta_deg            sideslip angle (deg)

        # Outputs:
            # cldr                Rolling Moment due to rudder coefficient
    #########################################################################################
        # Limit inputs to saturation range
        alpha_ltd_deg = max(-10, min(45, alpha_deg))
        beta_ltd_deg = max(-15, min(15, beta_deg))

        # Interpretation
        return self.cldrud_interp(alpha_ltd_deg, beta_ltd_deg)[0]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_cm(self, alpha_deg, elev_deg):
    #########################################################################################
        # F-16 Rolling Pitching Moemnt Coefficient

        # Inputs:
            # alpha_deg           angle of attack (deg)
            # elev_deg            elevator angle (deg)

        # Outputs:
            # cmi                 Pitching Moment Coeffecient
    #########################################################################################
        # Limit inputs to saturation range
        alpha_ltd_deg = max(-10, min(45, alpha_deg))
        elev_ltd_deg = max(-24, min(24, elev_deg))

        # Interpretation
        return self.cm_interp(alpha_ltd_deg, elev_ltd_deg)[0]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_cn(self, alpha_deg, beta_deg):
    #########################################################################################
        # F-16 Yawing Pitching Moemnt Coefficient

        # Inputs:
            # alpha_deg           angle of attack (deg)
            # beta_deg            sideslip angle (deg)

        # Outputs:
            # cni                 Yawing Moment Coefficient
    #########################################################################################
        # Limit inputs to saturation range
        alpha_ltd_deg = max(-10, min(45, alpha_deg))
        beta_ltd_deg = max(-30, min(30, beta_deg))

        # Interpretation
        cni = self.cn_interp(alpha_ltd_deg, abs(beta_ltd_deg))[0]
        if beta_deg < 0.0:
            cni = -cni
        return cni

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_cn_dail(self, alpha_deg, beta_deg):
    #########################################################################################
        # F-16 Yawing Moment Coefficient due to Aileron delfection Interpolation

        # Inputs:
            # alpha_deg           angle of attack (deg)
            # beta_deg            sideslip angle (deg)

        # Outputs:
            # clna                Rolling Moment due to aileron coefficient
    #########################################################################################
        # Limit inputs to saturation range
        alpha_ltd_deg = max(-10, min(45, alpha_deg))
        beta_ltd_deg = max(-15, min(15, beta_deg))

        # Interpretation
        return self.cnda_interp(alpha_ltd_deg, beta_ltd_deg)[0]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_cn_drud(self, alpha_deg, beta_deg):
    #########################################################################################
        # F-16 Yawing Moment Coefficient due to rudder delfection Interpolation

        # Inputs:
            # alpha_deg           angle of attack (deg)
            # beta_deg            sideslip angle (deg)

        # Outputs:
            # cndr                Rolling Moment due to rudder coefficient
    #########################################################################################
        # Limit inputs to saturation range
        alpha_ltd_deg = max(-10, min(45, alpha_deg))
        beta_ltd_deg = max(-15, min(15, beta_deg))

        # Interpretation
        return self.cndr_interp(alpha_ltd_deg, beta_ltd_deg)[0]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_damp(self, alpha_deg):
    #########################################################################################
        # Function to Compute the various damping derivatives for the f-16

        # Inputs:
            # alpha_deg               angle of attack (deg)

        # Outputs:
            # damp_coeff              np.array of computed damping coefficients
                      # [0]--CXq
                      # [1]--CYr
                      # [2]--CYp
                      # [3]--CZq
                      # [4]--Clr
                      # [5]--Clp
                      # [6]--Cmq
                      # [7]--Cnr
                      # [8]--Cnp
    #########################################################################################
        damp_coeff = np.zeros(9)
        for i in range(len(damp_coeff)):
            damp_coeff[i] = self.damp_interp[i](alpha_deg)

        return damp_coeff

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_del_power(self, p3, p1):
    #########################################################################################
        # F16 Power rate of change

        # Inputs:
            # p3          Scalar, actual power
            # p1          Scalar, power command

        # Outputs:
            # pdot        Scalar, Power Change
    #########################################################################################
        if p1 >= 50:
            if p3 >= 50:
                p2 = p1
                T = 5.0
            else:
                p2 = 60.0
                T = self.rtau(p2-p3)
        else:
            if p3 >= 50:
                p2 = 40.0
                T = 5.0
            else:
                p2 = p1
                T = self.rtau((p2-p3))

        return T*(p2-p3)

    def rtau(self, delta_p):
    #########################################################################################
        # Computes reciprocal time constant

        # Inputs:
            # delta_p       Scalar, difference in power

        # Outputs:
            # T             Scalar, reciprocal time constant
    #########################################################################################
        if delta_p <= 25:
            T = 1.0;
        elif delta_p >= 50:
            T = 0.1
        else:
            T = 1.9 - 0.036*delta_p
        return T

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def f16_tgear(self, throttle):
    #########################################################################################
        # Data and interpolation for F-16 power command vs throttle

        # Inputs:
            # throttle        scalar, throttle position (decimal percentage)

        # Outputs:
            # tgear           scalar, power command
    #########################################################################################
        if throttle <= 0.77:
            tgear = 64.94*throttle
        else:
            tgear = 217.38*throttle - 117.38
        return tgear

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def _f16_tables_and_interp(self):
    #########################################################################################
        # Sets up F16 tables and necessary interpolation objects
    #########################################################################################
        # Setup breakpoints
        htabl_ft = np.arange(0,6e4,1e4)
        machtabl = np.arange(0,1.2,0.2)

        """From F16_thrust.m"""
        # Thrust tables depending upon power.
        idle = [[1060, 670, 880, 1140, 1500, 1860],
                [635, 425, 690, 1010, 1330, 1700],
                [60, 25, 345, 755, 1130, 1525],
                [-1020, -710, -300, 350, 910, 1360],
                [-2700, -1900, -1300, -247, 600, 1100],
                [-3600, -1400, -595, -342, -200, 700]]
        idle = np.array(idle)

        mil = [[12680, 9150, 6200, 3950, 2450, 1400],
               [12680, 9150, 6313, 4040, 2470, 1400],
               [12610, 9312, 6610, 4290, 2600, 1560],
               [12640, 9839, 7090, 4660, 2840, 1660],
               [12390, 10176, 7750, 5320, 3250, 1930],
               [11680, 9848, 8050, 6100, 3800, 2310]]
        mil = np.array(mil)

        maxtab = [[20000, 15000, 10800, 7000, 4000, 2500],
                  [21420, 15700, 11225, 7323, 4435, 2600],
                  [22700, 16860, 12250, 8154, 5000, 2835],
                  [24240, 18910, 13760, 9285, 5700, 3215],
                  [26070, 21075, 15975, 11115, 6860, 3950],
                  [28886, 23319, 18300, 13484, 8642, 5057]]
        maxtab = np.array(maxtab)

        # Set-up mil, idle, and maxtab interpolaters
        self.mil_interp = interpolate.interp2d(htabl_ft, machtabl, mil)
        self.idle_interp = interpolate.interp2d(htabl_ft, machtabl, idle)
        self.maxtab_interp = interpolate.interp2d(htabl_ft, machtabl, maxtab)

        """From F16_cx.m"""
        alphtab_deg = np.arange(-10,50,5)
        elevtab_deg = np.arange(-24,36,12)
        # CX table
        cxtabl = [[-0.099, -0.081, -0.081, -0.063, -0.025, 0.044, 0.097, 0.113, 0.145, 0.167, 0.174, 0.166],
                  [-0.048, -0.038, -0.040, -0.021, 0.016, 0.083, 0.127, 0.137, 0.162, 0.177, 0.179, 0.167],
                  [-0.022, -0.020, -0.021, -0.004, 0.032, 0.094, 0.128, 0.130, 0.154, 0.161, 0.155, 0.138],
                  [-0.040, -0.038, -0.039, -0.025, 0.006, 0.062, 0.087, 0.085, 0.100, 0.110, 0.104, 0.091],
                  [-0.083, -0.073, -0.076, -0.072, -0.046, 0.012, 0.024, 0.025, 0.043, 0.053, 0.047, 0.040]]
        cxtabl = np.array(cxtabl)
        # Set-up cx interp object
        self.cx_interp = interpolate.interp2d(alphtab_deg, elevtab_deg, cxtabl)

        """From F16_cz.m"""
        alphtab_deg = np.arange(-10,50,5)
        cztabl = [0.770, 0.241, -0.100, -0.416, -0.731, -1.503, -1.366, -1.646, -1.917, -2.120, -2.248, -2.229]
        cztabl = np.array(cztabl)
        # Create cz interp object
        self.cz_interp = interpolate.interp1d(alphtab_deg, cztabl)

        """From F16_cl.m"""
        alphtab_deg = np.arange(-10,50,5)
        betatab_deg = np.arange(0,35,5)
        cltabl = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [-0.001, -0.004, -0.008, -0.012, -0.016, -0.019, -0.020, -0.020, -0.015, -0.008, -0.013, -0.015],
                  [-0.003, -0.009, -0.017, -0.024, -0.030, -0.034, -0.040, -0.037, -0.016, -0.002, -0.010, -0.019],
                  [-0.001, -0.010, -0.020, -0.030, -0.039, -0.044, -0.050, -0.049, -0.023, -0.006, -0.014, -0.027],
                  [0.000, -0.010, -0.022, -0.034, -0.047, -0.046, -0.059, -0.061, -0.033, -0.036, -0.035, -0.035],
                  [0.007, -0.010, -0.023, -0.034, -0.049, -0.046, -0.068, -0.071, -0.060, -0.058, -0.062, -0.059],
                  [0.009, -0.011, -0.023, -0.037, -0.050, -0.047, -0.074, -0.079, -0.091, -0.076, -0.077, -0.076]]
        cltabl = np.array(cltabl)
        # Set-up cl interp object
        self.cl_interp = interpolate.interp2d(alphtab_deg, betatab_deg, cltabl)

        """From F16_cl_dail.m"""
        alphtab_deg = np.arange(-10,50,5)
        betatab_deg = np.arange(-15,20,5)
        clda_tab = [[-0.041, -0.052, -0.053, -0.056, -0.050, -0.056, -0.082, -0.059, -0.042, -0.038, -0.027, -0.017],
                    [-0.041, -0.053, -0.053, -0.053, -0.050, -0.051, -0.066, -0.043, -0.038, -0.027, -0.023, -0.016],
                    [-0.042, -0.053, -0.052, -0.051, -0.049, -0.049, -0.043, -0.035, -0.026, -0.016, -0.018, -0.014],
                    [-0.040, -0.052, -0.051, -0.052, -0.048, -0.048, -0.042, -0.037, -0.031, -0.026, -0.017, -0.012],
                    [-0.043, -0.049, -0.048, -0.049, -0.043, -0.042, -0.042, -0.036, -0.025, -0.021, -0.016, -0.011],
                    [-0.044, -0.048, -0.048, -0.047, -0.042, -0.041, -0.020, -0.028, -0.013, -0.014, -0.011, -0.010],
                    [-0.043, -0.049, -0.047, -0.045, -0.042, -0.037, -0.003, -0.013, -0.010, -0.003, -0.007, -0.008]]
        clda_tab = np.array(clda_tab)
        # Creat cldail interp object
        self.cldail_interp = interpolate.interp2d(alphtab_deg, betatab_deg, clda_tab)


        """From F16_cl_dr.m"""
        alphtab_deg = np.arange(-10,50,5)
        betatab_deg = np.arange(-15,20,5)
        cldr_tab = [[0.005,  0.017,  0.014,  0.010, -0.005,  0.009,  0.019,  0.005, -0.000, -0.005, -0.011,  0.008],
                    [0.007,  0.016,  0.014,  0.014,  0.013,  0.009,  0.012,  0.005, 0.000,  0.004,  0.009,  0.007],
                    [0.013,  0.013,  0.011,  0.012,  0.011,  0.009,  0.008,  0.005, -0.002,  0.005,  0.003,  0.005],
                    [0.018,  0.015,  0.015,  0.014,  0.014,  0.014,  0.014,  0.015, 0.013,  0.011,  0.006,  0.001],
                    [0.015,  0.014,  0.013,  0.013,  0.012,  0.011,  0.011,  0.010, 0.008,  0.008,  0.007,  0.003],
                    [0.021,  0.011,  0.010,  0.011,  0.010,  0.009,  0.008,  0.010, 0.006,  0.005,  0.000,  0.001],
                    [0.023,  0.010,  0.011,  0.011,  0.011,  0.010,  0.008,  0.010, 0.006,  0.014,  0.020,  0.000]]
        cldr_tab = np.array(cldr_tab)
        # Create cldr interp object
        self.cldrud_interp = interpolate.interp2d(alphtab_deg, betatab_deg, cldr_tab)

        """From F16_cm.m"""
        alphtab_deg = np.arange(-10,50,5)
        elevtab_deg = np.arange(-24,36,12)
        cmtabl = [[0.205,  0.168,  0.186,  0.196,  0.213,  0.251,  0.245,  0.238, 0.252,  0.231,  0.198,  0.192],
                  [0.081,  0.077,  0.107,  0.110,  0.110,  0.141,  0.127,  0.119, 0.133,  0.108,  0.081,  0.093],
                  [-0.046, -0.020, -0.009, -0.005, -0.006,  0.010,  0.006, -0.001, 0.014,  0.000, -0.013,  0.032],
                  [-0.174, -0.145, -0.121, -0.127, -0.129, -0.102, -0.097, -0.113, -0.087, -0.084, -0.069, -0.006],
                  [-0.259, -0.202, -0.184, -0.193, -0.199, -0.150, -0.160, -0.167, -0.104, -0.076, -0.041, -0.005]]
        cmtabl = np.array(cmtabl)
        # Create cm interp object
        self.cm_interp = interpolate.interp2d(alphtab_deg, elevtab_deg, cmtabl)

        """From F16_cn.m"""
        alphtab_deg = np.arange(-10,50,5)
        betatab_deg = np.arange(0,35,5)
        cntabl = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0.018,  0.019,  0.018,  0.019,  0.019,  0.018,  0.013,  0.007, 0.004, -0.014, -0.017, -0.033],
                  [0.038,  0.042,  0.042,  0.042,  0.043,  0.039,  0.030,  0.017, 0.004, -0.035, -0.047, -0.057],
                  [0.056,  0.057,  0.059,  0.058,  0.058,  0.053,  0.032,  0.012, 0.002, -0.046, -0.071, -0.073],
                  [0.064,  0.077,  0.076,  0.074,  0.073,  0.057,  0.029,  0.007, 0.012, -0.034, -0.065, -0.041],
                  [0.074,  0.086,  0.093,  0.089,  0.080,  0.062,  0.049,  0.022, 0.028, -0.012, -0.002, -0.013],
                  [0.079,  0.090,  0.106,  0.106,  0.096,  0.080,  0.068,  0.030, 0.064,  0.015,  0.011, -0.001]]
        cntabl = np.array(cntabl)
        # Create cn interp object
        self.cn_interp = interpolate.interp2d(alphtab_deg, betatab_deg, cntabl)

        """From F16_cn_dail.m"""
        alphtab_deg = np.arange(-10,50,5)
        betatab_deg = np.arange(-15,20,5)

        cnda_tab = [[0.001, -0.027, -0.017, -0.013, -0.012, -0.016,  0.001,  0.017, 0.011,  0.017,  0.008,  0.016],
                    [0.002, -0.014, -0.016, -0.016, -0.014, -0.019, -0.021,  0.002, 0.012,  0.015,  0.015,  0.011],
                    [-0.006, -0.008, -0.006, -0.006, -0.005, -0.008, -0.005,  0.007, 0.004,  0.007,  0.006,  0.006],
                    [-0.011, -0.011, -0.010, -0.009, -0.008, -0.006,  0.000,  0.004, 0.007,  0.010,  0.004,  0.010],
                    [-0.015, -0.015, -0.014, -0.012, -0.011, -0.008, -0.002,  0.002, 0.006,  0.012,  0.011,  0.011],
                    [-0.024, -0.010, -0.004, -0.002, -0.001,  0.003,  0.014,  0.006, -0.001,  0.004,  0.004,  0.006],
                    [-0.022,  0.002, -0.003, -0.005, -0.003, -0.001, -0.009, -0.009, -0.001,  0.003, -0.002,  0.001]]
        cnda_tab = np.array(cnda_tab)
        # Create cnda interp object
        self.cnda_interp = interpolate.interp2d(alphtab_deg, betatab_deg, cnda_tab)

        """From F16_cn_dr.m"""
        alphtab_deg = np.arange(-10,50,5)
        betatab_deg = np.arange(-15,20,5)

        cndr_tab = [[-0.018, -0.052, -0.052, -0.052, -0.054, -0.049, -0.059, -0.051, -0.030, -0.037, -0.026, -0.013],
                    [-0.028, -0.051, -0.043, -0.046, -0.045, -0.049, -0.057, -0.052, -0.030, -0.033, -0.030, -0.008],
                    [-0.037, -0.041, -0.038, -0.040, -0.040, -0.038, -0.037, -0.030, -0.027, -0.024, -0.019, -0.013],
                    [-0.048, -0.045, -0.045, -0.045, -0.044, -0.045, -0.047, -0.048, -0.049, -0.045, -0.033, -0.016],
                    [-0.043, -0.044, -0.041, -0.041, -0.040, -0.038, -0.034, -0.035, -0.035, -0.029, -0.022, -0.009],
                    [-0.052, -0.034, -0.036, -0.036, -0.035, -0.028, -0.024, -0.023, -0.020, -0.016, -0.010, -0.014],
                    [-0.062, -0.034, -0.027, -0.028, -0.027, -0.027, -0.023, -0.023, -0.019, -0.009, -0.025, -0.010]]
        cndr_tab = np.array(cndr_tab)
        # Create cndr interp object
        self.cndr_interp = interpolate.interp2d(alphtab_deg, betatab_deg, cndr_tab)

        """From F16_damp.m"""
        alphtab_deg = np.arange(-10,50,5)

        damp_tab = [[-0.267, -0.110, 0.308, 1.340, 2.080, 2.910, 2.760, 2.050, 1.500, 1.490, 1.830, 1.210],
                    [0.882,  0.852, 0.876, 0.958, 0.962,  0.974,  0.819, 0.483, 0.590, 1.210, -0.493, -1.040],
                    [-0.108, -0.108, -1.880, 0.110, 0.258, 0.226, 0.344, 0.362, 0.611, 0.529, 0.298, -2.270],
                    [-8.800,-25.800,-28.900,-31.400,-31.200,-30.700,-27.700, -28.200, -29.000,-29.800,-38.300,-35.300],
                    [-0.126, -0.026, 0.063, 0.113,  0.208, 0.230, 0.319, 0.437, 0.680, 0.100, 0.447, -0.330],
                    [-0.360, -0.359, -0.443, -0.420, -0.383, -0.375, -0.329, -0.294, -0.230, -0.210, -0.120, -0.100],
                    [-7.210, -0.540, -5.230, -5.260, -6.110, -6.640, -5.690, -6.000, -6.200, -6.400, -6.600, -6.000],
                    [-0.380, -0.363, -0.378, -0.386, -0.370, -0.453, -0.550, -0.582, -0.595, -0.637, -1.020, -0.840],
                    [0.061,  0.052,  0.052, -0.012, -0.013, -0.024,  0.050,  0.150, 0.130,  0.158,  0.240,  0.150]]
        damp_tab = np.transpose(np.array(damp_tab))
        # Create damp interp objects
        self.damp_interp = [None]*9
        for i in range(9):
            self.damp_interp[i] = interpolate.interp1d(alphtab_deg, damp_tab[:,i])

    def _x_u_vector_to_dict(self,x,u):
        x_dict = dict(zip(self.x_names,x.flatten()))
        u_dict = dict(zip(self.u_names,u.flatten()))
        return x_dict, u_dict

if __name__ == "__main__":
    f16 = F16_Model()
    ux0 = np.array([0,-1.931,0,0,0.1485])
    input_conditions = np.array([500,2*np.pi/180,1])
    x,u = f16.build_x_u(ux0,input_conditions)
    print(f16.nonlinear_model(0,x,u))
