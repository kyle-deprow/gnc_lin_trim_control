#! /usr/bin/env
import numpy as np
import pprint
import control
import matplotlib.pyplot as plt

class sys_setup():
    def __init__(self):
        self.A = None
        self.B = None
        self.B2 = None
        self.C = None
        self.D = None
        self.D2 = None


class LQR_Controller():
    def __init__(self):
        pass

    def lqr_control(self,A,B,C,D,reg_array):
        n,m = B.shape
        sys_ss = control.ss(A,B,C,D)

        # Optimal Contol: Tracker
        Creg = np.zeros(C.shape[1])
        Creg[reg_array] = 1
        Dreg = sys_ss.D[2,:]
        nreg = 1

        # Build the servomechanizm state space model of the plant by appending an
        # integral error state for each regulating channel.
        Aw = np.zeros([n+nreg, n+nreg])
        # print(Aw)
        # print(Aw[0:nreg][0,1:])
        # print(Creg)
        Aw[0:nreg][0,1:] = Creg

        Aw[nreg:,nreg:] = sys_ss.A

        Bw = np.zeros([n+nreg, m])
        Bw[0:nreg,:] = Dreg
        Bw[nreg:,:] = sys_ss.B

        # # Specify range of scale factor to apply to the integral error state.
        # q = np.logspace(1,6,300)

        # % Specify the LQR matrices.
        Qw = np.zeros([n+nreg, n+nreg])
        Qw[4,4] = 3000

        # Rw = 1*np.eye(m)
        Rw = 0.01*np.eye(m)

        # % Define frequency range for performing frequency analysis.
        # w = np.logspace(-3,3,600)

        # % Form plant system.  Set output matrix for state feedback control.
        sys_p = sys_setup()
        sys_p.A = sys_ss.A
        sys_p.B = sys_ss.B
        sys_p.C = np.eye(n)
        sys_p.D = np.zeros([n,1])

        # % Pass in the characteristics of the system.  Note, this process
        # % returns plots.  This process assumes Dreg == 0.

        # Perform LQR Design # Set the integral error state weights to the desired scale factor
        for j in range(nreg):
            Qw[j,j] = 6000000

        # Solve the LQR problem.
        print(Aw)
        print(Bw)
        print(Qw)
        print(Rw)
        Kw, S, E = control.lqr(Aw, Bw, Qw, Rw)
        print(Kw)

        # Set the state feedback gain and the integral error gain.
        kI = Kw[:,:nreg][0]
        Kv = Kw[:,nreg:][0]

        # Build controller object to track reference theta
        sys_c = sys_setup
        sys_c.A = np.zeros([nreg,nreg])
        sys_c.B= Creg
        sys_c.B2= -np.eye(nreg)
        sys_c.C = -kI
        sys_c.D= -Kv
        sys_c.D2= np.zeros([m,nreg])

        # Form closed loop system.
        Acl,Bcl,Ccl,Dcl = self.closed_loop_system(sys_p, sys_c)

        # # Perform a step response.
        sys_cl = control.ss(Acl,Bcl,Ccl,Dcl)
        dtr = np.pi/180
        rtd = 180/np.pi
        step_mag = dtr
        sys_cl.A *= step_mag
        sys_cl.B *= step_mag
        timeoptimal_sec, yoptimal, x = control.step_response(sys_cl,return_x=True,transpose=True)

        # plt.figure()
        # plt.plot(timeoptimal_sec, yoptimal[:2])
        # plt.show()
        x = np.transpose(x)*np.pi/180

        # Reconstruct the control.
        # take first index of K matrix and append to end
        ind = np.roll(np.arange(nreg+n),-1)
        uoptimal = np.matmul(-Kw[0][ind],x)

        return uoptimal, yoptimal, timeoptimal_sec

#################################################################################
#    def perform_charting(self, sys_p, Aw, Bw, Creg, Qw, Rw, q, w):
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#################################################################################
    def closed_loop_system(self, plant_sys, controller_sys):
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Function computes the closed loop matrices for a particular plant and
    # controller.
    #
    #
    # Inputs:
    #   sys_p   structure, plant state space system
    #           Fields:
    #               .A     [nxn] matrix, State matrix
    #               .B     [nxm] matrix, Control matrix
    #               .C     [pxn] matrix, Output matrix
    #               .D     [pxm] matrix, Direct transmission matrix
    #   sys_c   structure, controller state space system
    #           Fields:
    #               .A     [ncxnc] matrix, Controller state matrix
    #               .B    [ncxmc1] matrix, Plant output matrix
    #               .B2    [ncxnr] matrix, Reference matrix
    #               .C     [mxnc] matrix, Controller output matrix
    #               .D    [mxmc1] matrix, Control-Plant output matrix
    #               .D2    [mxnr] matrix, Control-Reference matrix
    #
    #
    # Outputs:
    #   A   [nxn] matrix, Closed loop state matrix
    #   B   [nxm] matrix, Closed loop control matrix
    #   C   [pxn] matrix, Closed loop output matrix
    #   D   [pxm] matrix, Closed loop direct transmission matrix
    #
    #
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Plant Matrices
        Ap = plant_sys.A
        Bp = plant_sys.B
        Cp = plant_sys.C
        Dp = plant_sys.D

        # Controller Matrices
        Ac = controller_sys.A
        Bc1 = controller_sys.B
        Bc2 = controller_sys.B2
        Cc = controller_sys.C
        Dc1 = controller_sys.D
        Dc2 = controller_sys.D2

        # Intermediate matrices.
        Dc1Dp = np.matmul(Dc1,Dp)
        Z = np.subtract(np.eye(Dc1Dp.shape[0]),Dc1Dp)
        DpiZDc1 = Dp*np.linalg.inv(Z)*Dc1

        n = Ap.shape[0]
        m = Bp.shape[1]

        if (all(np.absolute(Ac[:]) < 1e-6) and all(abs(Bc1[:]) < 1e-6) and \
                all(abs(Bc2[:]) < 1e-6)):
            # No controller state defined.

            # Closed loop A matrix.
            A = Ap + Bp*np.linalg.inv(Z)*Dc1*Cp

            # Closed loop B matrix.
            B = Bp

            # Closed loop C matrix.
            C = Cp

            # Closed loop D matrix.
            D = Dp
        else:
            pass
            # Controller state defined.

            # Closed loop A matrix.

            A = np.zeros([n+m,n+m])
            A[:n,:n] = Ap + Bp*np.linalg.inv(Z)*Dc1*Cp
            A[:,n][:n] = (Bp*np.linalg.inv(Z)*Cc).flatten()

            A[n:,:n] = np.matmul(np.matmul(Bc1,(np.eye(DpiZDc1.shape[0])+DpiZDc1)),Cp)
            A[n,n] = Ac+np.matmul(np.matmul(np.matmul(Bc1,Dp),np.linalg.inv(Z)),Cc)

            # Closed loop B matrix.
            B = np.zeros([n+m,m])
            B[:n] = Bp*np.linalg.inv(Z)*Dc2
            B[n] = Bc2 + np.matmul(Bc1,Dp)*np.linalg.inv(Z)*Dc2

            # # Closed loop C matrix.
            C = np.zeros([n,n+m])
            C[:n,:n] = np.eye(DpiZDc1.shape[0])+DpiZDc1*Cp
            C[:n][:,n] = (Dp*np.linalg.inv(Z)*Cc).flatten()

            # # Closed loop D matrix.
            D = Dp*np.linalg.inv(Z)*Dc2

            return A,B,C,D


if __name__ == "__main__":
    from f16_model import F16_Model
    from gnc_lin_trim import Lin_Trim

    lintrim = Lin_Trim(F16_Model)
    # ic = [Vt_fps, AOA, altitude]
    input_conditions = np.array([250,-2*np.pi/180,5000])
    # initial guess for controls
    ux0 = np.array([0,-1.931,0,0,0.1485])
    # find trimmed x_dot, x, and u states
    # return is a dictionary with states labeled
    xd, x, u = lintrim.trim(input_conditions, ux0)[:3]
    # find linear matrices a and b
    a,b,c,d = lintrim.linearize_nonlinear_model(x, u, full_output=True)
    # Pull longitudinal and latitudinal(?) matrices from A
    Alon, Blon, Clon, Dlon = lintrim.model.build_long_state_matrices(a,b,c,d,full_output=True)
    Blon = Blon[:,1].reshape((Blon.shape[0],1))
    Dlon = Dlon[:,1].reshape((Dlon.shape[0],1))

    # M   = 10   # kg
    # B1  = 10   # Nm-sec
    # J1  = 15   # Nm-sec^2
    # J2  = 15   # Nm-sec^2
    # B2  = 100  # N-sec/m
    # K1  = 100  # Nm/rad
    # K2  = 130  # N/m
    # K3  = 160  # N/m
    # R   = 1    # m
    # A = np.array([[0, 1, 0, 0, 0, 0],
                  # [-K1/J1, 0, K1/J1, 0, 0 ,0],
                  # [0, 0, 0, 1, 0, 0],
                  # [K1/J2, 0, -(K1+R**2*K2)/J2, -B1/J2, R*K2/J2, 0],
                  # [0, 0, 0, 0, 0, 1],
                  # [0, 0, R*K2/M, 0, -(K2+K3)/M, -B2/M]])

    # B = np.array([0, 1/J1, 0, 0, 0, 0])
    # B.shape = (len(B),1)

    # n,m = B.shape
    # C = np.eye(n)

    # D = np.zeros([n,m])
    reg_array = np.array([2])
    controller = LQR_Controller()
    uoptimal, yoptimal, timeoptimal_sec = controller.lqr_control(Alon,Blon,Clon,Dlon,reg_array)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(timeoptimal_sec,yoptimal[:,3]*180/np.pi,'-b')
    ax.set_xlabel('Time(s)')
    ax.set_ylabel(r'$\theta$ degrees')





