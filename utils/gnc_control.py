#! /usr/bin/env
import numpy as np
import pprint
import control
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)

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

    def lqr_control(self,A,B,reg_array,Q,R,Q_end=None):
        n,m = B.shape
        C = np.eye(n)
        D = np.zeros([n,m])
        sys_ss = control.ss(A,B,C,D)

        # Optimal Contol: Tracker
        Creg = np.zeros(C.shape[1])
        Creg[reg_array] = 1
        Dreg = sys_ss.D[2,:]
        nreg = 1

        # Build the servomechanizm state space model of the plant by appending an
        # integral error state for each regulating channel.
        Aw = np.zeros([n+nreg, n+nreg])
        Aw[0:nreg][0,1:] = Creg

        Aw[nreg:,nreg:] = sys_ss.A

        Bw = np.zeros([n+nreg, m])
        Bw[0:nreg,:] = Dreg
        Bw[nreg:,:] = sys_ss.B

        # # Specify range of scale factor to apply to the integral error state.
        # q = np.logspace(1,6,300)

        # % Specify the LQR matrices.
        Qw = np.zeros([n+nreg, n+nreg])
        if Q_end:
            Qw[n+nreg-1,n+nreg-1] = Q_end

        # Rw = 1*np.eye(m)
        Rw = R*np.eye(m)

        # % Define frequency range for performing frequency analysis.
        # w = np.logspace(-3,3,600)


        # % Pass in the characteristics of the system.  Note, this process
        # % returns plots.  This process assumes Dreg == 0.

        # Perform LQR Design # Set the integral error state weights to the desired scale factor
        for j in range(nreg):
            Qw[j,j] = Q

        # Solve the LQR problem.
        Kw, S, E = control.lqr(Aw, Bw, Qw, Rw)
        # print(Kw)

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
        #Acl,Bcl,Ccl,Dcl = self.closed_loop_system(sys_p, sys_c)
        sys_cl = self.create_closed_loop_system(sys_ss,Kw,nreg)

        # Perform a step response.
        sys_cl = control.ss(Acl,Bcl,Ccl,Dcl)
        dtr = np.pi/180
        rtd = 180/np.pi
        step_mag = dtr
        timeoptimal_sec, yoptimal, x = control.step_response(sys_cl*0.1,T=np.arange(0,15,0.01),return_x=True,transpose=True)

        x = np.transpose(x)

        # Reconstruct the control.
        # take first index of K matrix and append to end
        ind = np.roll(np.arange(nreg+n),-1)
        uoptimal = np.matmul(-Kw[0][ind],x)

        return uoptimal, yoptimal, timeoptimal_sec

#################################################################################
    # def perform_charting(self, sys_ss, Aw, Bw, Creg, Qw, Rw, q, w):
        # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # PO = np.zeros([len(q),1])
        # tr = np.zeros_like(PO)
        # ts = np.zeros_like(PO)
        # wc = np.zeros_like(PO)
        # maxu = np.zeros_like(PO)
        # srmin = np.zeros_like(PO)
        # rdmin = np.zeros_like(PO)

        # # Number of channels to regulate
        # nreg = Creg.shape[0]

        # # Get the dimensions of the problem
        # n,m = sys_p.Bp.shape

        # # Define magnitude of step.
        # step_mag = 0.1

        # for i in range(len(q)):
            # # Set the integral error state weights to the current instance of the
            # # scale factor.
            # for j in range(nreg):
                # Qw[j,j] = q[i]

            # # Solve the LQR problem.
            # try:
                # Kw, S, E = control.lqr(Aw, Bw, Qw, Rw)
            # except:
                # PO[i] = np.nan
                # tr[i] = np.nan
                # ts[i] = np.nan
                # wc[i] = np.nan
                # rdmin[i] = np.nan
                # srmin[i] = np.nan
                # maxu[i] = np.nan
                # continue

            # sys_cl = self.create_closed_loop_system(sys_ss,Kw,nreg)

            # # Get the eigenvalues.
            # ee = eig(Acl)
            # if any(real(ee) > 0)
                # continue
            # end

            # # Perform frequency response from control (channel 1) to the ball
            # # position output (output 1).
            # freq_response = freq_analysis(sys_p, sys_c, w, 1, 3)

            # # Get the Loop Gain crossover frequency.
            # wc(i) = freq_response.input.LoopGainXover_rps

            # # Get the Return Difference, I + L, and its min magnitude.
            # rd = freq_response.input.ReturnDiff
            # rdmin(i) = min(abs(rd))

            # # Get the Stability Robustness, I + inv(L), and its min magnitdue.
            # sr = freq_response.input.StabRobust
            # srmin(i) = min(abs(sr))

            # [y, time_sec, x] = step(syscl*step_mag)

            # # Reconstruct the control and pick the biggest magnitude.
            # u = -Kw(:,[nreg+1:nreg+n,1:nreg])*x'
            # maxu(i) = max(abs(u))

            # if cond(Acl) > 50000
                # # Compute step information for the system.  Don't use stepinfo() as
                # # it does not deliver correct results for multi-outputs.
                # s = step_analysis(syscl, 1, 3, step_mag)

                # # Get the transient metrics.
                # PO(i) = s.PO
                # tr(i) = s.tr90_sec
                # ts(i) = s.ts98_sec
            # else
                # sout = stepinfo(syscl*step_mag)

                # # Get the transient metrics.
                # PO(i) = sout(3).Overshoot
                # tr(i) = sout(3).RiseTime
                # ts(i) = sout(3).SettlingTime
            # end
        # end


    def create_closed_loop_system(self,sys_ss,Kw,nreg):

        # % Form plant system.  Set output matrix for state feedback control.
        n = sys_ss.A.shape[0]
        sys_p = sys_setup()
        sys_p.A = sys_ss.A
        sys_p.B = sys_ss.B
        sys_p.C = np.eye(n)
        sys_p.D = np.zeros([n,1])

        # Set the state feedback gain and the integral error gain.
        kI = Kw[:,:nreg][0]
        Kv = Kw[:,nreg:][0]

        # Build controller object to track reference theta2.

        # Form closed loop system.
        Acl, Bcl, Ccl, Dcl = self.closed_loop_system(sys_p, sys_c)
        # Perform a step response.
        sys_cl = control.ss(Acl,Bcl,Ccl,Dcl)
        return sys_cl

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
    Alon, Blon = lintrim.model.build_long_state_matrices(a,b,c,d)
    Blon = Blon[:,1].reshape((Blon.shape[0],1))

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
    # reg_array = np.array([2])
    # controller = LQR_Controller()
    # Q = 30000
    # R = 0.1
    # Q_end = 10000
    # uoptimal, yoptimal, timeoptimal_sec = controller.lqr_control(Alon,Blon,reg_array,Q,R,Q_end)

    # fig = plt.figure()
    # ax1 = fig.add_subplot(2,1,1)
    # ax1.plot(timeoptimal_sec,yoptimal[:,2],'-b')
    # ax1.set_xlabel('Time(s)')
    # ax1.set_ylabel(r'$\theta$ radians')
    # ax2 = fig.add_subplot(2,1,2)
    # ax2.plot(timeoptimal_sec,uoptimal,'-b')
    # ax2.set_xlabel('Time(s)')
    # ax2.set_ylabel(r'$\delta_{Elevator}$ (degrees)')





