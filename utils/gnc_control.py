#! /usr/bin/env python3

import numpy as np
import pprint
import control
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)

class control_object():
    def __init__(self):
        self.sys_ss = None
        self.sys_p = None
        self.sys_c = None
        self.sys_cl = None
        self.Aw = None
        self.Bw = None
        self.Qw = None
        self.Rw = None
        self.Creg = None
        self.Dreg = None
        self.Kw = None
        self.freq_fig = None
        self.dtr = np.pi/180
        self.rtd = 180/np.pi
        self.nreg = None
        self.n = None
        self.m = None
        
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

    def create_control_object(self,A,B,reg_array,Q,R,Q_end=None): 

        ctl = control_object()
        n,m = B.shape
        ctl.n = n
        ctl.m = m

        C = np.eye(n)
        D = np.zeros([n,m])
        ctl.sys_ss = control.ss(A,B,C,D)

        # Create Plant system
        ctl.sys_p = sys_setup()
        ctl.sys_p.A = ctl.sys_ss.A
        ctl.sys_p.B = ctl.sys_ss.B
        ctl.sys_p.C = np.eye(n)
        ctl.sys_p.D = np.zeros([n,1])

        # Optimal Contol: Tracker
        ctl.Creg = np.zeros(C.shape[1])
        ctl.Creg[reg_array] = 1
        ctl.Dreg = ctl.sys_ss.D[2,:]
        nreg = len(reg_array)
        ctl.nreg = nreg

        # Build the servomechanizm state space model of the plant by appending an
        # integral error state for each regulating channel.
        ctl.Aw = np.zeros([n+nreg, n+nreg])
        ctl.Aw[0:nreg][0,1:] = ctl.Creg

        ctl.Aw[nreg:,nreg:] = ctl.sys_ss.A

        ctl.Bw = np.zeros([n+nreg, m])
        ctl.Bw[0:nreg,:] = ctl.Dreg
        ctl.Bw[nreg:,:] = ctl.sys_ss.B

        # % Specify the LQR matrices.
        ctl.Qw = np.zeros([n+nreg, n+nreg])
        if Q_end:
            ctl.Qw[n+nreg-1,n+nreg-1] = Q_end

        ctl.Rw = R*np.eye(m)

        # Perform LQR Design # Set the integral error state weights to the desired scale factor
        for j in range(nreg):
            ctl.Qw[j,j] = Q

        # Solve the LQR problem.
        ctl.Kw, S, E = control.lqr(ctl.Aw, ctl.Bw, ctl.Qw, ctl.Rw)

        # Set the state feedback gain and the integral error gain.
        kI = ctl.Kw[:,:nreg][0]
        Kv = ctl.Kw[:,nreg:][0]

        # Build controller object to track reference theta
        ctl.sys_c = sys_setup()
        ctl.sys_c.A = np.zeros([nreg,nreg])
        ctl.sys_c.B= ctl.Creg
        ctl.sys_c.B2= -np.eye(nreg)
        ctl.sys_c.C = -kI
        ctl.sys_c.D= -Kv
        ctl.sys_c.D2= np.zeros([m,nreg])

        # Form closed loop system.
        Acl,Bcl,Ccl,Dcl = self.closed_loop_system(ctl.sys_p, ctl.sys_c)

        # Perform a step response.
        ctl.sys_cl = control.ss(Acl,Bcl,Ccl,Dcl)

        return ctl

#################################################################################
    def perform_charting(self, ctl_object):
        q = np.logspace(1,6,300)
        w = np.logspace(-3,3,600)
        PO = np.zeros([len(q),1])
        tr = np.zeros_like(PO)
        ts = np.zeros_like(PO)
        wc = np.zeros_like(PO)
        maxu = np.zeros_like(PO)
        srmin = np.zeros_like(PO)
        rdmin = np.zeros_like(PO)
        nreg = ctl_object.Creg.shape[0]

        Qw = np.zeros_like(ctl_object.Qw)

        for i in range(len(q)):
            for j in range(nreg):
                Qw[j,j] = q[i]

            try:
                # Solve the LQR problem.
                ctl.Kw, S, E = control.lqr(ctl.Aw, ctl.Bw, ctl.Qw, ctl.Rw)
            except:
                PO[i] = np.nan
                tr[i] = np.nan
                ts[i] = np.nan
                wc[i] = np.nan
                rdmin[i] = np.nan
                srmin[i] = np.nan
                maxu[i] = np.nan
                continue

            # Set the state feedback gain and the integral error gain.
            kI = ctl.Kw[:,:nreg][0]
            Kv = ctl.Kw[:,nreg:][0]

            # Build controller object to track reference theta
            sys_c = sys_setup()
            sys_c.A = np.zeros([nreg,nreg])
            sys_c.B= Creg
            sys_c.B2= -np.eye(nreg)
            sys_c.C = -kI
            sys_c.D= -Kv
            sys_c.D2= np.zeros([m,nreg])
                
            # Form closed loop system.
            Acl,Bcl,Ccl,Dcl = self.closed_loop_system(ctl.sys_p, ctl.sys_c)

            # Get the eigenvalues
            ee = np.eig(Acl)
            if any(np.real(ee) > 0):
                continue

#################################################################################
    def step_response(self,ctl,t=10,t_step=0.01):
        timeoptimal_sec, yoptimal, x = control.step_response(ctl.sys_cl*0.1,T=np.arange(0,3,t_step),return_x=True,transpose=True)

        x = np.transpose(x)*np.pi/180

        # Reconstruct the control.
        # take first index of K matrix and append to end
        ind = np.roll(np.arange(ctl.nreg+ctl.n),-1)
        uoptimal = np.matmul(-ctl.Kw[0][ind],x)

        return uoptimal, yoptimal, timeoptimal_sec

#################################################################################

    def freq_analysis(self,ctl_object,w):
        nf = len(w)
        ctl_object.freq_data.input.Lin = np.zeros([nf,1])
        ctl_object.freq_data.input.ReturnDiff = np.zeros([nf,1])
        ctl_object.freq_data.input.StabRobust = np.zeros([nf,1])
        ctl_object.freq_data.input.LoopGainXover_rps = 0
        ctl_object.freq_data.input.GMRD = np.array([0,0])
        ctl_object.freq_data.input.PMRD = np.array([0,0])
        ctl_object.freq_data.input.PMSR = np.array([0,0])
        ctl_object.freq_data.input.Gnoise = np.zeros([nf,1])

        # initialize output freq analysis substructure
        ctl_object.freq_data.output.S = np.zeros([nf,1])
        ctl_object.freq_data.output.T = np.zeros([nf,1])

        # Initialize system substructure
        ctl_object.freq_data.system.OLEvalues = np.array([])
        ctl_object.freq_data.system.CLEvalues = np.array([])





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

    reg_array = np.array([2])
    controller = LQR_Controller()
    Q = 400000
    R = 0.01
    Q_end = 300
    ctl = controller.create_control_object(Alon,Blon,reg_array,Q,R,Q_end)

    uoptimal, yoptimal, timeoptimal_sec = controller.step_response(ctl)

    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax1.plot(timeoptimal_sec,yoptimal[:,2]*180/np.pi,'-b')
    ax1.set_xlabel('Time(s)')
    ax1.set_ylabel(r'$\theta$ radians')
    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(timeoptimal_sec,uoptimal,'-b')
    ax2.set_xlabel('Time(s)')
    ax2.set_ylabel(r'$\delta_{Elevator}$ (degrees)')





