#! /usr/bin/python3

import numpy as np
from scipy import optimize
from utils.f16_model import F16_Model

class Lin_Trim(F16_Model):
    def __init__(self):
    #########################################################################################
        # Version history:
            # Nov 9, 2018    initial release, inspired by .m code written by Dr. Buckholtz
        # Data taken from:
            # Aircraft Control and Simulation, 2nd Edition
            # by Stevens and Lewis
    #########################################################################################
        F16_Model.__init__(self)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def trim(self, input_conditions, initial_guess):
    #########################################################################################
        # this function trims the input model based on the nonlinear model given in the cost
        # function for prescribed input_conditions

        # inputs:
            # input_conditions -- [3x1] vector, Input conditions of trim
                # [0]: Vt_fps
                # [1]: alt_ft
                # [2]: alpha_rad

            # initial_guess -- [3x1] vector, initial x states
                # [0]: beta_rad
                # [1]: del_el_deg
                # [2]: del_rud_deg
                # [3]: del_ail_deg
                # [4]: throttle

        # outputs:
    #########################################################################################
        #self.set_input_conditions(input_conditions)
        minimum = optimize.fmin(self.cost_function, initial_guess, args=(input_conditions,)\
                                                   xtol=1e-10, ftol=1e-10, maxfun=5e04, maxiter=1e05)

        if minimum[4]:
            print('Trim algorithm DID NOT converge')
        else:
            print('Trim algorithm converged.')
            return minimum[0]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def cost_function(self, ux0, input_conditions):
    #########################################################################################
        # this function computes the cost used in assessing a vehicle state as
        # being in trim

        # inputs:

            # ux0 -- [5x1] vector, initial conditions of the controls

                # ux0(0): beta (rad)
                # ux0(1): elevator deflection (deg)
                # ux0(2): aileron deflection (deg)
                # ux0(3): rudder deflection (deg)
                # ux0(4): throttle (0-1)

            # input_conditions -- [3x1] vector, initial conditions of the states
                # [0]: V_fps
                # [0]: alt_ft
                # [0]: alpha_rad

        # outputs:
            # cost        scalar, converged cost of the trimmed system
            # xd          n-element vector, derivative of state vector
            # x           n-element vector, trimmed state vector
            # u           n-element vector, trimmed control vector
    #########################################################################################

        # compute power
        power = self.f16_tgear(ux0[4])

        # build the state vector:
        x = np.array([[input_conditions[0]],           # vt
                      [input_conditions[1]],           # aoa
                      [ux0[0]],                              # aos
                      [0],
                      [0],
                      [0],
                      [0],
                      [0],
                      [0],
                      [0],
                      [0],
                      [input_conditions[2]],      # altitude
                      [power]])                          # power


        # build the control vector:
        u = [ux0[4],                # throttle
             ux0[1],                # del_el
             ux0[3],                # del_aileron
             ux0[2]]                # del_rudder

        u = np.array(u)
        u[:, np.newaxis]            # convert to column vector
        print(u)

        # compute the state derivatives
        xd = self.nonlinear_model(0, x, u)[0]

        # compute cost of vector.  cost is a weighted sum of the pertinent states required
        # we will drive this to zero
        weight = np.array([1,100,100,10,10,10])

        cost = np.matmul(weight,(np.array([[xd[0]],[xd[1]],[xd[2]],[xd[6]],[xd[7]],[xd[8]]])**2))

        return cost

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def linearize_nonlinear_model(self, x, u, u1):
    #########################################################################################
        # this function linearizes the nonlinear model defined by state input x
        # control input u

    #########################################################################################

        # Get the size of the state and controls
        n = len(x)
        m = len(u)
        tol = 1e-6
        time = 0.

        # Create the Jacobian for the A matrix by numerically computing the partial
        # derivatives of the state equations with respect to each state.
        dx = 0.1*x
        for i in range(n):
            if dx[i]==0.0:
                dx[i] = 0.1

        last = np.zeros(n)
        a = np.zeros([n,n])
        for j in range(n):
            xt = x
            for i in range(20):
                xt[j] = x[j] + dx[j]
                xd1 = self.nonlinear_model(time, xt, np.array([[u[0]],[u[1]],[u[2]],[u1]]))
                xt[j] = x[j] - dx[j]
                xd2 = self.nonlinear_model(time, xt, np.array([[u[0]],[u[1]],[u[2]],[u1]]))
                a[:,j]= np.divide(np.subtract(xd1,xd2),(2*dx[j]))
                if np.amax(np.divide(np.absolute(np.subtract(a[:,j],last)),np.absolute(a[:,j] + 1e-12))) < tol:
                    break
                dx[j] = 0.5*dx[j]
                last = a[:,j]
            # column=j
            iteration = i
            if iteration==19:
                print('not converged on A, column: '+str(j))

        # Create the Jacobian for the B matrix by numerically computing the partial
        # derivatives of the state equations with respect to each input.
        du = 0.1*u
        for i in range(m):
            if du[i]==0.0:
                du[i]=0.1

        last = np.zeros(n)
        b = np.zeros([n,m])
        for j in range(m):
            usave = u
            for i in range(10):
                u[j] = usave[j] + du[j]
                xd1 = self.nonlinear_model(time, x, np.array([[u[0]],[u[1]],[u[2]],[u1]]))
                u[j] = usave[j] - du[j]
                xd2 = self.nonlinear_model(time, x, np.array([[u[0]],[u[1]],[u[2]],[u1]]))
                b[:,j]= np.divide(np.subtact(xd1,xd2),(2*du[j]))
                if np.amax(np.divide(np.absolute(np.subtract(b[:,j],last)), np.absolute(b[:,j] + 1e-12))) < tol:
                    break
                du[j] = 0.5*du[j]
                last = b[:,j]
            # column=j
            iteration = i
            if iteration==10:
                print('not converged on B, column: '+str(j))

        return a, b


if __name__ == "__main__":
    lintrim = Lin_Trim()
    lintrim.trim(np.array([500,2*np.pi/180,1]), np.array([0,-1.931,0,0,0.1485]))
