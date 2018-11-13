#! /usr/bin/python3

import numpy as np
from scipy import optimize
from f16_model import F16_Model

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
            # xd -- Dictionary containing derivative states and values
            # x -- Dictionary containing states and values
            # u -- Dictionary containing control and values
    #########################################################################################
        #self.set_input_conditions(input_conditions)
        output = optimize.fmin(self.cost_function, initial_guess, args=(input_conditions,),\
                                                   xtol=1e-10, ftol=1e-10, maxfun=5e04, maxiter=1e05, \
                                                   full_output = True)

        if output[4]:
            print('Trim algorithm DID NOT converge')
        else:
            print('Trim algorithm converged.')
            x,u = self.build_x_u(output[0],input_conditions)
            xd = self.nonlinear_model(0, x, u)
            return xd, x, u

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

        x, u = self.build_x_u(ux0, input_conditions)

        # compute the state derivatives
        xd = self.nonlinear_model(0, x, u)[0]

        # compute cost of vector.  cost is a weighted sum of the pertinent states required
        # we will drive this to zero
        weight = np.array([1,100,100,10,10,10])

        xd_array = np.array([xd['Airspeed_dot'],
                             xd['alpha_dot'],
                             xd['beta_dot'],
                             xd['Roll_Rate_dot'],
                             xd['Pitch_Rate_dot'],
                             xd['Yaw_Rate_dot']])

        cost = np.matmul(weight,(xd_array**2))
        return cost


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def linearize_nonlinear_model(self, x_dict, u_dict):
    #########################################################################################
        # this function linearizes the nonlinear model defined by state input x
        # control input u
    #########################################################################################

        # Convert Dictionaries to column vectors
        x = self._dict_to_colvector(x_dict)
        u = self._dict_to_colvector(u_dict)
        print(x)

        # Get the size of the state and controls
        n = len(x)
        xt = np.zeros([n,1])
        m = len(u)
        ut = np.zeros([m,1])
        tol = 1e-6
        time = 0.

        usave = np.zeros([m,1])
        usave[:] = u[:]

        # Set number of outputs
        p = 3

        # Create the Jacobian for the A matrix by numerically computing the partial
        # derivatives of the state equations with respect to each state.
        dx = 0.1*x
        for i in range(n):
            if dx[i][0] == 0.0:
                dx[i][0] = 0.1

        last = np.zeros([n,1])
        a = np.zeros([n,n])
        for j in range(n):
            xt[:] = x[:]
            for i in range(20):
                xt[j][0] = x[j][0] + dx[j][0]
                xd1 = self.nonlinear_model(time, xt, u)[0]
                xd1 = self._dict_to_colvector(xd1)
                xt[j][0] = x[j][0] - dx[j][0]
                xd2 = self.nonlinear_model(time, xt, u)[0]
                xd2 = self._dict_to_colvector(xd2)
                a[:,j]= np.divide(np.subtract(xd1,xd2),(2*dx[j][0]))[:,0]
                if np.amax(np.divide(np.absolute(np.subtract(a[:,j],last)),np.absolute(a[:,j] + 1e-12))) < tol:
                    break
                dx[j][0] = 0.5*dx[j][0]
                last = a[:,j]
            # column=j
            iteration = i
            if iteration==19:
                print('not converged on A, column: '+str(j))

        # Create the Jacobian for the B matrix by numerically computing the partial
        # derivatives of the state equations with respect to each input.
        du = 0.1*u
        for i in range(m):
            if du[i][0]==0.0:
                du[i][0]=0.1

        last = np.zeros(n)
        b = np.zeros([n,m])
        for j in range(m):
            ut[:] = u[:]
            for i in range(10):
                u[j][0] = ut[j][0] + du[j][0]
                xd1 = self.nonlinear_model(time, x, u)[0]
                xd1 = self._dict_to_colvector(xd1)
                u[j][0] = ut[j][0] - du[j][0]
                xd2 = self.nonlinear_model(time, x, u)[0]
                xd2 = self._dict_to_colvector(xd2)
                b[:,j]= np.divide(np.subtract(xd1,xd2),(2*du[j][0]))[:,0]
                if np.amax(np.divide(np.absolute(np.subtract(b[:,j],last)), np.absolute(b[:,j] + 1e-12))) < tol:
                    break
                du[j][0] = 0.5*du[j][0]
                last = b[:,j]
            # column=j
            iteration = i
            if iteration==10:
                print('not converged on B, column: '+str(j))

        # Create the Jacobian for the C matrix by numerically conputing the partial
        # derivative of the output with respect to each state.
        dx = 0.1*x
        for i in range(n):
            if dx[i][0] == 0.0:
                dx[i][0] = 0.1

        last = np.zeros([p,1])
        c = np.zeros([p,n])
        for j in range(n):
            xt[:] = x[:]
            for i in range(10):
                xt[j][0] = x[j][0] + dx[j][0]
                xd1,az1 = self.nonlinear_model(time, xt, u)[:2]
                xd1 = self._dict_to_colvector(xd1)
                xt[j][0] = x[j][0] - dx[j][0]
                xd2,az2 = self.nonlinear_model(time, xt, u)[:2]
                xd2 = self._dict_to_colvector(xd2)
                c[0][j] = (az1-az2)/(2*dx[j][0])
                if np.amax(np.divide(np.absolute(np.subtract(c[:,j],last)), np.absolute(c[:,j] + 1e-12))) < tol:
                    break
                dx[j][0] = 0.5*dx[j][0]
                last = c[:,j]
            # column=j
            iteration = i
            if iteration==10:
                print('Not Converged on C, column: ' + str(j))

        # Compute the Jacobian of the D matrix by numerically computing the partial
        # derivative of the output with respect to each input.
        du = 0.1*u
        for i in range(m):
           if du[i][0] == 0.0:
              du[i][0] = 0.1

        last = np.zeros([p,1])
        d = np.zeros([p,m])
        u[:] = usave[:]
        for j in range(m):
            ut[:] = u[:]
            for i in range(10):
                u[j][0] = ut[j][0] + du[j][0]
                xd1,az1,ay1 = self.nonlinear_model(time,x,u)
                xd1 = self._dict_to_colvector(xd1)
                u[j][0] = ut[j][0] - du[j][0]
                xd2,az2,ay2 = self.nonlinear_model(time,x,u)
                xd2 = self._dict_to_colvector(xd2)
                d[0][j] = (az1-az2)/(2*du[j][0])
                if np.amax(np.divide(np.absolute(np.subtract(d[:,j],last)), np.absolute(d[:,j] + 1e-12))) < tol:
                    break
                du[j][0] = 0.5*du[j][0]
                last = d[:,j]
           #column=j
            iteration = i
            if iteration==10:
                print('Not Converged on D, column: ' + str(j))

        return a,b,c,d

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def linearize_trim_state(self,x,u):
        a,b,c,d = self.linearize_nonlinear_model(x, u)
        return a,b,c,d

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def analyze_linear_state(self,a,b,c,d):
        self.analyze_state(a,b,c,d)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def _dict_to_colvector(self, dictionary):
        l = [[x] for x in dictionary.values()]
        return np.array(l)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == "__main__":
    lintrim = Lin_Trim()
    xd, x, u = lintrim.trim(np.array([500,2*np.pi/180,1]), np.array([0,-1.931,0,0,0.1485]))
    a,b,c,d = lintrim.linearize_trim_state(x, u)
    lintrim.analyze_linear_state(a,b,c,d)
