
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def analyze_state(self, a, b, c, d):
        xretain = np.array([0,1,2,3,4,5,6,7,8,12])
        A = a[xretain,:][:,xretain]
        B = b[xretain,:]
        c[1,7] = 180/np.pi
        c[2,1] = 180/np.pi
        C = c[:,xretain]
        D = d

        # Create the Longitudinal Modes Matrix, A, by only using the states
        # associated with Vt, AOA, theta and q.
        ind = self.long_ind
        Am = A[ind,:][:,ind]

        # Create the Longitudinal Input Matrix, B, by only using the states
        # associated with Vt, AOA, theta and q, and controls delT, dele.
        Bm = B[ind,:][:,tuple([0,1])]

        # Compute the eigenvalues and eigenvectors.
        w,v = np.linalg.eig(Am)

        # print(np.transpose(np.array(B[ind,:][:,1])))
        # print(C[1][ind])
        # print(D[2][2])
        B_in = B[ind,:][:,1]
        C_in = C[1][ind]
        D_in = D[2][2]

        B_in.shape = (4,1)
        C_in.shape = (1,4)

        num, den = signal.ss2tf(Am, B_in, C_in, D_in)
        G = signal.TransferFunction(num,den)
        #zpk = signal.tf2zpk(G.num,G.den)

        # short period approximation
        numa = np.array([1, 1.034])*-11.463
        numb = np.array([1, 2.516, 4.023])
        Ga = signal.TransferFunction(numa,numb)

        w,mag,phase = signal.bode(G)
        wa,maga,phasea = signal.bode(Ga)

        plt.figure(1)
        plt.subplot(211)
        plt.semilogx(w,mag,'b-',label='Actual')
        plt.semilogx(wa,maga,'b--',label='Aproximation')
        plt.legend()
        plt.subplot(212)
        plt.semilogx(w,phase,'b-',label='Actual')
        plt.semilogx(wa,phasea,'b--',label='Aproximation')
        plt.legend()
        plt.show()

        return Am,Bm,0,0

