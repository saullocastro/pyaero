"""
-----------------------------------------------------------------------
Computation of Matrices of Influence Coefficients

Author: Higor Luis Silva
-----------------------------------------------------------------------
"""
import numpy as np


class DLMParameters(object):
    def __init__(self):
        # Flight parameters
        self.M = 0.25
        # Mesh parameters
        self.sym = 0 # 1: symmetric ; -1: anti-symmetric; 0: no symmetric
        self.xboxes = 10 # Number of panels in x direction
        self.yboxes = 30 # Number of panels in y direction
        self.apex = np.array([0, 0, 0], dtype=np.float64)     # Global coordinates (x,y,z) - Apex
        self.inbLE = np.array([0, 0, 0], dtype=np.float64)  # Local coordinates (x,y,z) -  Inboard Leading Edge
        self.inbTE = np.array([1.0, 0, 0], dtype=np.float64)    # Local coordinats (x,y,z) - Inboard Trailing Edge
        self.outbTE = np.array([1.3, 3.5, 0], dtype=np.float64)  # Local coordinates (x,y,z) - Outboard Trailing Edge
        self.outbLE = np.array([0.5, 3.5, 0], dtype=np.float64) # Local coordinates (x,y,z) - Outboard Leading Edge
        # Kernel Interpolation Parameters
        self.eps = 1e-10
        self.bigp = 1e30
        self.bigm = -1e30
        self.degree = 4      # Interpolation order: 2 or 4


def calcD(wing, box, totalboxes, k, paramDLM=None):
    if paramDLM is None:
        paramDLM = DLMParameters()
    M = paramDLM.M
    eps = paramDLM.eps
    bigp = paramDLM.bigp
    bigm = paramDLM.bigm
    degree = paramDLM.degree

    b = wing.b
    tipo = wing.tipo

    n = 11
    c = 0.372
    a = an = np.array([.24186198,-2.7918027,24.991079,-111.59196,271.43549,-305.75288,\
    -41.18363,545.98537,-644.78155,328.72755,-64.279511],dtype=np.complex128)

    D = np.zeros((totalboxes, totalboxes), dtype=np.complex128)
    Ds = np.zeros((totalboxes, totalboxes), dtype=np.complex128)

    for rb in range(totalboxes):
        print("rb = ",rb)
        for sb in range(totalboxes):
            gamas = box.dihedral[sb]
            if tipo == 2:
                e = np.abs(box.z14o[sb]-box.z14i[sb])/2
            else:
                e = np.abs((box.y14o[sb]-box.y14i[sb])/(2*np.cos(gamas)))

            x0i = box.x34c[rb]-box.x14i[sb]
            y0i = box.y34c[rb]-box.y14i[sb]
            z0i = box.z34c[rb]-box.z14i[sb]

            x0o = box.x34c[rb]-box.x14o[sb]
            y0o = box.y34c[rb]-box.y14o[sb]
            z0o = box.z34c[rb]-box.z14o[sb]

            x0 = box.x34c[rb]-box.x14c[sb]
            y0 = box.y34c[rb]-box.y14c[sb]
            z0 = box.z34c[rb]-box.z14c[sb]

            x = box.x34c[rb]
            y = box.y34c[rb]
            z = box.z34c[rb]
            xiC = box.x14c[sb]
            etaC = box.y14c[sb]
            zetaC = box.z14c[sb]
            gamar = box.dihedral[rb]
            gamas = box.dihedral[sb]
            xbar = x-xiC
            ybar = (y-etaC)*np.cos(gamas)+(z-zetaC)*np.sin(gamas)
            zbar = (z-zetaC)*np.cos(gamas)-(y-etaC)*np.sin(gamas)
            lambdas = box.lambdas[sb]
            csi = box.x14c[sb]
            eta = box.y14c[sb]
            zeta = box.z14c[sb]
            xibar = csi-xiC
            etabar = (eta-etaC)*np.cos(gamas)+(zeta-zetaC)*np.sin(gamas)
            zetabar = (zeta-zetaC)*np.cos(gamas)-(eta-etaC)*np.sin(gamas)
            deltaxs = box.chord[sb]

            Ds[rb,sb] = Drss(x-xiC,y-etaC,z-zetaC,gamar,gamas,lambdas,e,deltaxs,paramDLM)

    for rb in range(totalboxes):
        print("rb = ",rb)
        for sb in range(totalboxes):
            gamas = box.dihedral[sb]
            if tipo == 2:
                e = np.abs(box.z14o[sb]-box.z14i[sb])/2
            else:
                e = np.abs((box.y14o[sb]-box.y14i[sb])/(2*np.cos(gamas)))

            x0i = box.x34c[rb]-box.x14i[sb]
            y0i = box.y34c[rb]-box.y14i[sb]
            z0i = box.z34c[rb]-box.z14i[sb]

            x0o = box.x34c[rb]-box.x14o[sb]
            y0o = box.y34c[rb]-box.y14o[sb]
            z0o = box.z34c[rb]-box.z14o[sb]

            x0 = box.x34c[rb]-box.x14c[sb]
            y0 = box.y34c[rb]-box.y14c[sb]
            z0 = box.z34c[rb]-box.z14c[sb]
            x = box.x34c[rb]
            y = box.y34c[rb]
            z = box.z34c[rb]
            xiC = box.x14c[sb]
            etaC = box.y14c[sb]
            zetaC = box.z14c[sb]
            gamar = box.dihedral[rb]
            gamas = box.dihedral[sb]
            xbar = x-xiC
            ybar = (y-etaC)*np.cos(gamas)+(z-zetaC)*np.sin(gamas)
            zbar = (z-zetaC)*np.cos(gamas)-(y-etaC)*np.sin(gamas)
            lambdas = box.lambdas[sb]
            csi = box.x14c[sb]
            eta = box.y14c[sb]
            zeta = box.z14c[sb]
            xibar = csi-xiC
            etabar = (eta-etaC)*np.cos(gamas)+(zeta-zetaC)*np.sin(gamas)
            zetabar = (zeta-zetaC)*np.cos(gamas)-(eta-etaC)*np.sin(gamas)
            deltaxs = box.chord[sb]

            zbar = zbar-zetabar

            if degree == 2:
                vec_x0 = np.array([x0i,x0,x0o],dtype=np.complex128)
                vec_y0 = np.array([y0i,y0,y0o],dtype=np.complex128)
                vec_z0 = np.array([z0i,z0,z0o],dtype=np.complex128)
                vec_KP1 = np.zeros([0,0,0])

                for station in range(3):

                    x0 = vec_x0[station]
                    y0 = vec_y0[station]
                    z0 = vec_z0[station]

                    beta2 = 1-M*M
                    r = np.sqrt(y0*y0+z0*z0)

                    T1 = np.cos(gamar-gamas)

                    if x0 == 0 and r == 0:
                        K1 = 0
                        K1s = 0
                    else:
                        R = np.sqrt(x0*x0+beta2*r*r)
                        k1 = k*r/b

                        if np.abs(r) < eps:
                            if x0.real > 0:
                                u1 = bigm
                            else:
                                u1 = bigp
                        else:
                            u1 = (M*R-x0)/(beta2*r)

                        K1s = 1+x0/R

                        if u1.real >= 0:
                            if u1.real >= bigp:
                                II1 = 0
                            else:
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*np.exp(-(di+1)*c*u1)/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                factor = 1-u1/np.sqrt(1+u1*u1)
                                i1_u1_k1 = np.exp(-1j*k1*u1)*(factor-1j*k1*II0)
                                II1 = i1_u1_k1

                        else:
                            if u1.real <= bigm:
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*1/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                i1_0_k1 = 1-1j*k1*II0
                                II1 = 2*i1_0_k1
                            else:
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*1/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                i1_0_k1 = 1-1j*k1*II0
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*np.exp((di+1)*c*u1)/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                factor = 1+u1/np.sqrt(1+u1*u1)
                                i1_mu1_k1 = np.exp(1j*k1*u1)*(factor-1j*k1*II0)
                                II1 = i1_mu1_k1
                                II1 = (2*np.real(i1_0_k1)-np.real(II1))+(np.imag(II1))*1j



                        if k == 0:
                            K1 = K1s
                        else:
                            K1 = II1+(M*r/R)*(np.exp(-1j*k1*u1)/np.sqrt(1+u1*u1))



                    vec_KP1[station] = (K1*np.exp(-1j*k*x0/b)-K1s)*T1


                KL = vec_KP1[0]
                KC = vec_KP1[1]
                KR = vec_KP1[2]

                A1 = 1/2*(-2*KC+KL+KR)/(e*e)
                B1 = -1/2*(KL-KR)/e
                C1 = KC

                if np.abs(zbar)<eps:
                    D_ = (C1-A1*e*e+B1*ybar+2*A1*e*ybar+(-e+ybar)*(B1+2*A1*ybar)*np.log(-e+ybar))/(-e+ybar)-(C1-A1*e*e+B1*ybar-2*A1*e*ybar+(e+ybar)*(B1+2*A1*ybar)*np.log(e+ybar))/(e+ybar)
                    D_ = (deltaxs/(8*np.pi))*D_
                else:
                    D_ = (4*A1*e*zbar+2*(C1+ybar*(B1+A1*ybar)-A1*zbar*zbar)*np.arctan((e-ybar)/zbar)+2*(C1+ybar*(B1+A1*ybar)-A1*zbar*zbar)*np.arctan((e+ybar)/zbar)+(B1+2*A1*ybar)*zbar*(np.log((e-ybar)*(e-ybar)+zbar*zbar)-np.log((e+ybar)*(e+ybar)+zbar*zbar)))/(2*zbar)
                    D_ = (deltaxs/(8*np.pi))*D_

            elif degree == 4:
                x0LM = (x0i+x0)/2
                y0LM = (y0i+y0)/2
                z0LM = (z0i+z0)/2
                x0RM = (x0o+x0)/2
                y0RM = (y0o+y0)/2
                z0RM = (z0o+z0)/2
                vec_x0 = np.array([x0i,x0LM,x0,x0RM,x0o], dtype=np.complex128)
                vec_y0 = np.array([y0i,y0LM,y0,y0RM,y0o], dtype=np.complex128)
                vec_z0 = np.array([z0i,z0LM,z0,z0RM,z0o], dtype=np.complex128)
                vec_KP1 = np.array([0,0,0,0,0], dtype=np.complex128)

                for station in range(5):

                    x0 = vec_x0[station]
                    y0 = vec_y0[station]
                    z0 = vec_z0[station]

                    beta2 = 1-M*M
                    r = np.sqrt(y0*y0+z0*z0)

                    T1 = np.cos(gamar-gamas)

                    if x0 == 0 and r == 0:
                        K1 = 0
                        K1s = 0
                    else:
                        R = np.sqrt(x0*x0+beta2*r*r)
                        k1 = k*r/b

                        if np.abs(r) < eps:
                            if x0.real > 0:
                                u1 = bigm
                            else:
                                u1 = bigp

                        else:
                            u1 = (M*R-x0)/(beta2*r)


                        K1s = 1+x0/R

                        if u1.real >= 0:
                            if u1.real >= bigp:
                                II1 = 0
                            else:
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*np.exp(-(di+1)*c*u1)/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                factor = 1-u1/np.sqrt(1+u1*u1)
                                i1_u1_k1 = np.exp(-1j*k1*u1)*(factor-1j*k1*II0)
                                II1 = i1_u1_k1

                        else:
                            if u1.real <= bigm:
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*1/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                i1_0_k1 = 1-1j*k1*II0
                                II1 = 2*i1_0_k1
                            else:
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*1/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                i1_0_k1 = 1-1j*k1*II0
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*np.exp((di+1)*c*u1)/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                factor = 1+u1/np.sqrt(1+u1*u1)
                                i1_mu1_k1 = np.exp(1j*k1*u1)*(factor-1j*k1*II0)
                                II1 = i1_mu1_k1
                                II1 = (2*np.real(i1_0_k1)-np.real(II1))+(np.imag(II1))*1j



                        if k == 0:
                            K1 = K1s
                        else:
                            K1 = II1+(M*r/R)*(np.exp(-1j*k1*u1)/np.sqrt(1+u1*u1))



                    vec_KP1[station] = (K1*np.exp(-1j*k*x0/b)-K1s)*T1


                KL = vec_KP1[0]
                KML = vec_KP1[1]
                KC = vec_KP1[2]
                KMR = vec_KP1[3]
                KR = vec_KP1[4]

                A1 = 2/3*(6*KC-4*KMR-4*KML+KL+KR)/e**4
                B1 = -2/3*(2*KMR-2*KML+KL-KR)/e**3
                C1 = -1/6/e**2*(-16*KMR-16*KML+KL+KR+30*KC)
                D1 = 1/6*(8*KMR-8*KML+KL-KR)/e
                E1 = KC

                if np.abs(zbar) < eps:
                    D_ = (2*e*(3*C1*(e*e-2*ybar*ybar)-3*(E1+D1*ybar-2*B1*e*e*ybar+3*B1*ybar**3)+A1*(e**4+8*e*e*ybar*ybar-12*ybar**4))+3*(e-ybar)*(e+ybar)*(D1+ybar*(2*C1+ybar*(3*B1+4*A1*ybar)))*(np.log(-e+ybar)-np.log(e+ybar)))/(3*(e-ybar)*(e+ybar))
                    D_ = (deltaxs/(8*np.pi))*D_
                else:
                    D_ = (6*(E1+ybar*(D1+ybar*(C1+ybar*(B1+A1*ybar)))-(C1+3*ybar*(B1+2*A1*ybar))*zbar*zbar+A1*zbar**4)*np.arctan((e-ybar)/zbar)+6*(E1+ybar*(D1+ybar*(C1+ybar*(B1+A1*ybar)))-(C1+3*ybar*(B1+2*A1*ybar))*zbar*zbar+A1*zbar**4)*np.arctan((e+ybar)/zbar)+zbar*(4*e*(3*C1+6*B1*ybar+A1*(e*e+9*ybar*ybar-3*zbar*zbar))+3*(D1+2*C1*ybar+3*B1*ybar*ybar+4*A1*ybar**3-(B1+4*A1*ybar)*zbar*zbar)*(np.log((e-ybar)**2+zbar*zbar)-np.log((e+ybar)**2+zbar*zbar))))/(6*zbar)
                    D_ = (deltaxs/(8*np.pi))*D_

            D1 = D_

            gamas = box.dihedral[sb]
            if tipo == 2:
                e = np.abs(box.z14o[sb]-box.z14i[sb])/2
            else:
                e = np.abs((box.y14o[sb]-box.y14i[sb])/(2*np.cos(gamas)))

            x0i = box.x34c[rb]-box.x14i[sb]
            y0i = box.y34c[rb]-box.y14i[sb]
            z0i = box.z34c[rb]-box.z14i[sb]

            x0o = box.x34c[rb]-box.x14o[sb]
            y0o = box.y34c[rb]-box.y14o[sb]
            z0o = box.z34c[rb]-box.z14o[sb]

            x0 = box.x34c[rb]-box.x14c[sb]
            y0 = box.y34c[rb]-box.y14c[sb]
            z0 = box.z34c[rb]-box.z14c[sb]
            x = box.x34c[rb]
            y = box.y34c[rb]
            z = box.z34c[rb]
            xiC = box.x14c[sb]
            etaC = box.y14c[sb]
            zetaC = box.z14c[sb]
            gamar = box.dihedral[rb]
            gamas = box.dihedral[sb]
            xbar = x-xiC
            ybar = (y-etaC)*np.cos(gamas)+(z-zetaC)*np.sin(gamas)
            zbar = (z-zetaC)*np.cos(gamas)-(y-etaC)*np.sin(gamas)
            lambdas = box.lambdas[sb]
            csi = box.x14c[sb]
            eta = box.y14c[sb]
            zeta = box.z14c[sb]
            xibar = csi-xiC
            etabar = (eta-etaC)*np.cos(gamas)+(zeta-zetaC)*np.sin(gamas)
            zetabar = (zeta-zetaC)*np.cos(gamas)-(eta-etaC)*np.sin(gamas)
            deltaxs = box.chord[sb]

            zbar = zbar-zetabar

            if degree == 2:
                vec_x0 = np.array([x0i,x0,x0o], dtype=np.complex128)
                vec_y0 = np.array([y0i,y0,y0o], dtype=np.complex128)
                vec_z0 = np.array([z0i,z0,z0o], dtype=np.complex128)
                vec_KP2 = np.array([0,0,0], dtype=np.complex128)

                for station in range(3):

                    x0 = vec_x0[station]
                    y0 = vec_y0[station]
                    z0 = vec_z0[station]

                    beta2 = 1-M*M
                    r = np.sqrt(y0*y0+z0*z0)

                    T2est = (z0*np.cos(gamar)-y0*np.sin(gamar))*(z0*np.cos(gamas)-y0*np.sin(gamas))

                    if x0==0 and r==0:
                        K2 = 0
                        K2s = 0
                    else:
                        R = np.sqrt(x0*x0+beta2*r*r)
                        k1 = k*r/b

                        if np.abs(r) < eps:
                            if x0.real > 0:
                                u1 = bigm
                            else:
                                u1 = bigp

                        else:
                            u1 = (M*R-x0)/(beta2*r)


                        K2s = -2-(x0/R)*(2+beta2*r*r/(R*R))


                        if u1.real >= 0:
                            if u1.real >= bigp:
                                II2 = 0
                            else:
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*np.exp(-(di+1)*c*u1)/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                JJ0 = 0
                                for di in range(n):
                                    factor1 = a[di]*np.exp(-(di+1)*c*u1)/((((di+1)*(di+1)*c*c)+(k1*k1))*(((di+1)*(di+1)*c*c)+(k1*k1)))
                                    factor2 = (di+1)*(di+1)*c*c-k1*k1+(di+1)*c*u1*(((di+1)*(di+1)*c*c)+(k1*k1))-1j*k1*(2*(di+1)*c+u1*(((di+1)*(di+1)*c*c)+(k1*k1)))
                                    JJ0 = JJ0+factor1*factor2

                                factor = 1-u1/np.sqrt(1+u1*u1)
                                i2_u1_k1 = np.exp(-1j*k1*u1)*((2+1j*k1*u1)*factor-u1/((1+u1*u1)*np.sqrt(1+u1*u1))-1j*k1*II0+k1*k1*JJ0)
                                i2_u1_k1 = i2_u1_k1/3
                                II2 = i2_u1_k1

                        else:
                            if u1.real <= bigm:
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*1/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                JJ0 = 0
                                for di in range(n):
                                    factor1 = a[di]*1/((((di+1)*(di+1)*c*c)+(k1*k1))*(((di+1)*(di+1)*c*c)+(k1*k1)))
                                    factor2 = (di+1)*(di+1)*c*c-k1*k1-1j*k1*(2*(di+1)*c)
                                    JJ0 = JJ0+factor1*factor2

                                i2_0_k1 = (2-1j*k1*II0+k1*k1*JJ0)
                                i2_0_k1 = i2_0_k1/3
                                II2 = 2*i2_0_k1
                            else:
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*1/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                JJ0 = 0
                                for di in range(n):
                                    factor1 = a[di]*1/((((di+1)*(di+1)*c*c)+(k1*k1))*(((di+1)*(di+1)*c*c)+(k1*k1)))
                                    factor2 = (di+1)*(di+1)*c*c-k1*k1-1j*k1*(2*(di+1)*c)
                                    JJ0 = JJ0+factor1*factor2

                                i2_0_k1 = (2-1j*k1*II0+k1*k1*JJ0)
                                i2_0_k1 = i2_0_k1/3
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*np.exp((di+1)*c*u1)/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                JJ0 = 0
                                for di in range(n):
                                    factor1 = a[di]*np.exp((di+1)*c*u1)/((((di+1)*(di+1)*c*c)+(k1*k1))*(((di+1)*(di+1)*c*c)+(k1*k1)))
                                    factor2 = (di+1)*(di+1)*c*c-k1*k1-(di+1)*c*u1*(((di+1)*(di+1)*c*c)+(k1*k1))-1j*k1*(2*(di+1)*c-u1*(((di+1)*(di+1)*c*c)+(k1*k1)))
                                    JJ0 = JJ0+factor1*factor2

                                factor = 1+u1/np.sqrt(1+u1*u1)
                                i2_mu1_k1 = np.exp(1j*k1*u1)*((2-1j*k1*u1)*factor+u1/((1+u1*u1)*np.sqrt(1+u1*u1))-1j*k1*II0+k1*k1*JJ0)
                                i2_mu1_k1 = i2_mu1_k1/3
                                II2 = (2*np.real(i2_0_k1)-np.real(i2_mu1_k1))+(np.imag(i2_mu1_k1))*1j



                        if k == 0:
                            K2 = K2s
                        else:
                            K2 = -3*II2-((1j*k1*M*M*r*r)/(R*R)*(np.exp(-1j*k1*u1)/np.sqrt(1+u1*u1)))-M*r/R*((1+u1*u1)*(beta2*r*r/(R*R))+2+M*r*u1/R)*(np.exp(-1j*k1*u1)/((1+u1*u1)*np.sqrt(1+u1*u1)))



                    vec_KP2[station] = (K2*np.exp(-1j*k*x0/b)-K2s)*T2est


                KL = vec_KP2[0]
                KC = vec_KP2[1]
                KR = vec_KP2[2]


                A2 = 1/2*(-2*KC+KL+KR)/(e*e)
                B2 = -1/2*(KL-KR)/e
                C2 = KC

                if np.ans(zbar) < eps:
                    D_ = (-2*(C2*(e**3+3*e*ybar*ybar)+e**3*(4*B2*ybar+A2*(3*e*e+ybar*ybar))))/(3*(e*e-ybar*ybar)**3)
                    D_ = (deltaxs/(8*np.pi))*D_
                else:
                    D_ = ((2*e*zbar*((e-ybar)*(e+ybar)*(C2+ybar*(B2+A2*ybar))+(C2-A2*e*e-B2*ybar-2*A2*ybar*ybar)*zbar*zbar-A2*zbar**4))/(e**4-2*e*e*(ybar-zbar)*(ybar+zbar)+(ybar*ybar+zbar*zbar)**2)+(C2+B2*ybar+A2*(ybar*ybar+zbar*zbar))*(np.arctan((e-ybar)/zbar)+np.arctan((e+ybar)/zbar)))/(2*zbar**3)
                    D_ = (deltaxs/(8*np.pi))*D_

            elif degree == 4:
                x0LM = (x0i+x0)/2
                y0LM = (y0i+y0)/2
                z0LM = (z0i+z0)/2
                x0RM = (x0o+x0)/2
                y0RM = (y0o+y0)/2
                z0RM = (z0o+z0)/2
                vec_x0 = np.array([x0i,x0LM,x0,x0RM,x0o], dtype=np.complex128)
                vec_y0 = np.array([y0i,y0LM,y0,y0RM,y0o], dtype=np.complex128)
                vec_z0 = np.array([z0i,z0LM,z0,z0RM,z0o], dtype=np.complex128)
                vec_KP2 = np.array([0,0,0,0,0], dtype=np.complex128)

                for station in range(5):

                    x0 = vec_x0[station]
                    y0 = vec_y0[station]
                    z0 = vec_z0[station]

                    beta2 = 1-M*M
                    r = np.sqrt(y0*y0+z0*z0)

                    T2est = (z0*np.cos(gamar)-y0*np.sin(gamar))*(z0*np.cos(gamas)-y0*np.sin(gamas))

                    if x0 == 0 and r == 0:
                        K2 = 0
                        K2s = 0
                    else:
                        R = np.sqrt(x0*x0+beta2*r*r)
                        k1 = k*r/b

                        if np.abs(r) < eps:
                            if x0.real > 0:
                                u1 = bigm
                            else:
                                u1 = bigp

                        else:
                            u1 = (M*R-x0)/(beta2*r)


                        K2s = -2-(x0/R)*(2+beta2*r*r/(R*R))



                        if u1.real >= 0:
                            if u1.real >= bigp:
                                II2 = 0
                            else:
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*np.exp(-(di+1)*c*u1)/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                JJ0 = 0
                                for di in range(n):
                                    factor1 = a[di]*np.exp(-(di+1)*c*u1)/((((di+1)*(di+1)*c*c)+(k1*k1))*(((di+1)*(di+1)*c*c)+(k1*k1)))
                                    factor2 = (di+1)*(di+1)*c*c-k1*k1+(di+1)*c*u1*(((di+1)*(di+1)*c*c)+(k1*k1))-1j*k1*(2*(di+1)*c+u1*(((di+1)*(di+1)*c*c)+(k1*k1)))
                                    JJ0 = JJ0+factor1*factor2

                                factor = 1-u1/np.sqrt(1+u1*u1)
                                i2_u1_k1 = np.exp(-1j*k1*u1)*((2+1j*k1*u1)*factor-u1/((1+u1*u1)*np.sqrt(1+u1*u1))-1j*k1*II0+k1*k1*JJ0)
                                i2_u1_k1 = i2_u1_k1/3
                                II2 = i2_u1_k1

                        else:
                            if u1.real <= bigm:
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*1/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                JJ0 = 0
                                for di in range(n):
                                    factor1 = a[di]*1/((((di+1)*(di+1)*c*c)+(k1*k1))*(((di+1)*(di+1)*c*c)+(k1*k1)))
                                    factor2 = (di+1)*(di+1)*c*c-k1*k1-1j*k1*(2*(di+1)*c)
                                    JJ0 = JJ0+factor1*factor2

                                i2_0_k1 = (2-1j*k1*II0+k1*k1*JJ0)
                                i2_0_k1 = i2_0_k1/3
                                II2 = 2*i2_0_k1
                            else:
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*1/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                JJ0 = 0
                                for di in range(n):
                                    factor1 = a[di]*1/((((di+1)*(di+1)*c*c)+(k1*k1))*(((di+1)*(di+1)*c*c)+(k1*k1)))
                                    factor2 = (di+1)*(di+1)*c*c-k1*k1-1j*k1*(2*(di+1)*c)
                                    JJ0 = JJ0+factor1*factor2

                                i2_0_k1 = (2-1j*k1*II0+k1*k1*JJ0)
                                i2_0_k1 = i2_0_k1/3
                                II0 = 0
                                for di in range(n):
                                    factor = a[di]*np.exp((di+1)*c*u1)/(((di+1)*(di+1)*c*c)+(k1*k1))
                                    II0 = II0+factor*((di+1)*c-1j*k1)

                                JJ0 = 0
                                for di in range(n):
                                    factor1 = a[di]*np.exp((di+1)*c*u1)/((((di+1)*(di+1)*c*c)+(k1*k1))*(((di+1)*(di+1)*c*c)+(k1*k1)))
                                    factor2 = (di+1)*(di+1)*c*c-k1*k1-(di+1)*c*u1*(((di+1)*(di+1)*c*c)+(k1*k1))-1j*k1*(2*(di+1)*c-u1*(((di+1)*(di+1)*c*c)+(k1*k1)))
                                    JJ0 = JJ0+factor1*factor2

                                factor = 1+u1/np.sqrt(1+u1*u1)
                                i2_mu1_k1 = np.exp(1j*k1*u1)*((2-1j*k1*u1)*factor+u1/((1+u1*u1)*np.sqrt(1+u1*u1))-1j*k1*II0+k1*k1*JJ0)
                                i2_mu1_k1 = i2_mu1_k1/3
                                II2 = (2*np.real(i2_0_k1)-np.real(i2_mu1_k1))+(np.imag(i2_mu1_k1))*1j



                        if k == 0:
                            K2 = K2s
                        else:
                            K2 = -3*II2-((1j*k1*M*M*r*r)/(R*R)*(np.exp(-1j*k1*u1)/np.sqrt(1+u1*u1)))-M*r/R*((1+u1*u1)*(beta2*r*r/(R*R))+2+M*r*u1/R)*(np.exp(-1j*k1*u1)/((1+u1*u1)*np.sqrt(1+u1*u1)))



                    vec_KP2[station] = (K2*np.exp(-1j*k*x0/b)-K2s)*T2est


                KL = vec_KP2[0]
                KML = vec_KP2[1]
                KC = vec_KP2[2]
                KMR = vec_KP2[3]
                KR = vec_KP2[4]


                A2 = 2/3*(6*KC-4*KMR-4*KML+KL+KR)/e**4
                B2 = -2/3*(2*KMR-2*KML+KL-KR)/e**3
                C2 = -1/6/e**2*(-16*KMR-16*KML+KL+KR+30*KC)
                D2 = 1/6*(8*KMR-8*KML+KL-KR)/e
                E2 = KC

                if np.abs(zbar) < eps:
                    D_ = (2*e*(3*A2*e**6-e*e*E2-e*e*(4*D2+9*B2*e*e)*ybar-3*(9*A2*e**4+E2)*ybar**2+8*B2*e*e*ybar**3+32*A2*e*e*ybar**4-3*B2*ybar**5-12*A2*ybar**6-C2*e*e*(3*e*e+ybar*ybar))+3*(B2+4*A2*ybar)*(e*e-ybar*ybar)**3*(np.log(-e+ybar)-np.log(e+ybar)))/(3*(e-ybar)**3*(e+ybar)**3)
                    D_ = (deltaxs/(8*np.pi))*D_
                else:
                    D_ = ((2*e*(2*A2*e**4*zbar*zbar+E2*(-ybar*ybar+zbar*zbar)+e*e*(E2+ybar*(D2+ybar*(C2+ybar*(B2+A2*ybar)))-(C2+3*B2*ybar+10*A2*ybar*ybar)*zbar*zbar+5*A2*zbar**4)-(ybar*ybar+zbar*zbar)*(D2*ybar+(ybar*ybar+zbar*zbar)*(C2+ybar*(B2+A2*ybar)-3*A2*zbar*zbar))))/((e*e-ybar*ybar)**2*zbar**2+2*(e*e+ybar*ybar)*zbar**4+zbar**6)+((E2+ybar*(D2+ybar*(C2+ybar*(B2+A2*ybar)))+(C2+3*ybar*(B2+2*A2*ybar))*zbar*zbar-3*A2*zbar**4)*np.arctan((e-ybar)/zbar)+(E2+ybar*(D2+ybar*(C2+ybar*(B2+A2*ybar)))+(C2+3*ybar*(B2+2*A2*ybar))*zbar*zbar-3*A2*zbar**4)*np.arctan((e+ybar)/zbar)+(B2+4*A2*ybar)*zbar**3*(np.log((e-ybar)**2+zbar*zbar)-np.log((e+ybar)**2+zbar*zbar)))/zbar**3)/2
                    D_ = (deltaxs/(8*np.pi))*D_

            D2 = D_

            D[rb,sb] = D[rb,sb]+(D1+D2)

            D[rb,sb] = D[rb,sb]+Ds[rb,sb]
            #print(rb,sb)

    return D


def Drss(x0, y0, z0, gamar, gamas, lambdas, e, deltaxs, paramDLM):
    M = paramDLM.M
    beta2 = 1 - M*M
    beta = np.sqrt(beta2)

    X0 = x0/beta
    tanlambda = np.tan(lambdas)/beta

    R = np.array([X0, y0, z0], dtype=np.complex128)
    ri = -e*np.array([tanlambda, np.cos(gamas), np.sin(gamas)], dtype=np.complex128)
    ro = e*np.array([tanlambda, np.cos(gamas), np.sin(gamas)], dtype=np.complex128)
    Ri = R-ri
    Ro = R-ro

    lamb = np.arctan(tanlambda)

    gamab = np.array([np.sin(lamb), np.cos(lamb)*np.cos(gamas), np.cos(lamb)*np.sin(gamas)], dtype=np.complex128)
    gamai = np.array([-1, 0, 0], dtype=np.complex128)
    gamao = np.array([1, 0, 0], dtype=np.complex128)

    Ri_ = Ri/np.sqrt(Ri[0]**2 + Ri[1]**2 + Ri[2]**2)
    Ro_ = Ro/np.sqrt(Ro[0]**2 + Ro[1]**2 + Ro[2]**2)

    costhetab = Ri_[0]*gamab[0]+Ri_[1]*gamab[1]+Ri_[2]*gamab[2]
    cosphib = Ro_[0]*gamab[0]+Ro_[1]*gamab[1]+Ro_[2]*gamab[2]
    db = Ri-gamab*np.sqrt(Ri[0]**2+Ri[1]**2+Ri[2]**2)*costhetab

    costhetai = 1
    cosphii = Ri_[0]*gamai[0]+Ri_[1]*gamai[1]+Ri_[2]*gamai[2]
    di = np.array([0, Ri[1], Ri[2]], dtype=np.complex128)

    costhetao = Ro_[0]*gamao[0]+Ro_[1]*gamao[1]+Ro_[2]*gamao[2]
    cosphio = -1
    do = np.array([0, Ro[1], Ro[2]], dtype=np.complex128)

    vetorial = np.array([gamab[1]*db[2]-gamab[2]*db[1], gamab[2]*db[0]-gamab[0]*db[2], gamab[0]*db[1]-gamab[1]*db[0]], dtype=np.complex128)

    if np.abs(costhetab) > 0.999 or np.abs(cosphib) > 0.999:
        if costhetab.real*cosphib.real < 0:
            Vb = np.array([0, 0, 0], dtype=np.complex128)
        else:
            Vb = (1/(4*np.pi))*vetorial*(0.5*np.abs(1/(Ri[0]**2+Ri[1]**2+Ri[2]**2))-1/(Ro[0]**2+Ro[1]**2+Ro[2]**2))
    else:
        Vb = (1/(4*np.pi*(db[0]**2+db[1]**2+db[2]**2)))*vetorial*(costhetab-cosphib)

    vetorial = np.array([gamai[1]*di[2]-gamai[2]*di[1],
                         gamai[2]*di[0]-gamai[0]*di[2],
                         gamai[0]*di[1]-gamai[1]*di[0]], dtype=np.complex128)
    modulo2 = di[0]**2+di[1]**2+di[2]**2

    Vi = (1/(4*np.pi*modulo2))*vetorial*(costhetai-cosphii)

    vetorial = np.array([gamao[1]*do[2]-gamao[2]*do[1],
                         gamao[2]*do[0]-gamao[0]*do[2],
                         gamao[0]*do[1]-gamao[1]*do[0]], dtype=np.complex128)
    modulo2 = do[0]**2 + do[1]**2 + do[2]**2

    Vo = (1/(4*np.pi*modulo2))*vetorial*(costhetao-cosphio)

    V = Vb+Vi+Vo

    D = V[1]*np.sin(gamar)-V[2]*np.cos(gamar)
    D = -D*deltaxs/2

    return D
