"""
-----------------------------------------------------------------------
Function for discretizing lifting surfaces (half on the right)

Author: Higor Luis Silva
-----------------------------------------------------------------------
"""
import numpy as np


class Box(object):
    def __init__(self, N):
        self.tipo = np.empty(N, dtype=np.int32)
        self.x14i = np.empty(N, dtype=np.complex128)
        self.y14i = np.empty(N, dtype=np.complex128)
        self.z14i = np.empty(N, dtype=np.complex128)
        self.x14o = np.empty(N, dtype=np.complex128)
        self.y14o = np.empty(N, dtype=np.complex128)
        self.z14o = np.empty(N, dtype=np.complex128)
        self.x14ci = np.empty(N, dtype=np.complex128)
        self.y14ci = np.empty(N, dtype=np.complex128)
        self.z14ci = np.empty(N, dtype=np.complex128)
        self.x14co = np.empty(N, dtype=np.complex128)
        self.y14co = np.empty(N, dtype=np.complex128)
        self.z14co = np.empty(N, dtype=np.complex128)
        self.x14c = np.empty(N, dtype=np.complex128)
        self.y14c = np.empty(N, dtype=np.complex128)
        self.z14c = np.empty(N, dtype=np.complex128)
        self.x34c = np.empty(N, dtype=np.complex128)
        self.y34c = np.empty(N, dtype=np.complex128)
        self.z34c = np.empty(N, dtype=np.complex128)
        self.inbLEx = np.empty(N, dtype=np.complex128)
        self.inbLEy = np.empty(N, dtype=np.complex128)
        self.inbLEz = np.empty(N, dtype=np.complex128)
        self.inbTEx = np.empty(N, dtype=np.complex128)
        self.inbTEy = np.empty(N, dtype=np.complex128)
        self.inbTEz = np.empty(N, dtype=np.complex128)
        self.outbLEx = np.empty(N, dtype=np.complex128)
        self.outbLEy = np.empty(N, dtype=np.complex128)
        self.outbLEz = np.empty(N, dtype=np.complex128)
        self.outbTEx = np.empty(N, dtype=np.complex128)
        self.outbTEy = np.empty(N, dtype=np.complex128)
        self.outbTEz = np.empty(N, dtype=np.complex128)
        self.dihedral = np.empty(N, dtype=np.complex128)
        self.lamb = np.empty(N, dtype=np.complex128)
        self.lambdas = np.empty(N, dtype=np.complex128)
        self.chord = np.empty(N, dtype=np.complex128)
        self.normal = np.empty((N, 3), dtype=np.complex128)
        self.area = np.empty(N, dtype=np.complex128)


def discretize(wing, paramDLM):
    xboxes, yboxes = paramDLM.xboxes, paramDLM.yboxes
    apex = paramDLM.apex
    inbLE = paramDLM.inbLE
    inbTE = paramDLM.inbTE
    outbTE = paramDLM.outbTE
    outbLE = paramDLM.outbLE

    totalboxes = 0

    if wing.tipo == 0:
        totalboxes = totalboxes + (xboxes*yboxes)
    elif wing.tipo == 1:
        totalboxes = totalboxes + (xboxes*yboxes)
    elif wing.tipo == 2:
        totalboxes = totalboxes + (xboxes*yboxes)

    tot = []
    TOTAL = 0

    totx = xboxes
    toty = yboxes

    xa = apex[0]; ya = apex[1]; za = apex[2]
    x1 = inbLE[0]; y1 = inbLE[1]; z1 = inbLE[2]
    x2 = inbTE[0]; y2 = inbTE[1]; z2 = inbTE[2]
    x3 = outbTE[0]; y3 = outbTE[1]; z3 = outbTE[2]
    x4 = outbLE[0]; y4 = outbLE[1]; z4 = outbLE[2]
    tot = xboxes*yboxes

    box = Box(TOTAL+tot)

    for j in range(TOTAL+tot):
        if wing.tipo != 2:
            dy = (y4-y1)/toty
            jj = j
            ui = jj-TOTAL
            box.tipo[jj] = wing.tipo
            passox = np.remainder(ui,totx)
            passoy = np.floor(ui/totx)
            yi = passoy*dy
            af = babf(yi, paramDLM)
            box.x14i[jj] = af[0]+passox*(af[1]-af[0])/totx+(af[1]-af[0])/totx/4
            box.y14i[jj] = ya+yi
            box.z14i[jj] = za+(yi-y1)/(y4-y1)*(z4-z1)

            yo = (passoy+1)*dy
            af = babf(yo, paramDLM)
            box.x14o[jj] = af[0]+passox*(af[1]-af[0])/totx+(af[1]-af[0])/totx/4
            box.y14o[jj] = ya+yo
            box.z14o[jj] = za+(yo-y1)/(y4-y1)*(z4-z1)

            yci = yi+(yo-yi)/4
            af = babf(yci, paramDLM)
            box.x14ci[jj] = af[0]+passox*(af[1]-af[0])/totx+(af[1]-af[0])/totx/4
            box.y14ci[jj] = ya+yci
            box.z14ci[jj] = za+(yci-y1)/(y4-y1)*(z4-z1)

            yco = yi+3*(yo-yi)/4
            af = babf(yco, paramDLM)
            box.x14co[jj] = af[0]+passox*(af[1]-af[0])/totx+(af[1]-af[0])/totx/4
            box.y14co[jj] = ya+yco
            box.z14co[jj] = za+(yci-y1)/(y4-y1)*(z4-z1)


            yc = yi+(yo-yi)/2
            af = babf(yc, paramDLM)
            box.x14c[jj] = af[0]+passox*(af[1]-af[0])/totx+(af[1]-af[0])/totx/4
            box.y14c[jj] = ya+yc
            box.z14c[jj] = za+(yc-y1)/(y4-y1)*(z4-z1)

            box.x34c[jj] = af[0]+passox*(af[1]-af[0])/totx+3*(af[1]-af[0])/totx/4
            box.y34c[jj] = ya+yc
            box.z34c[jj] = za+(yc-y1)/(y4-y1)*(z4-z1)

            af = babf(yi, paramDLM)
            box.inbLEx[jj] = af[0]+passox*(af[1]-af[0])/totx
            box.inbLEy[jj] = ya+yi
            box.inbLEz[jj] = za+(yi-y1)/(y4-y1)*(z4-z1)
            box.inbTEx[jj] = af[0]+passox*(af[1]-af[0])/totx+(af[1]-af[0])/totx
            box.inbTEy[jj] = ya+yi
            box.inbTEz[jj] = za+(yi-y1)/(y4-y1)*(z4-z1)

            af = babf(yo, paramDLM)
            box.outbTEx[jj] = af[0]+passox*(af[1]-af[0])/totx+(af[1]-af[0])/totx
            box.outbTEy[jj] = ya+yo
            box.outbTEz[jj] = za+(yo-y1)/(y4-y1)*(z4-z1)
            box.outbLEx[jj] = af[0]+passox*(af[1]-af[0])/totx
            box.outbLEy[jj] = ya+yo
            box.outbLEz[jj] = za+(yo-y1)/(y4-y1)*(z4-z1)

            box.dihedral[jj] = np.arctan((z4-z1)/(y4-y1))
            box.lamb[jj] = np.arctan((x4-x1)/(np.sqrt((y4-y1)**2+(z4-z1)**2)))
            box.lambdas[jj] = np.arctan((box.x14o[jj]-box.x14i[jj])/\
                      np.sqrt((box.y14o[jj]-box.y14i[jj])**2+(box.z14o[jj]-box.z14i[jj])**2))
            box.chord[jj] = 2*(box.x34c[jj]-box.x14c[jj])

            normal = -np.cross([box.inbTEx[jj]-box.inbLEx[jj],box.inbTEy[jj]-box.inbLEy[jj],box.inbTEz[jj]-box.inbLEz[jj]],[box.outbLEx[jj]-box.inbLEx[jj],box.outbLEy[jj]-box.inbLEy[jj],box.outbLEz[jj]-box.inbLEz[jj]])
            box.normal[jj] = normal/np.linalg.norm(normal)

            normal2 = np.cross([box.inbTEx[jj]-box.outbTEx[jj],box.inbTEy[jj]-box.outbTEy[jj],box.inbTEz[jj]-box.outbTEz[jj]],[box.outbLEx[jj]-box.outbTEx[jj],box.outbLEy[jj]-box.outbTEy[jj],box.outbLEz[jj]-box.outbTEz[jj]])
            box.area[jj] = 0.5*np.linalg.norm(normal)+0.5*np.linalg.norm(normal2)

    TOTAL = TOTAL + tot
    totalboxes = TOTAL
    return totalboxes, box


def babf(y, paramDLM):
    apex = paramDLM.apex
    inbLE = paramDLM.inbLE
    inbTE = paramDLM.inbTE
    outbTE = paramDLM.outbTE
    outbLE = paramDLM.outbLE

    # This function determines the coordinates of the leading and
    # trailing edge for each station
    xa = apex[0];   ya = apex[1];   za = apex[2]
    x1 = inbLE[0];  y1 = inbLE[1];  z1 = inbLE[2]
    x2 = inbTE[0];  y2 = inbTE[1];  z2 = inbTE[2]
    x3 = outbTE[0]; y3 = outbTE[1]; z3 = outbTE[2]
    x4 = outbLE[0]; y4 = outbLE[1]; z4 = outbLE[2]

    if y4 != y1:
        xLE = xa+x1+(x4-x1)*(y-y1)/(y4-y1)
        xTE = xa+x2+(x3-x2)*(y-y2)/(y3-y2)
    else:
        xLE = xa+x1+(x4-x1)*(y-z1)/(z4-z1)
        xTE = xa+x2+(x3-x2)*(y-z2)/(z3-z2)

    babf = [xLE, xTE]

    return babf
