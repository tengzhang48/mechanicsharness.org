import numpy as np
def G0(y):
    return (-2*y - (y**2 + 2*y + 1)*np.log((1 - y)/(y + 1)) + 1j*np.pi*(y**2 + 2*y + 1))/(y**2 + 2*y + 1)

def G1(y):
    return y*(2*y*(2*y + 3) - (np.log(1 - y) - np.log(y + 1))*(y**2 + 2*y + 1) + (-4 + 1j*np.pi)*(y**2 + 2*y + 1))/(y**2 + 2*y + 1)

def G2(y):
    return y**2*(-8*y**2 - 10*y + (np.log(-y**8) + 8)*(y**2 + 2*y + 1) - (y**2 + 2*y + 1)*np.log((1 - y)*(y + 1)**7))/(y**2 + 2*y + 1)

def G3(y):
    return y**2*(12*y**3 + 14*y**2 + y*(y**2 + 2*y + 1)*(-16*np.log(y) - 12 + 1j*np.pi) - (y**2 + 2*y + 1)*(y*np.log(1 - y) - 17*y*np.log(y + 1) + 8))/(y**2 + 2*y + 1)

def G4(y):
    return y**2*(-16*y**4 - 18*y**3 + y**2*(y**2 + 2*y + 1)*(32*np.log(y) + 16 + 1j*np.pi) - (y**2 + 2*y + 1)*(y**2*np.log(1 - y) + 31*y**2*np.log(y + 1) - 16*y + 4))/(y**2 + 2*y + 1)

def G5(y):
    return y**2*(20*y**5 + 22*y**4 + y**3*(y**2 + 2*y + 1)*(-48*np.log(y) - 20 + 1j*np.pi) - 1/3*(y**2 + 2*y + 1)*(3*y**3*np.log(1 - y) - 147*y**3*np.log(y + 1) + 96*y**2 - 24*y + 8))/(y**2 + 2*y + 1)

def G6(y):
    return y**2*(-24*y**6 - 26*y**5 + y**4*(y**2 + 2*y + 1)*(72*np.log(y) + 24 + 1j*np.pi) - 1/3*(y**2 + 2*y + 1)*(3*y**4*np.log(1 - y) + 213*y**4*np.log(y + 1) - 144*y**3 + 48*y**2 - 16*y + 6))/(y**2 + 2*y + 1)

def G7(y):
    return y**2*(28*y**7 + 30*y**6 + y**5*(y**2 + 2*y + 1)*(-96*np.log(y) - 28 + 1j*np.pi) - 1/15*(y**2 + 2*y + 1)*(15*y**5*np.log(1 - y) - 1455*y**5*np.log(y + 1) + 1080*y**4 - 360*y**3 + 160*y**2 - 60*y + 24))/(y**2 + 2*y + 1)

