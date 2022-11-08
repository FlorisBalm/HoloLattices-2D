from scipy.optimize import fsolve
import argparse
import numpy as np
parser=argparse.ArgumentParser()

parser.add_argument("--B", type=float)
parser.add_argument("--T", type=float)
parsed,_ = parser.parse_known_args()

B = parsed.B
T = parsed.T

def FTMU(mu,B,T):
    return (T - (12 - mu**2 - B**2*mu**4)/(16*np.pi*mu))

def FMUSimple(T):
    return 2*(-4*np.pi*T+np.sqrt(3+16*np.pi**2*T**2))

def SolveF(T,B):
    return fsolve(FTMU, FMUSimple(T),(B,T),xtol=1e-9)

print(SolveF(T,B)[0])
