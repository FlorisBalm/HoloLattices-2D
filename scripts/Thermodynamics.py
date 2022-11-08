#!/usr/bin/env python3
import sys
import h5py
import numpy as np

def FullThermodynamics(filename):

    with h5py.File(filename, "r") as infile:
        data = np.transpose(infile["result"][:], (3,2,1,0))
        params = infile["parameters"][:]

    mu = params[0]
    mu1=mu
    lx = params[4]
    ly = params[5]
    B = params[-2]
    dx=1/(data.shape[3]-1)
    #mu=mu1=0.861422
    
    dr   = (1./(12*dx))*(data[:,:,:,-5]*3 -data[:,:,:,-4]*16 + 36* data[:,:,:,-3] - 48*data[:,:,:,-2]+25*data[:,:,:,-1])

    dr3   = (1./(2*dx**3))*(data[:,:,:,-5]*3 -data[:,:,:,-4]*14 + 24* data[:,:,:,-3] - 18*data[:,:,:,-2]+5*data[:,:,:,-1])


    dr_bdy   = (-147*data[:,:,:,0]+360*data[:,:,:,1]-450*data[:,:,:,2]+400*data[:,:,:,3]-225*data[:,:,:,4]+72*data[:,:,:,5]-10*data[:,:,:,6])/(60*1.0*dx**1)

    dr3_bdy   = (((-17/4)*data[:,:,:,0]+(71/4)*data[:,:,:,1]-(59/2)*data[:,:,:,2]+(49/2)*data[:,:,:,3]-(41/4)*data[:,:,:,4] + (7/4)*data[:,:,:,5]))/(dx**3)


    #0    1   2   3    4    5   6   7    8    9   10  11   12  13  14  
    #{at, ax, ay, az, psi, Qtt, M, Q, Qzz, Qxz, Qyz, Qtx, Qty, Qtz,R}
    Txx =np.average( (1/mu)**3* (1 + (mu**2+B**2)/4 + (3*dr3_bdy[7,:,:] + dr3_bdy[6,:,:])/(8*6)))
    Tyy =np.average( (1/mu)**3* (1 + (mu**2+B**2)/4 + (3*dr3_bdy[7,:,:] - dr3_bdy[6,:,:])/(8*6)))

    Ttt =np.average( (1/mu)**3* (2 + mu**2/2 - 3*dr3_bdy[5,:,:]/(8*6)))

    Txy =np.average( (1/mu)**3* ((3*ly)/(4*lx)*dr3_bdy[14,:,:]))

    Rho = np.average((1/mu)*(data[0,:,:,0]-0.5*dr_bdy[0,:,:]))


    S = np.average((4*np.pi*data[7, :,:, -1]/mu**2))
    
    T = (12 - mu*mu - B*B)/(16*np.pi*mu)
    Mu = mu*data[0,:,:,0]

    InternalEnergy = np.average(((2 + (mu)**2/2 - 1*dr3_bdy[5,:,:]/16)/(mu**3)))
    FreeEnergyInUnitsOfMu= InternalEnergy - (T*S) - (np.average(Rho))
    
    return {"rho":Rho, "S":S, "T":T, "Mu": np.average(Mu), "Omega": (FreeEnergyInUnitsOfMu), "E":(InternalEnergy), "Txx": Txx, "Tyy":Tyy, "Ttt": Ttt, "Txy":Txy}

#
def ThermodynamicPotential(filename,printPots = False):
    return FullThermodynamics(filename, printPots)["Omega"]

if __name__ == '__main__':
    import argparse
    parser  = argparse.ArgumentParser()

    parser.add_argument("--Print", action="store_true", default=False)
    parser.add_argument("--output",type=str,  default=False)
    parser.add_argument("--omega", action="store_true", default=False)

    parser.add_argument("--files", type=str, nargs="+")
    parsed, _ = parser.parse_known_args()

    if parsed.omega:
        td=FullThermodynamics(parsed.files[0])
        print(td["Omega"])
        exit()

    T=P=A=0
    import pandas
    df = pandas.DataFrame()
    for arg in parsed.files:
        with h5py.File(arg, "r") as infile:
            params = infile["parameters"][:]
            mu=params[0]
            T = (12-mu**2)/(16*np.pi*mu)
            P=(2*np.pi*params[6])/(params[4]*mu)
            A=params[2]
        ftdDict = {k: np.array([v]) for k,v in FullThermodynamics(arg).items()}
        ftdDict["A"] = np.array([A])
        ftdDict["P"] = np.array([P])
        ftdDict["T"] = np.array([T])
        ftdDf = pandas.DataFrame(ftdDict)

        df = pandas.concat((df, ftdDf))


    if parsed.Print:
        from io import StringIO
        output=StringIO()

    else:
        output=parsed.output
    
    df.to_csv(output, sep="\t", index=False)
    if parsed.Print:
        output.seek(0)
        print(output.read())



