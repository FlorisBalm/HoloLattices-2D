#!/usr/bin/env python3
import sys
import h5py
import numpy as np
import yaml

def centred_finite_difference_weights(points, deriv_centre, max_deriv=2):                                                               
    """                                                                                                                                 
    Gets the finite difference weights for derivatives on the set of gridpoints, to a maximum                                           
    derivative order of max_deriv, centred at point x0. This point x0 can even be a non-gridpoint                                       
    Credit for this goes to some paper I found. TODO: find this paper for proper credit 
    :param points: Locations of gridpoints to evaluate derivative over               
    :param max_deriv: maximum order of derivative to take                            
    :param deriv_centre: CCentre of the derivative                                   
    :return: 0th, 1st, .., max_deriv derivatives on the gridpoints                   
    """                                                                              
    N = len(points) - 1                                                              
    c1 = 1.0                                                                         
    ret = [[0 for _ in range(N + 1)] for _ in range(max_deriv + 1)]                  
    c4 = points[0] - deriv_centre                                                    
    ret[0][0] = 1.0                                                                  
                                                                                     
    for i in range(1, N + 1):                                                        
        mn = min(i, max_deriv)                                                       
        c2 = 1.0                                                                     
        c5 = c4                                                                      
        c4 = points[i] - deriv_centre                                                
                                                                                     
        for j in range(i):                                                           
            c3 = points[i] - points[j]                                               
            c2 = c2 * c3                                                             
                                                                                     
            if j + 1 == i:                                                           
                for k in np.arange(mn, 0, -1):                                       
                    ret[k][i] = c1 * (k * ret[k - 1][i - 1] - c5 * ret[k][i - 1]) / c2
                ret[0][i] = -c1 * c5 * ret[0][i - 1] / c2                             
                                                                                      
            for k in np.arange(mn, 0, -1):                                            
                ret[k][j] = (c4 * ret[k][j] - k * ret[k - 1][j]) / c3                 
            ret[0][j] = c4 * ret[0][j] / c3                                           
        c1 = c2                                                                       
    return ret       

def FullThermodynamics(filename,extra_order=1):

    with h5py.File(filename, "r") as infile:
        data = infile["result"][:].transpose((3,2,1,0))
        params = infile["parameters"][:]

    mu= params[0]
    Q = params[1]
    ax = params[2]
    ay = params[3]
    G = params[4]
    npx = params[6]

    z_dimension = data.shape[-1]

    z_coords=np.linspace(0,1,z_dimension)

    z_coords = (1-np.cos(np.pi*z_coords))/2

         
    diff_order=6



    extra_spacing = extra_order

    dr_b_stencil = centred_finite_difference_weights(z_coords[:extra_spacing*diff_order:extra_spacing],z_coords[0], 3)
    dr_h_stencil = centred_finite_difference_weights(z_coords[-extra_spacing*diff_order::extra_spacing],z_coords[-1], 3)
    
    val_h = np.einsum("ijkl,l", data[:,:,:,-extra_spacing*diff_order::extra_spacing], dr_h_stencil[0])
    dr1_h = np.einsum("ijkl,l", data[:,:,:,-extra_spacing*diff_order::extra_spacing], dr_h_stencil[1])
    dr2_h = np.einsum("ijkl,l", data[:,:,:,-extra_spacing*diff_order::extra_spacing], dr_h_stencil[2])
    dr3_h = np.einsum("ijkl,l", data[:,:,:,-extra_spacing*diff_order::extra_spacing], dr_h_stencil[3])
    
    val_bdy = np.einsum("ijkl,l", data[:,:,:,:extra_spacing*diff_order:extra_spacing], dr_b_stencil[0])
    dr1_bdy = np.einsum("ijkl,l", data[:,:,:,:extra_spacing*diff_order:extra_spacing], dr_b_stencil[1])
    dr2_bdy = np.einsum("ijkl,l", data[:,:,:,:extra_spacing*diff_order:extra_spacing], dr_b_stencil[2])
    dr3_bdy = np.einsum("ijkl,l", data[:,:,:,:extra_spacing*diff_order:extra_spacing], dr_b_stencil[3])
    
    #0    1     2   3    4    5  6   
    #{at, phi Qtt,Qxx,Qyy, Qzz, Qxz}
    Ttt = (-2*(1+Q)**3 + 0.5*dr3_bdy[2] + 0.75*Q**2*dr1_bdy[1])/mu**3
    Txx = (   (1+Q)**3 + 0.5*dr3_bdy[3] + 0.75*Q**2*dr1_bdy[1])/mu**3
    Tyy = (   (1+Q)**3 + 0.5*dr3_bdy[4] + 0.75*Q**2*dr1_bdy[1])/mu**3

    EInternal = -Ttt

    TrTtt = Ttt+Txx+Tyy
    Txx2nd= dr2_bdy[3]
    
    Rho = (((1+Q)*val_bdy[0]-dr1_bdy[0])/(np.sqrt(3*Q*(1+Q))))
    S = (4*np.pi*np.sqrt((1+Q)**3*(1+val_h[3])*(1+val_h[4])))/mu**2
    
    T = np.sqrt(3)/(4*np.pi*np.sqrt(Q))
    Mu = val_bdy[0]

    Omega = EInternal - Mu*Rho - T*S
    
    return {
            "EInternal":np.average((EInternal)), "Txx":np.average(Txx), "Tyy":np.average(Tyy),
            "rho":np.average(Rho), "S":np.average(S), "T":T, "Mu": np.average(Mu),
            "Trace":np.average(TrTtt), "TRatio":np.average(TrTtt)/np.average(EInternal), "Omega":np.average(Omega),
            "murh":mu, "ax": ax, "Gmu":G/mu,"ay": ay, "2nd":np.average(Txx2nd)
            }

def ThermodynamicPotential(filename,printPots = False):
    return FullThermodynamics(filename, printPots)["Omega"]

if __name__ == '__main__':
    import argparse
    parser  = argparse.ArgumentParser()

    parser.add_argument("--output",type=str,  default=False)
    parser.add_argument("--omega", action="store_true", default=False)

    parser.add_argument("--files", type=str, nargs="+")
    parser.add_argument("--extra_order", type=int)
    parsed, _ = parser.parse_known_args()

    import pandas
    df = pandas.DataFrame()

    for arg in parsed.files:

        ftdDict = {k: np.array([v]) for k,v in FullThermodynamics(arg,extra_order=parsed.extra_order).items()}
        ftdDf = pandas.DataFrame(ftdDict)

        df = pandas.concat((df, ftdDf))

    if parsed.output:
        # This is to skip the first line
        output=parsed.output
    else:
        from io import StringIO
        output=StringIO()
    
    df.to_csv(output, sep="\t", index=False)

    if not parsed.output:
        output.seek(0)
        print(output.read())



