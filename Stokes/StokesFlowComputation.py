import logging
import argparse

import pandas

import tempfile
import re

import pathos.multiprocessing as multiprocess
# This makes warnings be sent to the log file instead of to the console output.
logging.captureWarnings(True)

import h5py

from HolographicLattices.DifferentialOperators.DifferentialOperator import *

from HolographicLattices.LinearSolver.SimpleLinearSolver import SimpleLinearSolver

from HolographicLattices.Options import Options
from HolographicLattices.SystemOnBackground import JacobianOnBackground
from HolographicLattices.SystemOnBackground.BackgroundOptions import BackgroundOptions
from HolographicLattices.SystemOnBackground.PerturbationObservables import PerturbationObservables
from HolographicLattices.Utilities import GridUtilities
from ThermodynamicsStokes import FullThermodynamics


def SolveLinearBackground(slice_background=None, **kwargs):
    logger = logging.getLogger(__name__)
    # Basic argparse for some settings. Can be extended, but not required.

    # Basic argparse for some settings. Can be extended, but not required.
    parser = argparse.ArgumentParser()

    # Need to choose one of the two, mandatory!
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--sparse", action="store_true",
                     help="approximate any spectral matrix with its finite difference counterpart in the jacobian")
    grp.add_argument("--dense", action="store_true",
                     help="Use fully spectral matrices for the jacobian. This may increase runtimes a lot!")

    parsed_options, _ = parser.parse_known_args()

    # Load all the options
    setupfile = kwargs["setupfile"]

    io_options = Options.IOOptions.load_from_file(setupfile)
    equation_options = Options.EquationOptions.load_from_file(setupfile)

    Lx = kwargs.get("Lx", 1)
    Ly = kwargs.get("Ly", 1)

    equation_options.grid_domains[0] = Lx
    equation_options.grid_domains[1] = Ly

    background_options = BackgroundOptions.load_from_file(setupfile)
    constant_options = Options.ConstantOptions.load_from_file(setupfile)
    for k, v in kwargs.items():
        if k in constant_options.constants:
            constant_options.constants[k] = v

    grid_information = GridUtilities.GridInformation(equation_options)

    if parsed_options.dense:
        differential_operator = MixedFDDSpectralDifferentialOperator.construct_differentiation_matrices(
            equation_options=equation_options, max_deriv=2
        )
        linear_solver = SimpleLinearSolver()
    else:
        assert parsed_options.sparse, "Somehow the arguments ot mixed up. This is bad, check this."
        differential_operator = MixedFDDSparseSpectralDifferentialOperator.construct_differentiation_matrices(
            equation_options=equation_options, max_deriv=2
        )
        import scipy
        import scipy.sparse.linalg
        linear_solver = SimpleLinearSolver()

    jacobian = JacobianOnBackground.JacobianOnBackground(constant_options=constant_options,
                                                         equation_options=equation_options,
                                                         background_options=background_options,
                                                         finite_difference_matrices=
                                                         differential_operator.get_matrix_representation())

    jacobian.load_from_folder(folder=io_options.coefficient_folder, periodicities=equation_options.grid_periodic)

    grid_volume = np.prod(equation_options.grid_sizes)

    # Only pick horizon
    if slice_background is not None:
        logger.info(f"Slicing background using {slice_background}")
        background_options.slice_background(slice_background)

    # This reshapes the field into what is expected of it for the rest of the computation
    background_options.background_field = background_options.background_field.reshape(
        (background_options.background_field.shape[0],
         np.prod(background_options.background_field.shape[1:])))

    bg_fields_derivatives = differential_operator(background_options.background_field)

    fields_zero = np.zeros((equation_options.num_eqs_of_motion, grid_volume), dtype=equation_options.field_dtype)

    fields_derivatives_zero = differential_operator(fields_zero)

    PDE_Linear = jacobian.evaluate(fields_derivatives_zero, bg_fields_derivatives, grid_information)
    rhs = jacobian.evaluate_rhs(fields_derivatives_zero, bg_fields_derivatives, grid_information)

    solution = linear_solver.solve(PDE_Linear, rhs)

    fields_and_derivatives = differential_operator(solution)

    observables = PerturbationObservables.parse_observables(io_options.observables_file,
                                                            constant_options=constant_options,
                                                            background_options=background_options,
                                                            equation_options=equation_options)


    observables_eval = observables.evaluate(fields_and_derivatives, bg_fields_derivatives, grid_information)


    return {**{i: np.average(j) for i, j in observables_eval.items()},
            **constant_options.constants,
            "filename": background_options.background_filename}

def run_2d(background, setup_file, output, **kwargs):

    logger = logging.getLogger(__name__)
    configuration = {
        "background": background,
        "setupfile": setup_file,
        "output": output,
        # this tells the background to only use the z = 1 data (i.e., the slice through the data at Nz-1)
        # which is required for the stokes flow computation. It makes sense to do here as this is relevant
        "slice_background": (slice(None), slice(None), slice(None), slice(-1, None)),
        "Lx":kwargs.get("Lx", 1),
        "Ly":kwargs.get("Ly", 1)
        }

    mixed_basis=False
    if mixed_basis:
        basis_0 = np.sqrt(0.5) * np.array([1,  1,  0,  0])
        basis_1 = np.sqrt(0.5) * np.array([1, -1,  0, -0])
        basis_2 = np.sqrt(0.5) * np.array([0,  0,  1,  1])
        basis_3 = np.sqrt(0.5) * np.array([0,  0,  1, -1])
        basis_mat = np.array([basis_0, basis_1, basis_2, basis_3])
    else:
        basis_mat = np.identity(4)
        basis_0 = basis_mat[0]
        basis_1 = basis_mat[1]
        basis_2 = basis_mat[2]
        basis_3 = basis_mat[3]

    ex, ey, zetax, zetay = basis_0

    result_basis0 = SolveLinearBackground(write_to_file=False,
                                          read_from_cmd=False,
                                          ex=ex, ey=ey, zetax=zetax, zetay=zetay,
                                          **configuration)

    ex, ey, zetax, zetay = basis_1
    result_basis1 = SolveLinearBackground(write_to_file=False,
                                          ex=ex, ey=ey, zetax=zetax, zetay=zetay,
                                          **configuration)

    ex, ey, zetax, zetay = basis_2
    result_basis2 = SolveLinearBackground(write_to_file=False,
                                          ex=ex, ey=ey, zetax=zetax, zetay=zetay,
                                          **configuration)

    ex, ey, zetax, zetay = basis_3
    result_basis3 = SolveLinearBackground(write_to_file=False,
                                          ex=ex, ey=ey, zetax=zetax, zetay=zetay,
                                          **configuration)

    # Choice of which result does not matter - all get the same name
    Thermodynamics = FullThermodynamics(result_basis0["filename"])

    with h5py.File(result_basis0["filename"], "r") as infile:
        data = infile["result"][:]
        params = infile["parameters"][:]
        mu=params[0]

    thermoelectric_conductivity_full = np.array([
        [result_basis0["Jx"], result_basis1["Jx"], result_basis2["Jx"], result_basis3["Jx"]],
        [result_basis0["Jy"], result_basis1["Jy"], result_basis2["Jy"], result_basis3["Jy"]],
        [result_basis0["Qx"], result_basis1["Qx"], result_basis2["Qx"], result_basis3["Qx"]],
        [result_basis0["Qy"], result_basis1["Qy"], result_basis2["Qy"], result_basis3["Qy"]]
    ]).dot(
        np.linalg.inv(
            np.array(
                basis_mat).transpose()))
    sigma_result = thermoelectric_conductivity_full[0:2, 0:2]
    talpha_result = thermoelectric_conductivity_full[0:2, 2:4] / mu
    talphabar_result = thermoelectric_conductivity_full[2:4, 0:2] / mu
    tkappa_result = thermoelectric_conductivity_full[2:4, 2:4] / mu ** 2
    rhoH = result_basis0["rhoH"]/mu**2
    S = result_basis0["Srh"]/mu**2


    #rho = np.average(Thermodynamics["rho"])
    #T = np.average(Thermodynamics["T"])
    #S = np.average(Thermodynamics["S"])

    Ttt = np.average(Thermodynamics["Ttt"])

    Txx = np.average(Thermodynamics["Txx"])
    Tyy = np.average(Thermodynamics["Tyy"])

    #Txy =  np.average(Thermodynamics["Txy"]) # Txy is symmetric
    #Tyx =  np.average(Thermodynamics["Tyx"]) # Txy is symmetric

    #Ttx =  np.average(Thermodynamics["Ttx"]) # Ttx is symmetric
    #Tty =  np.average(Thermodynamics["Tty"]) # Tty is symmetric

    #mu0 = mu
    #Mu = np.average(Thermodynamics["Mu"]) / mu0

    #SigmaQ0 = sigma_result - talpha_result @ np.linalg.inv(tkappa_result) @ talphabar_result

    data_ret = {
        "SigmaE11": sigma_result[0, 0], "SigmaE12": sigma_result[0, 1], "SigmaE21": sigma_result[1, 0],
        "SigmaE22": sigma_result[1, 1],
        "SigmaAlpha11": talpha_result[0, 0], "SigmaAlpha12": talpha_result[0, 1], "SigmaAlpha21": talpha_result[1, 0],
        "SigmaAlpha22": talpha_result[1, 1],
        "SigmaAlphaBar11": talphabar_result[0, 0], "SigmaAlphaBar12": talphabar_result[0, 1],
        "SigmaAlphaBar21": talphabar_result[1, 0], "SigmaAlphaBar22": talphabar_result[1, 1],
        "SigmaKappa11": tkappa_result[0, 0], "SigmaKappa12": tkappa_result[0, 1], "SigmaKappa21": tkappa_result[1, 0],
        "SigmaKappa22": tkappa_result[1, 1],
       "rhoH": rhoH,
       "S": S,
       "Ttt": Ttt, "Txx": Txx, "Tyy": Tyy
        #"T": T,
        #"S": S,
        #"Omega": np.average(Thermodynamics["Omega"]),
        #"Pressure": np.average(Thermodynamics["Pressure"]),
        #"EInternal": np.average(Thermodynamics["E"])}
        }

    return {k: np.array([v]) for k, v in data_ret.items()}



def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite_file", action="store_true", default=False)
    parser.add_argument("--setup_replace", type=str)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--log", type=str, default="INFO")
    parser.add_argument("--logfile", type=str, default="out.log")
    parser.add_argument("--backgrounds", metavar='N', type=str, nargs='+',
                        help='all backgrounds to run for')
    parser.add_argument("--ncores", type=int, default=1)
    parsed, _ = parser.parse_known_args()

    logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", filename=parsed.logfile,
                        level=parsed.log)
    logger = logging.getLogger(__name__)
    dataframe = pandas.DataFrame()
    data = None
    logger.info(f"Running using files files: {parsed.backgrounds}")
    with open(parsed.setup_replace, "r") as sources:
        lines = sources.readlines()

    class SolverObj():
        def __init__(self, setup_lines):
            self.lines = setup_lines

        def __call__(self, background):
            # with to have file deleted when the run is done
            with tempfile.NamedTemporaryFile(mode="w") as temp_file:
                # fix the background
                for line in self.lines:
                    temp_file.write(re.sub("REPLACE_BG", background, line))
                #    print(line)
    
                with h5py.File(background, "r") as infile_par:
                    parameters = infile_par["parameters"][:]
    
                    mu = parameters[0]
                    Q = parameters[1]
                    T = np.sqrt(3)/(4*np.pi*np.sqrt(Q))
                    Ax = parameters[2]
                    Ay = parameters[3]
                    Gx = parameters[4]/mu
                    Gy = parameters[5]/mu
                    nperiods = parameters[6]
    
                    Lx = 2*np.pi*nperiods/(Gx*mu)
                    Ly = 2*np.pi*nperiods/(Gy*mu)
    
                temp_file.flush()
    
                data = run_2d(setup_file=temp_file.name, output=None, background=background, Lx=Lx, Ly=Ly)

                data["T"]  = np.array([T ])
                data["Ax"] = np.array([Ax])
                data["Ay"] = np.array([Ay])
                data["Gx"] = np.array([Gx])
                data["Gy"] = np.array([Gy])
                data["mu"] = np.array([mu])
                data["Q"]  = np.array([Q ])
                print("step done")

                return pandas.DataFrame(data)

    solver_obj = SolverObj(lines)
    import sys
    print("cpu count is ", multiprocess.cpu_count())

    with multiprocess.Pool(parsed.ncores) as pool:
        result = pool.map(solver_obj, parsed.backgrounds)

    dataframe = pandas.concat(result)
    dataframe.reset_index(drop=True, inplace=True)


    import os

    if os.path.isfile(parsed.output_file) and not parsed.overwrite_file:
        import uuid

        output_csv_filename = str(uuid.uuid4())
        print(f"Warning! File already exists. Saving to {output_csv_filename}")
    else:
        output_csv_filename = parsed.output_file
    print(dataframe)
    print("Saving to file: ", output_csv_filename)
    dataframe.to_csv(output_csv_filename, sep="\t", index=False)

if __name__ == "__main__":
    main()
