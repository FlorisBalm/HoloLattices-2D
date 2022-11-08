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
from util.ThermodynamicsPetsc import FullThermodynamics


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

def run_2d(background, setup_file, output, basis_angle=0, **kwargs):

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

    cosT = np.cos(basis_angle)
    sinT = np.sin(basis_angle)

    rot_mat = np.array([[cosT, -sinT],[sinT, cosT]])

    basis_2d = np.array([[1,0],[0,1]])

    basis_zero = np.zeros(2)

    basis_2d_rot_x = rot_mat.dot(basis_2d[0])
    basis_2d_rot_y = rot_mat.dot(basis_2d[1])

    basis_elec = np.hstack((basis_2d_rot_x, basis_zero))

    basis_therm = np.hstack((basis_zero ,basis_2d_rot_x)) 
    ex, ey, zetax, zetay = basis_elec
    result_elec = SolveLinearBackground(write_to_file=False,

                                          read_from_cmd=False,
                                          ex=ex, ey=ey, zetax=zetax, zetay=zetay,
                                          **configuration)

    ex, ey, zetax, zetay = basis_therm
    result_therm = SolveLinearBackground(write_to_file=False,

                                          ex=ex, ey=ey, zetax=zetax, zetay=zetay,
                                          **configuration)


    # Choice of which result does not matter - all get the same name
    Thermodynamics = FullThermodynamics(result_elec["filename"])

    with h5py.File(result_elec["filename"], "r") as infile:
        data = infile["result"][:]
        params = infile["parameters"][:]

    mu = params[0]

    thermoelectric_conductivity_full = np.array([
        [result_elec["Jx"], result_elec["Jy"]],
        [result_elec["Qx"], result_elec["Qy"]],
        [result_therm["Jx"], result_therm["Jy"]],
        [result_therm["Qx"], result_therm["Qy"]]
    ])

    cond_result = np.zeros_like(thermoelectric_conductivity_full)

    for i in range(4):
        inner_prod = np.dot(thermoelectric_conductivity_full[i], basis_2d_rot_x)
        CondLong = inner_prod*basis_2d_rot_x
        CondPerp = thermoelectric_conductivity_full[i] - CondLong
        cond_result[i] = np.array([np.linalg.norm(CondLong), np.linalg.norm(CondPerp)])

    cond_result[1:] = cond_result[1:]/mu
    cond_result[3] = cond_result[3]/mu

    sigma_result = cond_result[0]
    talpha_result =  cond_result[1]
    talphabar_result =  cond_result[2]
    tkappa_result =  cond_result[3]

    rhoH = result_elec["rhoH"]/mu**2
    S = result_elec["Srh"]/mu**2

    td_observables = {k:np.average(v) for k,v in Thermodynamics.items()}

    data_ret = {

        "SigmaEL": sigma_result[0], "SigmaET": sigma_result[1],
        "SigmaAlphaL": talpha_result[0], "SigmaAlphaT": talpha_result[1],
        "SigmaAlphaBarL": talphabar_result[0], "SigmaAlphaBarT": talphabar_result[1],
        "SigmaKappaL": tkappa_result[0], "SigmaKappaT": tkappa_result[1],
        "basis_angle":basis_angle,
        "rhoH": rhoH,
        **td_observables
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
    parser.add_argument("--basis_angle", type=float, default=0)
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
    
                data = run_2d(basis_angle=parsed.basis_angle, setup_file=temp_file.name, output=None, background=background, Lx=Lx, Ly=Ly)

                data["T"]  = np.array([T ])
                data["Ax"] = np.array([Ax])
                data["Ay"] = np.array([Ay])
                data["Gx"] = np.array([Gx])
                data["Gy"] = np.array([Gy])
                data["mu"] = np.array([mu])
                data["Q"]  = np.array([Q ])
                print(f"Step done: {background}",flush=True)

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
