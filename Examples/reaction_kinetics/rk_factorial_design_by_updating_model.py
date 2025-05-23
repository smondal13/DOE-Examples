import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.contrib.parmest.experiment import Experiment
import random
import idaes.core.solvers.get_solver
from pyomo.common.dependencies import numpy as np
from pyomo.contrib.doe import DesignOfExperiments
import pyomo.environ as pyo

# Import module
from reaction_kinetics_cls import ReactionOrder

# Data as a dictionary
bromine_data = {
    "Time": [
        0.00,
        2.25,
        4.50,
        6.33,
        8.00,
        10.25,
        12.00,
        13.50,
        15.60,
        17.85,
        19.60,
        27.00,
        30.00,
        38.00,
        41.00,
        45.00,
        47.00,
        57.00,
        63.00,
    ],
    "CA": [
        0.3335,
        0.2965,
        0.2660,
        0.2450,
        0.2255,
        0.2050,
        0.1910,
        0.1794,
        0.1632,
        0.1500,
        0.1429,
        0.1160,
        0.1053,
        0.0830,
        0.0767,
        0.0705,
        0.0678,
        0.0553,
        0.0482,
    ],
}

FIM_list = []


def run_reactor_doe():
    # Create a ReactorOrder object
    experiment = ReactionOrder(data=bromine_data, nfe=30, ncp=3)

    # Use the determinant (D-optimality) with scaled sensitivity matrix

    # solver = pyo.SolverFactory("ipopt")
    # solver.options["linear_solver"] = "MUMPS"

    # Create the DesignOfExperiments object
    # Factorial design
    CA0_values = np.linspace(0.01, 0.5, 10)
    for conc in CA0_values:
        experiment.data["CA"][0] = conc
        if experiment.model is not None:
            print("updating model ")
            experiment.update_model(CA0=conc)
        # experiment.CA0 = conc
        doe_obj = DesignOfExperiments(
            experiment,
            scale_nominal_param_value=True,
            # solver=solver,
            tee=True,
        )
        FIM = doe_obj.compute_FIM()
        print(FIM)
        FIM_list.append(FIM)


run_reactor_doe()
