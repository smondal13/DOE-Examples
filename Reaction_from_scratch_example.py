# continuation from reaction_order.py
from  pyomo.common.dependencies import numpy as np
from reaction_order import ReactionOrder
from pyomo.contrib.doe import DesignOfExperiments
import pyomo.environ as pyo

# the data for DOE
data = {
    "t" : [0, 2.25, 4.50],
    "CA" : [0.3335, 0.2965, 0.2660],
    "k" : 0.12,
    "n" : 1.6,
    "t_range" : [0, 63]
}
# actual k = 0.1, n = 1.54

def run_reactor_doe():
    data_ex = data

    # Create a ReactorOrder object
    experiment = ReactionOrder(data = data_ex, nfe = 10, ncp = 3)

    # Use the determinant (D-optimality) with scaled sensitivity matrix
    objective_option = "determinant"
    scale_nominal_param_value = True

    # Create the DesignOfExperiment object
    doe_obj = DesignOfExperiments(
        experiment,
        objective_option= objective_option,
        scale_nominal_param_value=scale_nominal_param_value
    )

    # Make design ranges to compute the full factorial design
    design_ranges = {"CA[0]" : [0.1, 5, 10]}

    # Compute the full factorial design  with the sequential FIM calculation
    doe_obj.compute_FIM_full_factorial(design_ranges = design_ranges, method = "sequential")

    # plot the results
    doe_obj.draw_factorial_figure(
        sensitivity_design_variables= ["CA[0]"],
        title_text = "Reaction Order and Rate Constant",
        xlabel_text= "CA0",
        log_scale= False,
        figure_file_name= "Reaction_order_example"
    )

    # End of Sensitivity Analysis
    #################################

    #################################
    # Begin optimal DoE
    doe_obj.run_doe()

    # Print out a results summary
    print("OPTIMAL EXPERIMENT VALUES: ")
    print(
        f"\tInitial Concentration: {doe_obj.results['Experiment Design'][0]: 0.2f}"
    )

    if __name__ == "__main__":
        run_reactor_doe()
