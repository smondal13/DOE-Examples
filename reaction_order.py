import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar, Simulator
from pyomo.contrib.parmest.experiment import Experiment

class ReactionOrder(Experiment):
    """
    This example data is imported Hill chapter 3 reaction engineering
    Example 3.1 p.37 "ILLUSTRATION 3.1 Use of a Differential Method to Determine a
    Pseudo Reaction Rate Expression for the Iodine-Catalyzed Bromination of m-Xylene"
    """
    """
    Arguments
    ---------
    data: object containing vital experimental information
    nfe: number of finite elements
    ncp: number of collocation points for the finite elements
    """

    def __init__(self, data, nfe, ncp):
        self.data = data
        self.nfe = nfe
        self.ncp = ncp
        self.model = None

        # End of constructor definition
        ########################################################################

    def get_labeled_model(self):
        if self.model is None:
            self.create_model()
            self.finalize_model()
            self.label_experiment()
            return self.model

    ########################################################################
    # Create flexible model without data
    def create_model(self):
        m = self.model = pyo.ConcreteModel()

        ########################################################################
        # Define model variables
        # time
        m.t = ContinuousSet(bounds=(0, 63))  # the documentation says bounds as tuple

        # concentration of A
        m.CA = pyo.Var(m.t, domain=pyo.NonNegativeReals)  # within is an alias of domain


        # Parameters
        """
        rate constant and order of reaction (they are parameters, 
        but we will treat them as Var since we seek their values)
        """
        m.k = pyo.Var(domain=pyo.PositiveReals)  # k is the rate constant. I do not want the rate to be 0.
        # that's why I used PositiveReals instead of NonNegativeReals
        m.n = pyo.Var(domain=pyo.NonNegativeReals)  # n is the order of the rxn

        # Differential variable
        m.dCAdt = DerivativeVar(m.CA, wrt=m.t)

        # End of variable def
        ########################################################################

        ########################################################################
        # Equation Definition
        # Expression for rate
        @m.Constraint(m.t)
        def CA_rxn_rate(m, t):
            return m.dCAdt[t] == -m.k * m.CA[t] ** m.n
            # In the example k was k[t], since there was T[t] in the expression of k.
            # I don't have T[t] here, since it's isothermal, and so we will have a fixed k and n.
            # Therefore, I only have m.k and m.n

        # End of Expression definition
        ########################################################################

    def finalize_model(self):
        m = self.model
        ## I have not set the control points
        ##########################################################################
        ###################################
        ########################
        ###################################
        ##########################################################################

        # Set initial concentration value
        m.CA[0].value = self.data["CA"][0]  # do I need a separate CA0 in the data?
        # or can we just use the first element of "CA" in the data?

        # update model time `t` with time range
        m.t.update(self.data["t_range"])
        # may need to add control points here????????   #136

        # Fix the unknown parameter values
        m.k.fix(self.data["k"])  # fixing the rate constant
        m.n.fix(self.data["n"])  # fixing the reaction order

        # Add lower bound of CA[0]
        m.CA[0].setlb(0.1)
        # need to set upper bound???????????????
        ### Control points again #line148

        # Discretizing the model
        discr = pyo.TransformationFactory("dae.collocation")
        discr.apply_to(m, nfe=self.nfe, ncp=self.ncp, wrt=m.t)

        # End of model finalization
        ####################################################

    ####################################################
    # Labeling the experiment
    def label_experiment(self):
        m = self.model

        # Set measurement labels
        m.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        # Add CA to experiment outputs
        m.experiment_outputs.update((m.CA[t], None) for t in self.data["t"])

        # Add measurement error. We will assume a constant error of 0.01 M
        m.measurement_error = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        conc_measure_error = 0.01  # 0.01 M
        m.measurement_error.update((m.CA[t], conc_measure_error) for t in self.data["t"])

        # Identify design variables
        m.experiment_inputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        # Addd experimental input label for initial concentration
        m.experiment_inputs[m.CA[m.t.first()]] = None

        # Add unknown parameter labels
        m.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.unknown_parameters.update((p, pyo.value(p)) for p in [m.k, m.n])

        # End of model labeling
        ####################################################


