#!/usr/bin/env python
# coding: utf-8

# ## Importing packages

# In[5]:


import pyomo.environ as pyo
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.contrib.parmest.experiment import Experiment


# ## Dataset for the model

# In[17]:


CA0 = 0.3335
# Data as a dictionary
bromine_data = {
    "Time": [0.00, 2.25, 4.50, 6.33, 8.00, 10.25, 12.00, 13.50, 15.60, 17.85,
                   19.60, 27.00, 30.00, 38.00, 41.00, 45.00, 47.00, 57.00, 63.00],
    "CA": [0.3335, 0.2965, 0.2660, 0.2450, 0.2255, 0.2050, 0.1910, 0.1794, 0.1632, 
           0.1500, 0.1429, 0.1160, 0.1053, 0.0830, 0.0767, 0.0705, 0.0678, 0.0553, 0.0482]
}


# ## Helper class to make the input data as a `dict`

# In[1]:


# Creating a helper class that builds a dictionary to help initialization
def helper(my_array, time):
    assert len(my_array) == len(time), "Dimension of the entries DO NOT MATCH"
    data2 = {}
    for k, t in enumerate(time):
        if my_array[k] is not None:
            data2[t] = my_array[k]
        else:
            data2[t] = 0
    return data2


# ## Class for ``ParmEst``

# In[52]:


class ReactionOrder(Experiment):
    """
    This example data is imported Hill chapter 3 reaction engineering
    Example 3.1 p.37 "ILLUSTRATION 3.1 Use of a Differential Method to Determine a
    Pseudo Reaction Rate Expression for the Iodine-Catalyzed Bromination of m-Xylene"    
    """
    
    def __init__(self, data, CA0,  theta_initial = None, nfe, ncp):
        """
        Arguments:
        data : object containing vital experimental information. dictionary, pandas dataframe
        theta_initial : initial values of the parameters theta. dtype: dict
        nfe : number of number of finite elements
        ncp : number of collocation points
        """
        self.data = data
        self.nfe = nfe  # number of finite elements
        self.ncp = ncp  # number of collocation points
        self.model = None   

        # Defining the parameter values
        if theta_initial is None:
            self.theta_initial = {
                "k" = 0.10, # rate constant (L/mol)^0.5 /min
                "n" = 1.5
            }
        else:
            self.theta_initial = theta_initial
        
    # End of constructor definition
    ########################################################################

    def get_labeled_model(self):
        # get_labeled_model is a mandatory method in both parmest and DOE.
        if self.model is None:
            self.create_model()
            self.finalize_model()
            self.label_experiment()
        return self.model

    ########################################################################
    # Create flixible model without data
    def create_model(self):
        """
        Here, we will create different variables, parameters, and constraints
        """
        #m = self.model = reaction_model(self.data)
        m = self.model = pyo.ConcreteModel()

        ########################################################################
        # Define model variables
        # time
        m.t = ContinuousSet(bounds=(0, tf))
        
        # concentration of A and B
        m.CA = pyo.Var(m.t, domain = pyo.NonNegativeReals)
        m.CB = pyo.Var(m.t, domain = pyo.NonNegativeReals)

        m.k = pyo.Var(domain = pyo.NonNegativeReals) 
        m.n = pyo.Var(domain = pyo.Reals)       
        
        #m.k = pyo.Var(initialize = theta_initial["k"]) # don't need to specify here. we can, but not necessary
        #m.n = pyo.Var(initialize = theta_initial["n"])
        
        # Differential variable
        m.dCAdt = DerivativeVar(m.CA, wrt = m.t)
        # m.dCBdt = DerivativeVar(m.CB, wrt = m.t)
        
        # Expression for rate
        @m.Constraint(m.t)
        def CA_rxn_rate(m, t):
            return m.dCAdt[t] == - m.k * (m.CA[t])**m.n
        
        
        return m        
        # End of Expression definiton
        ########################################################################

    def finalize_model(self):
        """
        Finalizing the model. Here, we will set the experimental conditions (e.g, initial conditions) and discretize the model
        It makes a solvable model.
        """
        m = self.model
                
        m.k.fix(self.theta_initial["k"])  # fixing the rate constant
        m.n.fix(self.theta_initial["n"])  # fixing the reaction order

        # Add lower bound of CA[0]
        m.CA[0].value = self.data["CA"][0]  # Setting the first value of CA as CA0
        m.CA[0].setlb(0.1)
        # Add upper bound of CA[0]
        m.CA[0].setub(10)  # let's set the upper bound to 10M
        
        # Control points again 
        m.t_control = self.data["Time"]

        # Discretizing the model
        discr = pyo.TransformationFactory("dae.collocation")
        discr.apply_to(m, nfe = self.nfe, ncp = self.ncp, wrt = m.t)
        # or nfe = len(self.data.Time) - 1 ???????????????????????

        # End of model finalization
        ####################################################

        ####################################################
        # Labeling the experiment
    def label_experiment(self):
        """
        The model is updated with outputs, inputs, errors and unknown parameters
        This makes the model labeled with full experiment
        """
        m = self.model

        # Set measurement labels
        m.experiment_outputs = pyo.Suffix(direction = pyo.Suffix.LOCAL)
        # Add CA to experiment outputs
        m.experiment_outputs.update((m.CA[t], None) for t in m.t_control)
        

        # Add measurement error. Let's assume an constant error of 3% CA. 
        # My next plan is to create a random error between 3% and 10% and see the result
        m.measurement_error = pyo.Suffix(direction = pyo.Suffix.LOCAL)
        m.measurement_error.update((m.CA[t], 0.03)] for t in self.data["Time"])

        # Identify design variables
        m.experiment_inputs = pyo.Suffix(direction = pyo.Suffix.LOCAL)
        # Addd experimental input label for initial concentration
        m.experiment_inputs(m.CA[m.t.first()]) = None

        # Add unknown parameter labels
        m.unknown_parameters = pyo.Suffix(direction = pyo.Suffix.LOCAL)
        m.unknown_parameters.update((p, pyo.value(p)) for p in [m.k, m.n])

        # End of model labeling
        ####################################################   


# ## Checking my model with `ipopt`

# In[16]:


tf = 63
CA0 = 0.3335  # initial
CB0 = 0
k = 0.1 
n = 1.54

m = pyo.ConcreteModel()

########################################################################
# Define model variables
# time
m.t = ContinuousSet(bounds=(0, tf))

# concentration of A and B
m.CA = pyo.Var(m.t, domain = pyo.NonNegativeReals)#, initialize= 1e-6)
m.CB = pyo.Var(m.t, domain = pyo.NonNegativeReals)#, initialize= 1e-6)

m.k = pyo.Param(initialize = k)
m.n = pyo.Param(initialize = n)

# Differential variable
m.dCAdt = DerivativeVar(m.CA, wrt = m.t)
# m.dCBdt = DerivativeVar(m.CB, wrt = m.t)

# Expression for rate
@m.Constraint(m.t)
def CA_rxn_rate(m, t):
    return m.dCAdt[t] == - m.k * (m.CA[t]+1e-12)**m.n

m.CA[0].fix(CA0)

discr = pyo.TransformationFactory('dae.collocation')
discr.apply_to(m, nfe = 20, ncp=10)

solver = pyo.SolverFactory("ipopt")
solver.solve(m, tee=True)

# print
ca = [pyo.value(m.CA[t]) for t in m.t]
print(ca)


# In[ ]:




