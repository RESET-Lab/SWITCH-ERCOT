from __future__ import print_function
from __future__ import division
from pyomo.environ import *
import os
import pandas as pd
import numpy as np


def define_arguments(argparser):
    group = argparser.add_argument_group(__name__)
    group.add_argument('--obj-min', default=False, action='store_true',
        dest='obj_min',
        help=("This is intended to change the objective function between minimize and maximize based on the MGA script")
    )

    # scenario information
    argparser.add_argument(
        "--mga-tech", default="", 
        dest='mga_tech',
        help="Name of the technology to be either maximized or minimized in the MGA module"
    )

def define_components(mod):

    #Input parameters
    #mod.mga_tech = Param()
    #mod.obj_min = Param(mod.TECH_TYPE, within=Boolean)
    mod.mga_slack = Param() #Provided in the script and then added to a csv input file
    mod.opt_cost = Param()  #Script pulls this number from the output of the first run and copies it to a csv input file

    #Constraints:
    def Mga_Upper_Cost_rule(m):
        rule = ( (1+m.mga_slack)*m.opt_cost >= m.SystemCost )
        return rule

    def Mga_Lower_Cost_rule(m):
        rule = (m.SystemCost >= m.opt_cost)
        return rule

    #Constraints
    mod.Mga_Lower_Cost = Constraint(rule=Mga_Lower_Cost_rule)
    mod.Mga_Upper_Cost = Constraint(rule=Mga_Upper_Cost_rule)

def define_dynamic_components(mod):

    #Need to define new objective function
    #Create using HSJ method

    mod.DispatchSum = Expression(
        mod.GEN_TPS,
        rule = lambda m, n, t: m.DispatchGen[n, t] * m.tp_weight[t]
    )

    mod.MinTech = Expression(
        rule = lambda m: sum( m.DispatchSum[n,t] for n,t in m.GEN_TPS if m.gen_tech[n] in m.options.mga_tech)
    )
    if mod.options.obj_min:
        mod.Min_Tech = Objective(
            rule=lambda m: m.MinTech,
            sense=minimize)
    else:
        mod.Min_Tech = Objective(
            rule=lambda m: m.MinTech,
            sense=maximize)

        #will need to deactivate the cost minimizing function (pyomo does not support solvers that can use multiple objectice functions)
        #Note I have seen other Pyomo models use multiple objective functions, but they do so by turning off pyomo and solving with other means (see Temoa)
    mod.Minimize_System_Cost.deactivate()

    mod.Min_Tech.activate()

def load_inputs(mod, switch_data, inputs_dir):

    #single input file with just two values, one for slack and one for number of iterations
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'mga.csv'),
        autoselect=True,
        param=(mod.mga_slack, mod.opt_cost))
