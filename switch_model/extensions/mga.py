from __future__ import print_function
from __future__ import division
from pyomo.environ import *
import os
import pandas as pd
import numpy as np

"""

"""

def define_arguments(argparser):
    group = argparser.add_argument_group(__name__)
    group.add_argument('--obj-min', default=False, action='store_true',
        dest='obj_min',
        help=("Change the objective function between minimize and maximize based on the MGA script")
    )

def define_components(mod):

    #Input parameters
    mod.TECHS = Set()
    mod.is_storage = Param(mod.TECHS, within=Boolean)
    mod.is_mga = Param(mod.TECHS, within=Boolean)
    mod.mga_tech = Set(initialize = mod.TECHS, filter = lambda m, g: m.is_mga[g])
    mod.SystemCostNPVPeriod = Param(mod.PERIODS)
    mod.mga_slack = Param() #Provided in the script and then added to a csv input file
    mod.opt_cost = Param()  #Script pulls this number from the output of the first run and copies it to a csv input file
    mod.penalty_multiplier = Param()


def define_dynamic_components(mod):

    # Define constraints based on minimum system costs in each investment period
    def Mga_Cost_rule(m, p):
        slack_bound = (1+m.mga_slack)*m.SystemCostNPVPeriod[p]
        rule = inequality(m.SystemCostNPVPeriod[p], m.SystemCostPerPeriod[p], slack_bound)
        return rule

    mod.Mga_Cost = Constraint(mod.PERIODS, rule=Mga_Cost_rule)


    def MinTech_calc(m):
        total = 0
        for g in m.mga_tech:
            if m.is_storage[g]:
                total += sum( m.BuildGen[n, p] + m.BuildStorageEnergy[n,p] for n,p in m.GEN_BLD_YRS if m.gen_tech[n] == g )
            else:
                total += sum( m.BuildGen[n, p] for n,p in m.GEN_BLD_YRS if m.gen_tech[n] == g )
        return total


    mod.MinTech = Expression(
        rule = MinTech_calc 
    )
    
    # Maximize or minimize based on user input
    if mod.options.obj_min:
        mod.Min_Tech = Objective(
            rule=lambda m: m.MinTech,
            sense=minimize)
    else:
        mod.Min_Tech = Objective(
            rule=lambda m: m.MinTech,
            sense=maximize
        )

    #will need to deactivate the cost minimizing function (pyomo does not support solvers that can use multiple objective functions)
    #Note I have seen other Pyomo models use multiple objective functions, but they do so by turning off pyomo and solving with other means (see Temoa)
    mod.Minimize_System_Cost.deactivate()

    mod.Min_Tech.activate()


def load_inputs(mod, switch_data, inputs_dir):

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'mga.csv'),
        optional=True,
        autoselect=True,
        optional_params=['penalty_multiplier'],
        param=(mod.mga_slack, mod.opt_cost, mod.penalty_multiplier))
    
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'mga_costs.csv'),
        optional=True,
        autoselect=True,
        index=mod.PERIODS,
        param=(mod.SystemCostNPVPeriod))
    
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'mga_tech.csv'),
        optional=True,
        autoselect=True,
        index=mod.TECHS,
        param=(mod.is_storage, mod.is_mga))

def post_solve(instance, outdir):
    """
    Export storage build information to storage_builds.csv, and storage
    dispatch info to storage_dispatch.csv

    """

    import switch_model.reporting as reporting

    #Misc output table I can use for debugging
    #reporting.write_table(
    #    instance, instance.mga_tech, instance.TIMEPOINTS,
    #    output_file=os.path.join(outdir, "excess_cap.csv"),
    #    headings=("gen", 'tp', "required_cap", 'capacity_in_tp', 'excess_capacity', "max"),
    #    values=lambda m, g, tp: (
    #        g, tp, m.RequiredCapacityByTech[g, tp], m.test[g,tp], m.ExcessCapacity[g,tp], m.MaxExcessCapacity
    #    ))
