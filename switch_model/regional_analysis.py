from __future__ import print_function
from __future__ import division
from pyomo.environ import *
import os
import pandas as pd
import numpy as np

def define_components(mod):

    #Input parameters
    mod.TECHS = Set()
    mod.system_min = Param(mod.TECHS)
    mod.health_multiplier = Param(mod.TECHS) 
    mod.equity_index = Param(mod.LOAD_ZONES)
    mod.health_index = Param(mod.LOAD_ZONES)
    mod.cost_scal = Param()  #Script pulls this number from the output of the first run and copies it to a csv input file
    mod.equity_scal = Param() #Provided in the script and then added to a csv input file
    mod.health_scal = Param() 
    mod.alpha = Param()

def define_dynamic_components(mod):

    # Constraining capacity so ensure that any amount built in the system wide model is also built here
    def MinTech_rule(m, g):
        total = sum( m.BuildGen[n, p] for n,p in m.GEN_BLD_YRS if m.gen_tech[n] == g )
        return total >= m.system_min[g]


    mod.MinTech = Constraint( mod.TECHS,
        rule = MinTech_rule #lambda m: sum( m.BuildGen[n, p] for n,p in m.GEN_BLD_YRS if m.gen_tech[n] in m.mga_tech )
    )

    def MaxTech_rule(m, g):
        total = sum( m.BuildGen[n, p] for n,p in m.GEN_BLD_YRS if m.gen_tech[n] == g )
        return total <= 50*m.system_min[g]


    mod.MaxTech = Constraint( mod.TECHS,
        rule = MaxTech_rule #lambda m: sum( m.BuildGen[n, p] for n,p in m.GEN_BLD_YRS if m.gen_tech[n] in m.mga_tech )
    )

    #Do some kind of equity calculation
    mod.SystemEquity = Expression(
        rule = lambda m: sum( m.BuildGen[n, p]*m.equity_index[m.gen_load_zone[n]]**2 for n,p in m.GEN_BLD_YRS )# / sum( m.BuildGen[n, p] for n,p in m.GEN_BLD_YRS )
    )

    mod.SystemHealth = Expression(
        rule = lambda m: sum( m.BuildGen[n, p]*m.health_index[m.gen_load_zone[n]]**2*m.health_multiplier[m.gen_tech[n]] for n,p in m.GEN_BLD_YRS )# / sum( m.BuildGen[n, p] for n,p in m.GEN_BLD_YRS )
    )

    mod.MultiObj = Objective(
        rule=lambda m: m.alpha*m.cost_scal*m.SystemCost + (1-m.alpha)*(m.equity_scal*m.SystemEquity + m.health_scal*m.SystemHealth)
    )
        #will need to deactivate the cost minimizing function (pyomo does not support solvers that can use multiple objective functions)
        #Note I have seen other Pyomo models use multiple objective functions, but they do so by turning off pyomo and solving with other means (see Temoa)
    mod.Minimize_System_Cost.deactivate()

    mod.MultiObj.activate()


def load_inputs(mod, switch_data, inputs_dir):

    #single input file with just two values, one for slack and one for number of iterations
    
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'scaling.csv'),
        optional=True,
        autoselect=True,
        param=(mod.cost_scal, mod.equity_scal, mod.health_scal, mod.alpha))
    
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'regional_tech.csv'),
        optional=True,
        autoselect=True,
        index=mod.TECHS,
        param=(mod.system_min, mod.health_multiplier))
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'equity_data.csv'),
        optional=True,
        autoselect=True,
        index=mod.LOAD_ZONES,
        param=(mod.equity_index, mod.health_index))

def post_solve(instance, outdir):
    """
    Export storage build information to storage_builds.csv, and storage
    dispatch info to storage_dispatch.csv

    """

    import switch_model.reporting as reporting

    #Misc output table I can use for debugging
    reporting.write_table(
        instance,
        output_file=os.path.join(outdir, "metrics.csv"),
        headings=("SystemCost", 'SystemEquity', "SystemHealth"),
        values=lambda m: (
            m.SystemCost, m.SystemEquity, m.SystemHealth
        ))

    #instance.OverBuild_Penalty.display()
