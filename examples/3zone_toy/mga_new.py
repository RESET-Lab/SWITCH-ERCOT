from __future__ import print_function
from __future__ import division
from pyomo.environ import *
import os
import pandas as pd
import numpy as np

#Currently incorporating this into the model using a script which adds it in later
#Maybe a post-solve function would be a way to bake this into SWITCH directly? Temoa, EXPANSE, and PyPSA may be useful
    #PyPSA does a scripting approach using Snakemake, which they use to manage all their runs I believe
    #Temoa may still be helpful, but the only instance I found of them incorporating the iterations into the code directly was when they were not using
        #the Pyomo solver. Unclear if that has anything to do with this, but I want to understand how they are solving without Pyomo before I try and copy that
        #However, if we need to devise a way to turn off pyomo then that may be beyond the scope of this current attempt
    #EXPANSE, still have not looked yet

def define_components(m):

    #Input parameters
    m.TECH_TYPE = Set()
    m.mga_slack = Param() #Provided in the script and then added to a csv input file
    m.opt_cost = Param()  #Script pulls this number from the output of the first run and copies it to a csv input file
    m.prev_gen_tech = Param(m.TECH_TYPE)
    m.old_coefficients = Param(m.TECH_TYPE)
    m.hsj_coefficients = Var(m.TECH_TYPE)#, mutable=True)
    #m.dispatch_energy = Param(m.TECH_TYPE)

    #Constraints:
    def Mga_Upper_Cost_rule(m):
        rule = ((1+m.mga_slack)*m.opt_cost >= m.SystemCost)
        return rule

    def Mga_Lower_Cost_rule(m):
        rule = (m.SystemCost >= m.opt_cost)
        return rule

    #Constraints, shouldnt need to index since this is only being applied to total cost
    m.Mga_Lower_Cost = Constraint(rule=Mga_Lower_Cost_rule)
    m.Mga_Upper_Cost = Constraint(rule=Mga_Upper_Cost_rule)

    def Update_Coeff_rule(m,g):
        method = "integer"
        #Checking the dispatch of each tech type from previous iteration, adding to the coefficient
        if method == "integer":
            if m.prev_gen_tech[g] >= 1e-9:
                return m.hsj_coefficients[g] == m.old_coefficients[g] + 1
            else:
                return m.hsj_coefficients[g] == m.old_coefficients[g]
        else: #this is for the normalized method
            total_dispatch = sum(m.prev_gen_tech[n] for n in m.TECH_TYPE)
            return m.hsj_coefficients[g] == m.old_coefficients[g] + m.prev_gen_tech[g] / total_dispatch

    m.Update_Coeff = Constraint(m.TECH_TYPE, rule=Update_Coeff_rule)


def define_dynamic_components(m):

    #Need to define new objective function
    #Create using HSJ method

    def hop_skip_jump(m, g):

        #Future idea: Update coefficients in post_solve based on the results of that run and then export
        #If we do that we probably dont need to export the breakdown by technology anymore, which would save a lot of headache
        #Wait but actually that info is helpful to have as a reference, but you could also get that manually from output if you really wanted it
        dispatch_normalized_dat = [{
            "generation_project": n,
            "TECH_TYPE": m.gen_tech[n],
            #"period": m.tp_period[t],
            "dispatch_tech": (
                #m.DispatchGen[n, t] * m.tp_weight_in_year[t] / 1000),
                m.DispatchGen[n, t] * m.tp_weight[t] / 1000)
                #m.DispatchGen[n, t] * m.tp_duration_hrs[t] / 1000)

        } for n, t in m.GEN_TPS ]
        dispatch_full_df = pd.DataFrame(dispatch_normalized_dat)
        annual_dispatch = dispatch_full_df.groupby(
            ['TECH_TYPE']).sum()

        #print(annual_dispatch)
        dispatch_by_tech = annual_dispatch['dispatch_tech']
        return m.hsj_coefficients[g]*dispatch_by_tech[g]
        #return sum(m.hsj_coefficients[g]*dispatch_by_tech[g] for g in m.TECH_TYPE)

    #m.DispatchPerTech = Expression(
    #    m.TECH_TYPE,
    #    rule = hop_skip_jump)
    m.DispatchPerTech = Expression(
        m.TECH_TYPE,
        rule = hop_skip_jump)

    m.HopSkipJump = Expression(
        rule = lambda m : sum(m.DispatchPerTech[g] for g in m.TECH_TYPE)
    )

    m.Min_HSJ = Objective(
        rule=lambda m: m.HopSkipJump,
        sense=minimize)

    #will need to deactivate the cost minimizing function (pyomo does not support solvers that can use multiple objectice functions)
    #Note I have seen other Pyomo models use multiple objective functions, but they do so by turning off pyomo and solving with other means (see Temoa)
    m.Minimize_System_Cost.deactivate()

    m.Min_HSJ.activate()
    """
    m.HopSkipJump = Expression(
        rule = hop_skip_jump
    )

    m.Minimize_System_Cost.deactivate()

    m.Min_HSJ = Objective(
        rule=lambda m: m.HopSkipJump,
        sense=minimize)
    """
    #will need to deactivate the cost minimizing function (pyomo does not support solvers that can use multiple objectice functions)
    #Note I have seen other Pyomo models use multiple objective functions, but they do so by turning off pyomo and solving with other means (see Temoa)


    #m.Min_HSJ.activate()




def load_inputs(m, switch_data, inputs_dir):

    #single input file with just two values, one for slack and one for number of iterations
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'mga.csv'),
        autoselect=True,
        param=(m.mga_slack, m.opt_cost))

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'prev_run_gen.csv'),
        autoselect=True,
        index=m.TECH_TYPE,
        param=(m.prev_gen_tech, m.old_coefficients))

    #print(m.prev_gen_tech)

def post_solve(instance, outdir):
    m = instance

    dispatch_normalized_dat = [{
        "generation_project": g,
        "TECH_TYPE": instance.gen_tech[g],
        "period": instance.tp_period[t],
        "prev_gen_tech": value(
            m.DispatchGen[g, t] * m.tp_weight[t] / 1000),
    } for g, t in m.GEN_TPS ]
    dispatch_full_df = pd.DataFrame(dispatch_normalized_dat)

    coefficients = [{
        'TECH_TYPE' : g,
        "hsj_coefficients" : value(m.hsj_coefficients[g])}
        for g in m.TECH_TYPE]
    coefficients_df = pd.DataFrame(coefficients)
    coefficients_df = coefficients_df.set_index("TECH_TYPE")

    #print(coefficients_df)
    #test_1 = dispatch_full_df.groupby(
    #    ['gen_tech', "gen_load_zone", "gen_energy_source", "period"]
    #).sum()

    annual_dispatch = dispatch_full_df.groupby(
        ['TECH_TYPE']
    ).sum()
    dispatch_coeff = annual_dispatch.join(coefficients_df)
    #dispatch_coeff = dispatch_coeff.rename(columns={'gen_tech' : 'TECH_TYPE', 'Energy_GWh_typical_yr' : 'prev_gen_tech'})
    dispatch_coeff.to_csv(
        os.path.join(outdir, "dispatch_coeff.csv"),
        columns=["prev_gen_tech", "hsj_coefficients"]
    )

    #input_data = [{
    #    'TECH_TYPE' : g,
    #    'prev_gen_tech' : m.prev_gen_tech[g],
    #    'coefficients' : m.hsj_coefficients[g]
    #} for g in m.TECH_TYPE]

    #validation_output = pd.DataFrame(input_data)
    #validation_output=validation_output.set_index("TECH_TYPE")
    #print(validation_output)
    #validation_output.to_csv(
    #    os.path.join(outdir, "confirm_inputs.csv"),
    #)
