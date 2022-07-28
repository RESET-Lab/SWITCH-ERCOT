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
#Is the HSJ method going to be the best going forward? Hard to say, but there are extensions of this work that try to develop a more comprehensive approach to
#defining the near-optimal solution space. Maybe they are worth looking into more in the future
def define_components(mod):

    #Input parameters
    mod.TECH_TYPE = Set(dimen=1)
    mod.mga_slack = Param() #Provided in the script and then added to a csv input file
    mod.opt_cost = Param()  #Script pulls this number from the output of the first run and copies it to a csv input file
    mod.prev_gen_tech = Param(mod.TECH_TYPE)
    #mod.old_coefficients = Param(mod.TECH_TYPE)
    #mod.hsj_coefficients = Var(mod.TECH_TYPE, initialize=0)
    mod.hsj_coefficients = Param(mod.TECH_TYPE, mutable=True)
    #mod.dispatch_energy = Param(mod.TECH_TYPE)

    #Constraints:
    def Mga_Upper_Cost_rule(m):
        rule = ( (1+m.mga_slack)*m.opt_cost >= m.SystemCost )
        return rule

    def Mga_Lower_Cost_rule(m):
        rule = (m.SystemCost >= m.opt_cost)
        return rule

    #Constraints, shouldnt need to index since this is only being applied to total cost
    mod.Mga_Lower_Cost = Constraint(rule=Mga_Lower_Cost_rule)
    mod.Mga_Upper_Cost = Constraint(rule=Mga_Upper_Cost_rule)
    """
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

    mod.Update_Coeff = Constraint(mod.TECH_TYPE, rule=Update_Coeff_rule)
    """
def define_dynamic_components(mod):

        #Need to define new objective function
        #Create using HSJ method

    def hop_skip_jump(m, g):

        method = "normalized"
        #Checking the dispatch of each tech type from previous iteration, adding to the coefficient
        if method == "integer":
            if m.prev_gen_tech[g] >= 1e-9:
                m.hsj_coefficients[g] += 1
        else: #this is for the normalized method
            total_dispatch = sum(m.prev_gen_tech[n] for n in m.TECH_TYPE)
            m.hsj_coefficients[g] += m.prev_gen_tech[g] / total_dispatch

        #Future idea: Update coefficients in post_solve based on the results of that run and then export
        #If we do that we probably dont need to export the breakdown by technology anymore, which would save a lot of headache
        #Wait but actually that info is helpful to have as a reference, but you could also get that manually from output if you really wanted it

        dispatch_normalized_dat = [{
            "generation_project": n,
            "TECH_TYPE": m.gen_tech[n],
            #"period": m.tp_period[t],
            "gen_tech": (
                #m.DispatchGen[n, t] * m.tp_weight_in_year[t] / 1000),
                m.DispatchGen[n, t] * m.tp_weight[t] / 1000),

        } for n, t in m.GEN_TPS ]
        dispatch_full_df = pd.DataFrame(dispatch_normalized_dat)
        annual_dispatch = dispatch_full_df.groupby(
            ['TECH_TYPE']
        ).sum()
        #print(annual_dispatch)
        m.dispatch_by_tech = annual_dispatch['gen_tech']

        #rule = 1

        return m.hsj_coefficients[g]*m.dispatch_by_tech[g]

    mod.DispatchPerTech = Expression(
        mod.TECH_TYPE,
        rule = hop_skip_jump)

    mod.HopSkipJump = Expression(
        rule = lambda m : sum(m.DispatchPerTech[g] for g in m.TECH_TYPE)
        )

    mod.Min_HSJ = Objective(
        rule=lambda m: m.HopSkipJump,
        sense=minimize)

        #will need to deactivate the cost minimizing function (pyomo does not support solvers that can use multiple objectice functions)
        #Note I have seen other Pyomo models use multiple objective functions, but they do so by turning off pyomo and solving with other means (see Temoa)
    mod.Minimize_System_Cost.deactivate()

    mod.Min_HSJ.activate()

def load_inputs(mod, switch_data, inputs_dir):

    #single input file with just two values, one for slack and one for number of iterations
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'mga.csv'),
        autoselect=True,
        param=(mod.mga_slack, mod.opt_cost))

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'prev_run_gen.csv'),
        autoselect=True,
        index=mod.TECH_TYPE,
        param=(mod.prev_gen_tech, mod.hsj_coefficients))

    #print(m.prev_gen_tech)

def post_solve(instance, outdir):
    #m = instance
    '''
    Okay here is what I need to do:
    1. Find the total dispatch (in GWh) of each technology at each load zone for every timepoint (done)
    2. Subtract the CCS energy load from all technologies where that is relevant (done)
        2a. Need to add in the load zones and timepoints, will be important for indexing (done)
    3. Calculate the losses along each transmission line at each timepoint (done)
    4. Assign those losses to the load zone that is the origin of the TX dispatch decision (done)
    5. Calculate the proportion of each technology's dispatch to the total dispatch (done)
    6. Distribute TX Losses across each technology based on those proportions (done)
    7. Sum total dispatch of each technology across load zones and timepoints (done)

    '''

    test = [{
        "generation_project": g,
        "timepoint": t,
        "load_zone": instance.gen_load_zone[g],
        "TECH_TYPE": instance.gen_tech[g],
        "period": instance.tp_period[t],
        "total_dispatch_GWh": value(
            instance.DispatchGen[g, t] * instance.tp_weight[t] / 1000),
        "ccs_value": (instance.gen_ccs_energy_load[g] \
            if (hasattr(instance, 'gen_ccs_energy_load') and g in instance.CCS_EQUIPPED_GENS) else 0),
        "dispatch_after_ccs": value(
            instance.DispatchGen[g, t] * instance.tp_weight[t] / 1000)
    } for g, t in instance.GEN_TPS ]
    test_df = pd.DataFrame(test)
    test_df["ccs_load"] = test_df["total_dispatch_GWh"]*test_df["ccs_value"]
    test_df.loc[test_df["total_dispatch_GWh"]>1e-9, "dispatch_after_ccs"] = (test_df["total_dispatch_GWh"]-test_df["ccs_load"])

    total_dispatch = test_df.groupby(["load_zone", "timepoint"]).sum()
    total_dispatch = total_dispatch.rename(columns={"dispatch_after_ccs":"zonal_dispatch"})
    total_dispatch = total_dispatch.drop(["period", "total_dispatch_GWh", "ccs_load"], axis=1)
    test_df = test_df.merge(total_dispatch, on=["load_zone", "timepoint"])
    test_df["proportion"] = test_df["dispatch_after_ccs"]/test_df["zonal_dispatch"]
    #test_df.to_csv(
    #    os.path.join(outdir, "test.csv"))

    tx_test = [{
        "line": tx,
        "lz1": instance.trans_lz1[tx],
        "lz2": instance.trans_lz2[tx],
        "efficiency": instance.trans_efficiency[tx]
    } for tx in instance.TRANSMISSION_LINES]
    tx_df = pd.DataFrame(tx_test)
    #tx_df.to_csv(
    #    os.path.join(outdir, "txtest.csv"))

    tx_disp = [{
        "lz1": m,
        "lz2": n,
        "timepoint":t,
        "tp_weight": instance.tp_weight[t],
        "dispatch": value(instance.DispatchTx[m,n,t]),
        "total_dispatch": value(instance.DispatchTx[m,n,t])*instance.tp_weight[t]/1000
    } for m,n,t in instance.TRANS_TIMEPOINTS]
    tx_disp_test = pd.DataFrame(tx_disp)


    trans1 = tx_disp_test.merge(tx_df, left_on=["lz1", "lz2"], right_on=["lz1", "lz2"])
    trans2 = tx_disp_test.merge(tx_df, left_on=["lz2", "lz1"], right_on=["lz1", "lz2"])
    trans2 = trans2.rename(columns={"lz1_x":"lz1", "lz2_x":"lz2"})
    trans_total=trans1.append(trans2)
    trans_total = trans_total.drop(["lz1_y","lz2_y"], axis=1)
    trans_total["tx_losses"] = trans_total["total_dispatch"]*(1-trans_total["efficiency"])

    trans_zonal_losses = trans_total.loc[:,["lz1", "timepoint", "tx_losses"]].groupby(["lz1", "timepoint"]).sum()

    mga_dataset = test_df.merge(trans_zonal_losses, left_on=["load_zone", "timepoint"], right_on=["lz1", "timepoint"])
    mga_dataset["prev_gen_tech"] = (mga_dataset["dispatch_after_ccs"]-mga_dataset["proportion"]*mga_dataset["tx_losses"])*(1-instance.distribution_loss_rate)
    mga_dataset.to_csv(
        os.path.join(outdir, "mga_dataset.csv"))

    """
    dispatch_normalized_dat = [{
        "generation_project": g,
        "TECH_TYPE": instance.gen_tech[g],
        "period": instance.tp_period[t],
        "prev_gen_tech": value(
            instance.DispatchGen[g, t] * instance.tp_weight[t] / 1000),
    } for g, t in instance.GEN_TPS ]
    dispatch_full_df = pd.DataFrame(dispatch_normalized_dat)
    """

    coefficients = [{
        'TECH_TYPE' : g,
        "hsj_coefficients" : value(instance.hsj_coefficients[g])}
        for g in instance.TECH_TYPE]
    coefficients_df = pd.DataFrame(coefficients)
    coefficients_df = coefficients_df.set_index("TECH_TYPE")

    #print(coefficients_df)
    #test_1 = dispatch_full_df.groupby(
    #    ['gen_tech', "gen_load_zone", "gen_energy_source", "period"]
    #).sum()

    annual_dispatch = mga_dataset.loc[:,["TECH_TYPE","prev_gen_tech"]].groupby(
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
