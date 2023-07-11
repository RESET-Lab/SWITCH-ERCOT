from __future__ import absolute_import
from switch_model.reporting import write_table
import pandas as pd
import os
from pyomo.environ import *

# from .util import get


def define_components(m):
    """
    Incorporate the effect of the production tax credit
    """
    m.credit_years = Set(
        #input_file="ptc_values.csv",
        dimen=2#,
        #input_optional=True,
    )
    m.ptc_inputs = Param(
        m.credit_years,
        #input_file="ptc_values.csv",
        #input_column="ptc_value",
        default=0,
        domain=NonNegativeReals#,
        #doc="Production Tax Credit (PTC) for given technology by period. Data in $/MWh",
    )
    m.itc_inputs = Param(
        m.credit_years,
        #input_file="ptc_values.csv",
        #input_column="ptc_value",
        default=0,
        domain=NonNegativeReals#,
        #doc="Production Tax Credit (PTC) for given technology by period. Data in $/MWh",
    )

    # Create a set that has build capacity constrained by year (both caps of the PTC).
    # The two caps of the PTC are that the generator must be built prior to 2035 in
    # order to recieve credit, and that all generators will no longer recieve the PTC
    # starting in 2040.
    m.ptc_eligible_yrs = Set(
        m.GENERATION_PROJECTS,
        m.PERIODS,
        ordered=False,
        initialize=lambda m, g, period: set(
            bld_yr
            for bld_yr in m.BLD_YRS_FOR_GEN_PERIOD[g, period]
            if 2025 <= bld_yr < 2035 and period < 2040 #projects built in 2030 would be able to receive the credit in 2040, I think? Check this later
        ),
    )
    # Calculate the total eligible PTC capacity per period
    m.PTC_Capacity = Expression(
        m.GENERATION_PROJECTS,
        m.PERIODS,
        rule=lambda m, g, period: sum(
            m.BuildGen[g, bld_yr] for bld_yr in m.ptc_eligible_yrs[g, period]
        ),
    )

    # Same as PTC_Capacity but per timepoint
    m.PTC_CapacityInTP = Expression(
        m.GEN_TPS, rule=lambda m, g, t: m.PTC_Capacity[g, m.tp_period[t]]
    )

    # Create PTC variable that will either return the PTC Capacity or the DispatchGen
    # whichever is minimum.
    m.PTC = Var(m.GEN_TPS, domain=NonNegativeReals)

    m.PTC_lower_bound = Constraint(
        m.GEN_TPS, rule=lambda m, g, t: m.PTC[g, t] <= m.PTC_CapacityInTP[g, t]
    )
    m.PTC_upper_bound = Constraint(
        m.GEN_TPS, rule=lambda m, g, t: m.PTC[g, t] <= m.DispatchGen[g, t]
    )

    # Calculate PTC
    m.PTC_per_tp = Expression(
        m.TIMEPOINTS,
        rule=lambda m, t: sum(
            -m.PTC[g, t] * m.ptc_inputs[m.tp_period[t], m.gen_tech[g]]
            for g in m.GENS_IN_PERIOD[m.tp_period[t]]
            if m.gen_tech[g] in set([item[1] for item in m.credit_years.data()])
            and m.tp_period[t] < 2040
        ),
    )
    m.Cost_Components_Per_TP.append("PTC_per_tp")

    # Calculate ITC
    m.ITC_per_period = Expression(
        m.PERIODS,
        rule = lambda m, p: sum(
            -m.itc_inputs[p, m.gen_tech[g]] * m.GenCapitalCosts[g,p]
            for (g, p) in m.NEW_GEN_BLD_YRS #m.GENS_IN_PERIOD[p]
            if m.gen_tech[g] in set([item[1] for item in m.credit_years.data()])
            #and g not in 
            and p < 2040
        )
    )
    m.Cost_Components_Per_Period.append('ITC_per_period')

def load_inputs(mod, switch_data, inputs_dir):

    #Input file has values for the PTC and the percent of the PTC each project will receive when built in a given period
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'tax_credits.csv'),
        autoselect=True,
        index=mod.credit_years,
        param=(mod.ptc_inputs, mod.itc_inputs))
    

    #switch_data.load_aug(
    #    filename=os.path.join(inputs_dir, 'ptc_projects.csv'),
    #    autoselect=True,
    #    index=mod.GENERATION_PROJECTS,
    #    param=(mod.ptc_eligible))

# Exported files:
#         PTC.csv - Total ptc value aggregated per period
def post_solve(m, outdir):
    """ Work in progress"""
    #df = pd.DataFrame(
    #    {
    #        "GENERATOR_ID": g,
    #        "timestamp": value(m.tp_timestamp[t]),
    #        "PTC_Capacity_MW": value(m.PTC_CapacityInTP[g, t]),
    #        "Online_Capacity_MW": value(m.GenCapacityInTP[g, t]),
    #    }
    #    for t in m.TIMEPOINTS
    #    for g in m.GENS_IN_PERIOD
    #    if m.gen_tech[g] in set([item[1] for item in m.credit_years.data()])
    #)
    #write_table(m, output_file=os.path.join(outdir, "PTC.csv"), df=df, index=False)
