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
        help=("Change the objective function between minimize and maximize based on the MGA script")
    )

    # scenario information
    #argparser.add_argument(
    #    "--mga-tech", default="", 
    #    dest='mga_tech',
    #    help="Name of the technology to be either maximized or minimized in the MGA module"
    #)

def define_components(mod):

    #Input parameters
    mod.TECHS = Set()
    mod.is_storage = Param(mod.TECHS, within=Boolean)
    mod.is_mga = Param(mod.TECHS, within=Boolean)
    mod.mga_tech = Set(initialize = mod.TECHS, filter = lambda m, g: m.is_mga[g])
    mod.SystemCostNPVPeriod = Param(mod.PERIODS)
    mod.SystemCostRealPeriod = Param(mod.PERIODS)
    #mod.obj_min = Param(mod.TECH_TYPE, within=Boolean)
    mod.mga_slack = Param() #Provided in the script and then added to a csv input file
    mod.opt_cost = Param()  #Script pulls this number from the output of the first run and copies it to a csv input file
    mod.penalty_multiplier = Param()
    #mod.penalty_cost = Param(default=1000000000)

    mod.MGA_TPS = Set(dimen=2, initialize = lambda m: ((g, t) for g in m.mga_tech for t in m.TIMEPOINTS))
    #mod.MGA_ZONES = Set(dimen=2, initialize = lambda m: ((g, t) for g in m.mga_tech for t in m.TIMEPOINTS))


def define_dynamic_components(mod):

    # I need to define a function that will essentially penalize any unused capacity that the MGA module generates in that final period
    # calculate capacity required to achieve peak generation
    # define tuning parameter, adjust as needed (start with 1.1?)
    # Change the sign of the tuning parameter depending on if we are minimizing or maximizing (negative if maximizing) - just kidding its not really a problem for minimizing 
    # so just dont apply it in that case
    #mod.DispatchByTech = Expression(mod.gen_tech, mod.TIMEPOINTS,
    #                                rule = lambda m, gen, t: sum(m.DispatchGen[g, t] for g in m.GENERATION_PROJECTS if m.gen_tech[g] == gen)
    #)
    def Mga_Cost_rule(m, p):
        slack_bound = (1+m.mga_slack)*m.SystemCostNPVPeriod[p]
        rule = inequality(m.SystemCostNPVPeriod[p], m.SystemCostPerPeriod[p], slack_bound)
        return rule
    #)
    #def Mga_Upper_Cost_rule(m):
    #    slack_bound = (1+m.mga_slack)*sum(m.SystemCostRealPeriod[p] for p in m.PERIODS)
    #    rule = inequality(sum(m.SystemCostRealPeriod[p] for p in m.PERIODS), sum(m.SystemCostPerPeriod[p]/m.bring_annual_costs_to_base_year[p]for p in m.PERIODS), slack_bound)
     #   return rule
    #Constraints
    #mod.Mga_Upper_Cost = Constraint(rule=Mga_Upper_Cost_rule)
    mod.Mga_Cost = Constraint(mod.PERIODS, rule=Mga_Cost_rule)

    def required_capacity(m, g, t):
        cf_sum = 0
        length = 0
        #g = 'Gas_Combined_Cycle'
        for proj in m.GENERATION_PROJECTS:
            #print('enter')
            if m.gen_tech[proj] == g: 
                if proj in m.VARIABLE_GENS:
                    cf_sum += m.gen_max_capacity_factor[proj, t]
                else:
                    #print('enter')
                    cf_sum += .01
                length += 1
        average_cf = cf_sum / length
        if average_cf <= 0:
            average_cf = 1#1
        #x = sum(m.zone_demand_mw[z, t] for z in m.LOAD_ZONES) / average_cf
        #print(average_cf)
        return sum(m.zone_demand_mw[z, t] for z in m.LOAD_ZONES) / average_cf

            #if t in m.TPS_FOR_GEN[proj] and m.gen_tech[proj]==g:
                #if (g in m.VARIABLE_GENS) and inequality(0, m.gen_max_capacity_factor[g, t]):
            #    if (g in m.VARIABLE_GENS) and (0 < m.gen_max_capacity_factor[g, t]):
                    #return m.DispatchGen[g, t] / m.gen_max_capacity_factor[g, t]
                    #return sum(m.zone_demand_mw[z, t] for z in m.LOAD_ZONES) / m.gen_max_capacity_factor[g, t]
            #    elif g in m.VARIABLE_GENS:
                    #return sum(m.zone_demand_mw[z, t] for z in m.LOAD_ZONES)
                    #return m.GenCapacityInTP[g, t]#0
            #    else:
                    #return sum(m.zone_demand_mw[z, t] for z in m.LOAD_ZONES)
                    #return m.DispatchGen[g, t]
        #else: 
        #    return 0
        #m.total_cap_requirement = []
        #if (m.gen_tech[g] == gen) and (t in m.TPS_FOR_GEN[g]):
        #    return sum(m.DispatchGen[g,t] for g in )
        #    if g in m.VARIABLE_GENS:
        #        if m.gen_max_capacity_factor[g, t] > 0:
        #            m.total_cap_requirement.append(m.DispatchGen[g, t] / m.gen_max_capacity_factor[g, t])
        #    else:
        #        m.total_cap_requirement.append(m.DispatchGen[g, t])
        #return sum(m.total_cap_requirement)
    
    #mod.RequiredCapacityByTech = Expression(mod.GENERATION_TECHNOLOGIES, mod.TIMEPOINTS,
    #                rule = lambda m, gen, t: sum(required_capacity(m, g, t) for g in m.GENERATION_PROJECTS if m.gen_tech[g]==gen)
    #)
    #mod.RequiredCapacityByTech = Expression(mod.GENERATION_TECHNOLOGIES, mod.TIMEPOINTS,
    mod.RequiredCapacityByTech = Expression(mod.mga_tech, mod.TIMEPOINTS,                                        
    #mod.RequiredCapacityByTech = Expression(mod.GEN_TPS,
                    rule = required_capacity
    )
    """
    def overbuild_penalty_rule(m):
        scalar=1.1
        #if not m.options.obj_min:
        #    scalar *= -1
        max_capacity_requirement = 0
        for t in m.TIMEPOINTS:
            if inequality(value(max_capacity_requirement), value(m.RequiredCapacityByTech[m.options.mga_tech, t])) :
                max_capacity_requirement = m.RequiredCapacityByTech[m.options.mga_tech, t]
                max_t = t
            
        #peak_cap = max(m.RequiredCapacityByTech[m.options.mga_tech, t] for t in m.TIMEPOINTS)
        excess_capacity = sum(m.GenCapacityInTP[g, max_t] for g in m.GENERATION_PROJECTS if (m.gen_tech[g] in m.options.mga_tech) and (max_t in m.TPS_FOR_GEN[g])) - max_capacity_requirement
        return excess_capacity*scalar
    """
    mod.MaxExcessCapacity = Var(within=NonNegativeReals)
    mod.test = Expression(mod.MGA_TPS, 
                                    #rule = lambda m, gen, t: sum(m.GenCapacityInTP[g, t] for g in m.GENERATION_PROJECTS if m.gen_tech[g] == gen and (g,t) in m.GEN_TPS) - m.RequiredCapacityByTech[gen, t])
                                    rule = lambda m, gen, t: sum(m.GenCapacityInTP[g, t] for g in m.GENERATION_PROJECTS if (m.gen_tech[g] == gen and (g,t) in m.GEN_TPS)))
    mod.ExcessCapacity = Expression(mod.MGA_TPS, 
                                    #rule = lambda m, gen, t: sum(m.GenCapacityInTP[g, t] for g in m.GENERATION_PROJECTS if m.gen_tech[g] == gen and (g,t) in m.GEN_TPS) - m.RequiredCapacityByTech[gen, t])
                                    rule = lambda m, gen, t: sum(m.GenCapacityInTP[g, t] for g in m.GENERATION_PROJECTS if (m.gen_tech[g] == gen and (g,t) in m.GEN_TPS)) - m.RequiredCapacityByTech[gen, t])
    #mod.OverBuild_Penalty = Constraint(mod.GENERATION_TECHNOLOGIES, mod.TIMEPOINTS,
    mod.OverBuild_Penalty = Constraint(mod.MGA_TPS,                                   
    #    rule = lambda m, g, t: m.MaxExcessCapacity >= m.ExcessCapacity[g, t] 
        rule = lambda m, g, t: inequality(None, m.ExcessCapacity[g, t], m.MaxExcessCapacity)
    )

    #mod.ExcessCapacityCost = Expression(mod.TIMEPOINTS, rule = lambda m, t:
    #                                    (sum(m.GenCapacityInTP[g, t] for g in m.GENERATION_PROJECTS if (m.gen_tech[g] in m.mga_tech and (g,t) in m.GEN_TPS)) - sum(m.RequiredCapacityByTech[g, t] for g in m.mga_tech))*m.penalty_cost)
                                        #sum(m.ExcessCapacity[g,t] for g in m.mga_tech)*m.penalty_cost)
    
    #mod.Cost_Components_Per_TP.append('ExcessCapacityCost')
    #mod.Test_Penalty = Constraint(
    #    rule = lambda m: m.MaxExcessCapacity >= 10.0
    #)
    def MinTech_calc(m):
        total = 0
        for g in m.mga_tech:
            if m.is_storage[g]:
                total += sum( m.BuildGen[n, p] + m.BuildStorageEnergy[n,p] for n,p in m.GEN_BLD_YRS if m.gen_tech[n] == g )
            else:
                total += sum( m.BuildGen[n, p] for n,p in m.GEN_BLD_YRS if m.gen_tech[n] == g )
        return total


    mod.MinTech = Expression(
        rule = MinTech_calc #lambda m: sum( m.BuildGen[n, p] for n,p in m.GEN_BLD_YRS if m.gen_tech[n] in m.mga_tech )
    )
    
    if mod.options.obj_min:
        mod.Min_Tech = Objective(
            rule=lambda m: m.MinTech,
            sense=minimize)
    else:
        mod.Min_Tech = Objective(
            #rule=lambda m: m.MinTech - m.MaxExcessCapacity*m.penalty_multiplier,
            rule=lambda m: m.MinTech,
            sense=maximize
        )

        #will need to deactivate the cost minimizing function (pyomo does not support solvers that can use multiple objective functions)
        #Note I have seen other Pyomo models use multiple objective functions, but they do so by turning off pyomo and solving with other means (see Temoa)
    mod.Minimize_System_Cost.deactivate()

    mod.Min_Tech.activate()


def load_inputs(mod, switch_data, inputs_dir):

    #single input file with just two values, one for slack and one for number of iterations
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
    reporting.write_table(
        instance, instance.mga_tech, instance.TIMEPOINTS,
        output_file=os.path.join(outdir, "excess_cap.csv"),
        headings=("gen", 'tp', "required_cap", 'capacity_in_tp', 'excess_capacity', "max"),
        values=lambda m, g, tp: (
            g, tp, m.RequiredCapacityByTech[g, tp], m.test[g,tp], m.ExcessCapacity[g,tp], m.MaxExcessCapacity
        ))

    #instance.OverBuild_Penalty.display()
