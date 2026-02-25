# Copyright (c) 2015-2019 The Switch Authors. All rights reserved.
# Licensed under the Apache License, Version 2.0, which is in the LICENSE file.

"""
Defines transmission build-outs.
"""

import os, collections
from pyomo.environ import *
from switch_model.financials import capital_recovery_factor as crf
import pandas as pd

dependencies = 'switch_model.timescales', 'switch_model.balancing.load_zones',\
    'switch_model.financials'

def define_components(mod):
    #######################
    ######## Build ########
    #######################
    mod.CCS_PROJECTS = Set()
    mod.ccs_tech = Param(mod.CCS_PROJECTS)
    mod.CCS_TECHNOLOGIES = Set(initialize=lambda m:
        {m.ccs_tech[g] for g in m.CCS_PROJECTS}
    )
    mod.ccs_energy_source = Param(mod.CCS_PROJECTS,
        validate=lambda m,val, g: val in m.ENERGY_SOURCES or val == "multiple")
    mod.ccs_load_zone = Param(mod.CCS_PROJECTS, within=mod.LOAD_ZONES)
    mod.ccs_max_age = Param(mod.CCS_PROJECTS, within=PositiveIntegers)
    mod.min_data_check('CCS_PROJECTS', 'ccs_tech', 'ccs_energy_source',
        'ccs_load_zone', 'ccs_max_age')

    """Construct CCS_* indexed sets efficiently with a
    'construction dictionary' pattern: on the first call, make a single
    traversal through all generation projects to generate a complete index,
    use that for subsequent lookups, and clean up at the last call."""
    def CCS_IN_ZONE_init(m, z):
        if not hasattr(m, 'CCS_IN_ZONE_dict'):
            m.CCS_IN_ZONE_dict = {_z: [] for _z in m.LOAD_ZONES}
            for g in m.CCS_PROJECTS:
                m.CCS_IN_ZONE_dict[m.ccs_load_zone[g]].append(g)
        result = m.CCS_IN_ZONE_dict.pop(z)
        if not m.CCS_IN_ZONE_dict:
            del m.CCS_IN_ZONE_dict
        return result
    mod.CCS_IN_ZONE = Set(
        mod.LOAD_ZONES,
        initialize=CCS_IN_ZONE_init
    )

    def CCS_BY_TECHNOLOGY_init(m, t):
        if not hasattr(m, 'CCS_BY_TECH_dict'):
            m.CCS_BY_TECH_dict = {_t: [] for _t in m.CCS_TECHNOLOGIES}
            for g in m.CCS_PROJECTS:
                m.CCS_BY_TECH_dict[m.ccs_tech[g]].append(g)
        result = m.CCS_BY_TECH_dict.pop(t)
        if not m.CCS_BY_TECH_dict:
            del m.CCS_BY_TECH_dict
        return result
    mod.CCS_BY_TECHNOLOGY = Set(
        mod.CCS_TECHNOLOGIES,
        initialize=CCS_BY_TECHNOLOGY_init
    )

    mod.CAPACITY_LIMITED_CCS = Set(within=mod.CCS_PROJECTS)
    mod.ccs_capacity_limit_mtons = Param(
        mod.CAPACITY_LIMITED_CCS, within=NonNegativeReals)

    mod.ccs_uses_fuel = Param(
        mod.CCS_PROJECTS,
        initialize=lambda m, g: (
            m.ccs_energy_source[g] in m.FUELS
                or m.ccs_energy_source[g] == "multiple"))
    mod.NON_FUEL_BASED_CCS = Set(
        initialize=mod.CCS_PROJECTS,
        filter=lambda m, g: not m.ccs_uses_fuel[g])
    mod.FUEL_BASED_CCS = Set(
        initialize=mod.CCS_PROJECTS,
        filter=lambda m, g: m.ccs_uses_fuel[g])

    mod.ccs_capture_rate_mwh_per_ton = Param(
        mod.CCS_PROJECTS,
        within=NonNegativeReals)
    mod.ccs_fuel_mmbtu_per_ton = Param(
        mod.CCS_PROJECTS,
        within=NonNegativeReals)
    mod.FUELS_FOR_CCS = Set(mod.FUEL_BASED_CCS,
        initialize=lambda m, g:  m.ccs_energy_source[g])

    mod.CCS_BLD_YRS = Set(
        dimen=2,
        validate=lambda m, g, bld_yr: (
            (g, bld_yr) in m.CCS_PROJECTS * m.PERIODS))


    def ccs_build_can_operate_in_period(m, g, build_year, period):
        if build_year in m.PERIODS:
            online = m.period_start[build_year]
        else:
            online = build_year
        retirement = online + m.ccs_max_age[g]
        return (
            online <= m.period_start[period] < retirement
        )
        # This is probably more correct, but is a different behavior
        # mid_period = m.period_start[period] + 0.5 * m.period_length_years[period]
        # return online <= m.period_start[period] and mid_period <= retirement

    # The set of periods when a project built in a certain year will be online
    mod.PERIODS_FOR_CCS_BLD_YR = Set(
        mod.CCS_BLD_YRS,
        within=mod.PERIODS,
        ordered=True,
        initialize=lambda m, g, bld_yr: set(
            period for period in m.PERIODS
            if ccs_build_can_operate_in_period(m, g, bld_yr, period)))
    # The set of build years that could be online in the given period
    # for the given project.
    mod.BLD_YRS_FOR_CCS_PERIOD = Set(
        mod.CCS_PROJECTS, mod.PERIODS,
        initialize=lambda m, c, period: set(
            bld_yr for (ccs, bld_yr) in m.CCS_BLD_YRS
            if ccs == c and
               ccs_build_can_operate_in_period(m, c, bld_yr, period)))
    # The set of periods when a generator is available to run
    mod.PERIODS_FOR_CCS = Set(
        mod.CCS_PROJECTS,
        initialize=lambda m, c: [p for p in m.PERIODS if len(m.BLD_YRS_FOR_CCS_PERIOD[c, p]) > 0]
    )

    def bounds_BuildCCS(model, c, bld_yr):
        if(c in model.CAPACITY_LIMITED_GENS):
            # This does not replace Max_Build_Potential because
            # Max_Build_Potential applies across all build years.
            return (0, model.ccs_capacity_limit_mtons[c])
        else:
            return (0, None)
    mod.BuildCCS = Var(
        mod.CCS_BLD_YRS,
        within=NonNegativeReals,
        bounds=bounds_BuildCCS)

    # note: in pull request 78, commit e7f870d..., CCS_PERIODS
    # was mistakenly redefined as CCS_PROJECTS * PERIODS.
    # That didn't directly affect the objective function in the tests
    # because most code uses CCS_TPS, which was defined correctly.
    # But it did have some subtle effects on the main Hawaii model.
    # It would be good to have a test that this set is correct,
    # e.g., assertions that in the 3zone_toy model,
    # ('C-Coal_ST', 2020) in m.CCS_PERIODS and ('C-Coal_ST', 2030) not in m.CCS_PERIODS
    # and 'C-Coal_ST' in m.CCS_IN_PERIOD[2020] and 'C-Coal_ST' not in m.CCS_IN_PERIOD[2030]
    mod.CCS_PERIODS = Set(
        dimen=2,
        initialize=lambda m:
            [(c, p) for c in m.CCS_PROJECTS for p in m.PERIODS_FOR_CCS[c]])

    mod.CCSCapacity = Expression(
        mod.CCS_PROJECTS, mod.PERIODS,
        rule=lambda m, c, period: sum(
            m.BuildCCS[c, bld_yr]
            for bld_yr in m.BLD_YRS_FOR_CCS_PERIOD[c, period]))

    mod.Max_CCS_Potential = Constraint(
        mod.CAPACITY_LIMITED_CCS, mod.PERIODS,
        rule=lambda m, c, p: (
            m.ccs_capacity_limit_mtons[g] >= m.CCSCapacity[c, p]))

    # The following components enforce minimum capacity build-outs.
    # Note that this adds binary variables to the model.
    mod.ccs_min_build_capacity = Param (mod.CCS_PROJECTS,
        within=NonNegativeReals, default=0)
    mod.NEW_CCS_WITH_MIN_BUILD_YEARS = Set(
        initialize=mod.CCS_BLD_YRS,
        filter=lambda m, g, p: (
            m.ccs_min_build_capacity[g] > 0))
    mod.BuildMinCCSCap = Var(
        mod.NEW_CCS_WITH_MIN_BUILD_YEARS,
        within=Binary)
    mod.Enforce_Min_CCS_Lower = Constraint(
        mod.NEW_CCS_WITH_MIN_BUILD_YEARS,
        rule=lambda m, g, p: (
            m.BuildMinCCSCap[g, p] * m.ccs_min_build_capacity[g]
            <= m.BuildCCS[g, p]))

    # Define a constant for enforcing binary constraints on project capacity
    # The value of 100 GW should be larger than any expected build size. For
    # perspective, the world's largest electric power plant (Three Gorges Dam)
    # is 22.5 GW. I tried using 1 TW, but CBC had numerical stability problems
    # with that value and chose a suboptimal solution for the
    # discrete_and_min_build example which is installing capacity of 3-5 MW.
    mod._ccs_max_cap_for_binary_constraints = 10**5
    mod.Enforce_Min_CCS_Upper = Constraint(
        mod.NEW_CCS_WITH_MIN_BUILD_YEARS,
        rule=lambda m, g, p: (
            m.BuildCCS[g, p] <= m.BuildMinCCSCap[g, p] *
                mod._ccs_max_cap_for_binary_constraints))

    # Costs

    mod.ccs_overnight_cost = Param(
        mod.CCS_BLD_YRS,
        within=NonNegativeReals)
    mod.ccs_fixed_om = Param(
        mod.CCS_BLD_YRS,
        within=NonNegativeReals)
    mod.ccs_variable_om = Param(
        mod.CCS_BLD_YRS,
        within=NonNegativeReals)
    
    mod.min_data_check('ccs_overnight_cost', 'ccs_fixed_om', 'ccs_variable_om')

    # Derived annual costs
    mod.ccs_capital_cost_annual = Param(
        mod.CCS_BLD_YRS,
        initialize=lambda m, g, bld_yr: (
            m.ccs_overnight_cost[g, bld_yr] *
            crf(m.interest_rate, m.ccs_max_age[g])))

    mod.CCSCapitalCosts = Expression(
        mod.CCS_PROJECTS, mod.PERIODS,
        rule=lambda m, g, p: sum(
            m.BuildCCS[g, bld_yr] * m.ccs_capital_cost_annual[g, bld_yr]
            for bld_yr in m.BLD_YRS_FOR_CCS_PERIOD[g, p]))
    mod.CCSFixedOMCosts = Expression(
        mod.CCS_PROJECTS, mod.PERIODS,
        rule=lambda m, g, p: sum(
            m.BuildCCS[g, bld_yr] * m.ccs_fixed_om[g, bld_yr]
            for bld_yr in m.BLD_YRS_FOR_CCS_PERIOD[g, p]))
    # Summarize costs for the objective function. Units should be total
    # annual future costs in $base_year real dollars. The objective
    # function will convert these to base_year Net Present Value in
    # $base_year real dollars.
    mod.TotalCCSFixedCosts = Expression(
        mod.PERIODS,
        rule=lambda m, p: sum(
            m.CCSCapitalCosts[g, p] + m.CCSFixedOMCosts[g, p]
            for g in m.CCS_PROJECTS))
    mod.Cost_Components_Per_Period.append('TotalCCSFixedCosts')

    ######################
    ###### Dispatch ######
    ######################

    def period_active_ccs_rule(m, period):
        if not hasattr(m, 'period_active_ccs_dict'):
            m.period_active_ccs_dict = collections.defaultdict(set)
            for (_g, _period) in m.CCS_PERIODS:
                m.period_active_ccs_dict[_period].add(_g)
        result = m.period_active_ccs_dict.pop(period)
        if len(m.period_active_ccs_dict) == 0:
            delattr(m, 'period_active_ccs_dict')
        return result
    mod.CCS_IN_PERIOD = Set(mod.PERIODS, initialize=period_active_ccs_rule,
        doc="The set of projects active in a given period.")

    mod.TPS_FOR_CCS = Set(
        mod.CCS_PROJECTS,
        within=mod.TIMEPOINTS,
        initialize=lambda m, g: (
            tp for p in m.PERIODS_FOR_CCS[g] for tp in m.TPS_IN_PERIOD[p]
        )
    )

    def init(m, ccs, period):
        try:
            d = m._TPS_FOR_CCS_IN_PERIOD_dict
        except AttributeError:
            d = m._TPS_FOR_CCS_IN_PERIOD_dict = dict()
            for _ccs in m.CCS_PROJECTS:
                for t in m.TPS_FOR_CCS[_ccs]:
                    d.setdefault((_ccs, m.tp_period[t]), set()).add(t)
        result = d.pop((ccs, period), set())
        if not d:  # all gone, delete the attribute
            del m._TPS_FOR_CCS_IN_PERIOD_dict
        return result
    mod.TPS_FOR_CCS_IN_PERIOD = Set(
        mod.CCS_PROJECTS, mod.PERIODS,
        within=mod.TIMEPOINTS, initialize=init)

    mod.CCS_TPS = Set(
        dimen=2,
        initialize=lambda m: (
            (g, tp)
                for g in m.CCS_PROJECTS
                    for tp in m.TPS_FOR_CCS[g]))
    mod.FUEL_BASED_CCS_TPS = Set(
        dimen=2,
        initialize=lambda m: (
            (g, tp)
                for g in m.FUEL_BASED_CCS
                    for tp in m.TPS_FOR_CCS[g]))
    mod.CCS_TP_FUELS = Set(
        dimen=3,
        initialize=lambda m: (
            (g, t, f)
                for (g, t) in m.FUEL_BASED_CCS_TPS
                    for f in m.FUELS_FOR_CCS[g]))

    mod.CCSCapacityInTP = Expression(
        mod.CCS_TPS,
        rule=lambda m, g, t: m.CCSCapacity[g, m.tp_period[t]])
    mod.DispatchCCS = Var(
        mod.CCS_TPS,
        within=NonNegativeReals)

    mod.CCSFuelUse = Expression(
            mod.CCS_TPS,
            rule=lambda m, g, t: m.DispatchCCS[g, t]*m.ccs_fuel_mmbtu_per_ton[g])

    # mod.CCSFuelUseRate = Var(
    #     mod.CCS_TP_FUELS,
    #     within=NonNegativeReals,
    #     doc=("Other modules constraint this variable based on DispatchCCS and "
    #          "module-specific formulations of unit commitment and heat rates."))
    mod.CCSElectricityInTP = Expression(
        mod.CCS_TPS,
        rule=lambda m, g, t: m.DispatchCCS[g, t]*m.ccs_capture_rate_mwh_per_ton[g])
    mod.ZoneTotalCCSDemand = Expression(
        mod.LOAD_ZONES, mod.TIMEPOINTS,
        rule=lambda m, z, t: \
            sum(m.CCSElectricityInTP[p, t]
                for p in m.CCS_IN_ZONE[z]
                if (p, t) in m.CCS_TPS))
    mod.Zone_Power_Withdrawals.append('ZoneTotalCCSDemand')

    mod.AnnualCapture = Expression(mod.PERIODS,
        rule=lambda m, period: sum(
            -1*m.DispatchCCS[g, t] * m.tp_weight_in_year[t]
            for (g, t) in m.CCS_TPS
            if m.tp_period[t] == period),
        doc="The system's annual emissions, in metric tonnes of CO2 per year.")

    mod.System_Emissions.append('AnnualCapture')

    mod.CCSVariableOMCostsInTP = Expression(
        mod.TIMEPOINTS,
        rule=lambda m, t: sum(
            m.DispatchCCS[g, t] * m.ccs_variable_om[g, m.tp_period[t]] + m.CCSFuelUse[g, t]*m.fuel_cost[m.ccs_load_zone[g], m.ccs_energy_source[g], m.tp_period[t]]
            for g in m.CCS_IN_PERIOD[m.tp_period[t]]),
        doc="Summarize costs for the objective function")
    mod.Cost_Components_Per_TP.append('CCSVariableOMCostsInTP')

    ################
    # DAC Tax Credit
    ################

    mod.dac_credit_years = Set(
        dimen=2
    )

    mod.dac_credit_value = Param(
        mod.dac_credit_years,
        default=0,
        domain=NonNegativeReals
    )

    # Calculate Carbon Capture tax credit
    mod.dac_eligible_yrs = Set(
        mod.CCS_PROJECTS,
        mod.PERIODS,
        ordered=False,
        initialize=lambda m, g, period: set(
            bld_yr
            for bld_yr in m.BLD_YRS_FOR_CCS_PERIOD[g, period]
            if 2025 <= bld_yr < 2035 and (period-bld_yr) < 12 
        ),
    )
    # Calculate the total eligible PTC capacity per period
    mod.DAC_Capacity = Expression(
        mod.CCS_PROJECTS,
        mod.PERIODS,
        rule=lambda m, g, period: sum(
            m.BuildCCS[g, bld_yr] for bld_yr in m.dac_eligible_yrs[g, period]
        ),
    )

    # Same as PTC_Capacity but per timepoint
    mod.DAC_CapacityInTP = Expression(
        mod.CCS_TPS, rule=lambda m, g, t: m.DAC_Capacity[g, m.tp_period[t]]
    )

    # Create PTC variable that will either return the PTC Capacity or the DispatchGen
    # whichever is minimum.
    mod.DAC_credit = Var(mod.CCS_TPS, domain=NonNegativeReals)

    mod.DAC_credit_lower_bound = Constraint(
        mod.CCS_TPS, rule=lambda m, g, t: m.DAC_credit[g, t] <= m.DAC_CapacityInTP[g, t]
    )
    mod.DAC_credit_upper_bound = Constraint(
        mod.CCS_TPS, rule=lambda m, g, t: m.DAC_credit[g, t] <= m.DispatchCCS[g, t]
    )

    # Calculate PTC
    mod.DAC_credit_per_tp = Expression(
        mod.TIMEPOINTS,
        rule=lambda m, t: sum(
            -m.DAC_credit[g, t] * m.dac_credit_value[m.tp_period[t], m.ccs_tech[g]] 
            for g in m.CCS_IN_PERIOD[m.tp_period[t]]
            if m.ccs_tech[g] in set([item[1] for item in m.dac_credit_years.data()])
            and m.tp_period[t] < 2040
        ),
    )
    mod.Cost_Components_Per_TP.append("DAC_credit_per_tp")

def load_inputs(mod, switch_data, inputs_dir):
    """

    Import data describing project builds. The following files are
    expected in the input directory.

    ccs_project_info.csv has mandatory and optional columns. The
    operations.ccs_dispatch module will also look for additional columns in
    this file. You may drop optional columns entirely or mark blank
    values with a dot '.' for select rows for which the column does not
    apply. Mandatory columns are:
        "ccs_tech","ccs_load_zone","ccs_capture_rate_mwh_per_ton",
                         "ccs_capacity_limit_mtons","ccs_energy_source",
                         "ccs_min_build_capacity","ccs_max_age"

    The following file is mandatory, because it sets cost parameters for
    both existing and new project buildouts:

    ccs_build_costs.csv
        CCS_PROJECT, build_year, ccs_overnight_cost, ccs_fixed_om, ccs_variable_om

    """
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'ccs_projects_info.csv'),
        auto_select=True,
        optional_params=["ccs_tech","ccs_load_zone","ccs_capture_rate_mwh_per_ton",
                         "ccs_capacity_limit_mtons","ccs_energy_source", "ccs_fuel_mmbtu_per_ton",
                         "ccs_min_build_capacity","ccs_max_age"],
        index=mod.CCS_PROJECTS,
        param=(mod.ccs_tech, mod.ccs_load_zone, mod.ccs_capture_rate_mwh_per_ton,
               mod.ccs_capacity_limit_mtons, mod.ccs_energy_source,	mod.ccs_fuel_mmbtu_per_ton, 
               mod.ccs_min_build_capacity,	mod.ccs_max_age)
    )

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'dac_credit.csv'),
        autoselect=True,
        index=mod.dac_credit_years,
        param=(mod.dac_credit_value))
    # Construct sets of capacity-limited, ccs-capable and unit-size-specified
    # projects. These sets include projects for which these parameters have
    # a value
    if 'ccs_capacity_limit_mtons' in switch_data.data():
        switch_data.data()['CAPACITY_LIMITED_CCS'] = {
            None: list(switch_data.data(name='ccs_capacity_limit_mtons').keys())}
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'ccs_build_costs.csv'),
        auto_select=True,
        index=mod.CCS_BLD_YRS,
        param=(mod.ccs_overnight_cost, mod.ccs_fixed_om, mod.ccs_variable_om))
    # read FUELS_FOR_MULTIFUEL_GEN from ccs_multiple_fuels.dat if available
    multi_fuels_path = os.path.join(inputs_dir, 'ccs_multiple_fuels.dat')
    if os.path.isfile(multi_fuels_path):
        switch_data.load(filename=multi_fuels_path)