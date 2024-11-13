"""
This module defines flexible dispatch for geothermal. It builds on top of generic
generators, adding components for deciding how much energy to build into
storage, when to store energy, energy accounting, etc.
"""

from pyomo.environ import *
import os, collections
from switch_model.financials import capital_recovery_factor as crf

dependencies = 'switch_model.timescales', 'switch_model.balancing.load_zones',\
    'switch_model.financials', 'switch_model.energy_sources.properties', \
    'switch_model.generators.core.build', 'switch_model.generators.core.dispatch'

def define_components(mod):
    """
    dghfgfgh
    """

    mod.gen_is_flexible = Param(mod.GENERATION_PROJECTS, within=Boolean, default=False)
    mod.FLEXIBLE_GENS = Set(
        initialize=mod.GENERATION_PROJECTS,
        filter=lambda m, g: m.gen_is_flexible[g])
    mod.FLEXIBLE_GEN_PERIODS = Set(
        within=mod.GEN_PERIODS,
        initialize=lambda m: [(g, p) for g in m.FLEXIBLE_GENS for p in m.PERIODS_FOR_GEN[g]]
    )
    # I don't think I need this, just use the heat rate? NO efficiency is more appropriate conceptually
    mod.flexible_gen_efficiency = Param(
        mod.FLEXIBLE_GENS,
        within=PercentFraction,
        default=1.0)
    # TODO: rename to gen_charge_to_discharge_ratio?
    mod.flex_store_to_release_ratio = Param(
        mod.FLEXIBLE_GENS,
        within=NonNegativeReals,
        default=1.0)
    # mod.flex_gen_max_cycles_per_year = Param(
    #     mod.FLEXIBLE_GENS,
    #     within=NonNegativeReals,
    #     default=float('inf'))

    # TODO: build this set up instead of filtering down, to improve performance
    mod.FLEXIBLE_GEN_BLD_YRS = Set(
        dimen=2,
        initialize=mod.GEN_BLD_YRS,
        filter=lambda m, g, bld_yr: g in m.FLEXIBLE_GENS)
    mod.gen_flex_energy_overnight_cost = Param(
        mod.FLEXIBLE_GEN_BLD_YRS,
        within=NonNegativeReals)

    mod.BuildFlexibleEnergy = Var(
        mod.FLEXIBLE_GEN_BLD_YRS,
        within=NonNegativeReals)
    # Summarize capital costs of energy storage for the objective function
    # Note: A bug in to 2.0.0b3 - 2.0.5, assigned costs that were several times
    # too high
    # Do I need this for geothermal? Isn't energy capacity just the ground lol
    mod.FlexibleEnergyCapitalCost = Expression(
        mod.FLEXIBLE_GENS, mod.PERIODS,
        rule=lambda m, g, p: sum(
            m.BuildFlexibleEnergy[g, bld_yr]
            * m.gen_flex_energy_overnight_cost[g, bld_yr]
            * crf(m.interest_rate, m.gen_max_age[g])
            for bld_yr in m.BLD_YRS_FOR_GEN_PERIOD[g, p]
        )
    )
    mod.FlexibleEnergyFixedCost = Expression(
        mod.PERIODS,
        rule=lambda m, p: sum(
            m.FlexibleEnergyCapitalCost[g, p] for g in m.FLEXIBLE_GENS
        )
    )
    mod.Cost_Components_Per_Period.append('FlexibleEnergyFixedCost')

    # 2.0.0b3 code:
    # mod.FlexibleEnergyInstallCosts = Expression(
    # mod.PERIODS,
    # rule=lambda m, p: sum(m.BuildFlexibleEnergy[g, bld_yr] *
    #            m.gen_Flexible_energy_overnight_cost[g, bld_yr] *
    #            crf(m.interest_rate, m.gen_max_age[g])
    #            for (g, bld_yr) in m.Flexible_GEN_BLD_YRS))

    mod.FlexibleEnergyCapacity = Expression(
        mod.FLEXIBLE_GENS, mod.PERIODS,
        rule=lambda m, g, period: sum(
            m.BuildFlexibleEnergy[g, bld_yr]
            for bld_yr in m.BLD_YRS_FOR_GEN_PERIOD[g, period]
        )
    )

    mod.FLEXIBLE_GEN_TPS = Set(
        dimen=2,
        initialize=lambda m: (
            (g, tp)
                for g in m.FLEXIBLE_GENS
                    for tp in m.TPS_FOR_GEN[g]))

    mod.ChargeFlexible = Var(
        mod.FLEXIBLE_GEN_TPS,
        within=NonNegativeReals)

    # use fixed energy/power ratio (# hours of capacity) when specified
    # mod.Enforce_Fixed_Energy_Flexible_Ratio = Constraint(
    #     mod.FLEXIBLE_GEN_BLD_YRS,
    #     rule=lambda m, g, y:
    #         Constraint.Skip if m.gen_storage_energy_to_power_ratio[g] == float("inf") # no value specified
    #         else
    #         (m.BuildFlexibleEnergy[g, y] == m.gen_storage_energy_to_power_ratio[g] * m.BuildGen[g, y])
    # )

    # wondering if we need different inputs for "charging" or "discharging" in this instance

    def Charge_Flexible_Upper_Limit_rule(m, g, t):
        return m.ChargeFlexible[g,t] <= \
            m.DispatchUpperLimit[g, t] * m.flex_store_to_release_ratio[g]
    mod.Charge_Flexible_Upper_Limit = Constraint(
        mod.FLEXIBLE_GEN_TPS,
        rule=Charge_Flexible_Upper_Limit_rule)

    mod.FlexStateOfCharge = Var(
        mod.FLEXIBLE_GEN_TPS,
        within=NonNegativeReals)

    def Track_Flex_State_Of_Charge_rule(m, g, t):
        return m.FlexStateOfCharge[g, t] == \
            m.FlexStateOfCharge[g, m.tp_previous[t]] + \
            (m.ChargeFlexible[g, t] * m.flexible_gen_efficiency[g] -
             m.DispatchGen[g, t]) * m.tp_duration_hrs[t]
    mod.Flex_Track_State_Of_Charge = Constraint(
        mod.FLEXIBLE_GEN_TPS,
        rule=Track_Flex_State_Of_Charge_rule)

    def Flex_State_Of_Charge_Upper_Limit_rule(m, g, t):
        return m.FlexStateOfCharge[g, t] <= \
            m.FlexibleEnergyCapacity[g, m.tp_period[t]]
    mod.State_Of_Charge_Upper_Limit = Constraint(
        mod.FLEXIBLE_GEN_TPS,
        rule=Flex_State_Of_Charge_Upper_Limit_rule)



def load_inputs(mod, switch_data, inputs_dir):
    """

    Import storage parameters. Optional columns are noted with a *.

    generation_projects_info.csv
        GENERATION_PROJECT, ...
        gen_storage_efficiency, gen_store_to_release_ratio*,
        gen_storage_energy_to_power_ratio*, gen_storage_max_cycles_per_year*

    gen_build_costs.csv
        GENERATION_PROJECT, build_year, ...
        gen_storage_energy_overnight_cost

    """

    # TODO: maybe move these columns to a storage_gen_info file to avoid the weird index
    # reading and avoid having to create these extra columns for all projects;
    # Alternatively, say that these values are specified for _all_ projects (maybe with None
    # as default) and then define FLEXIBLE_GENS as the subset of projects for which
    # gen_storage_efficiency has been specified, then require valid settings for all
    # FLEXIBLE_GENS.
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'generation_projects_info.csv'),
        auto_select=True,
        optional_params=['gen_is_flexible'],
        param=(mod.gen_is_flexible))
    # Base the set of storage projects on storage efficiency being specified.
    # TODO: define this in a more normal way

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'gen_build_costs.csv'),
        auto_select=True,
        param=(mod.gen_flex_energy_overnight_cost))
        #

def post_solve(instance, outdir):
    """
    Export storage build information to storage_builds.csv, and storage
    dispatch info to storage_dispatch.csv
    """
    import switch_model.reporting as reporting

    """
    New version of the write_table functions to split storage build and capacity (copied from Switch WECC)
    """
    # # Write how much is built each build year for each project to storage_builds.csv
    # reporting.write_table(
    #     instance, instance.FLEXIBLE_GEN_BLD_YRS,
    #     output_file=os.path.join(outdir, "storage_builds.csv"),
    #     headings=("generation_project", "build_year", "load_zone",
    #               "IncrementalPowerCapacityMW", "IncrementalEnergyCapacityMWh"),
    #     values=lambda m, g, bld_yr: (
    #         g, bld_yr, m.gen_load_zone[g],
    #         m.BuildGen[g, bld_yr], m.BuildFlexibleEnergy[g, bld_yr],
    #         ))
    # # Write the total capacity for each project at each period to storage_capacity.csv
    # reporting.write_table(
    #     instance, instance.FLEXIBLE_GEN_PERIODS,
    #     output_file=os.path.join(outdir, "storage_capacity.csv"),
    #     headings=("generation_project", "period", "load_zone",
    #               "OnlinePowerCapacityMW", "OnlineEnergyCapacityMWh"),
    #     values=lambda m, g, p: (
    #         g, p, m.gen_load_zone[g],
    #         m.GenCapacity[g, p], m.FlexibleEnergyCapacity[g, p])
    # )
    """
    Original write_table function that was breaking
    """
    """
    reporting.write_table(
        instance, instance.FLEXIBLE_GEN_BLD_YRS,
        output_file=os.path.join(outdir, "storage_builds.csv"),
        headings=("generation_project", "period", "load_zone",
                  "IncrementalPowerCapacityMW", "IncrementalEnergyCapacityMWh",
                  "OnlinePowerCapacityMW", "OnlineEnergyCapacityMWh" ),
        values=lambda m, g, bld_yr: (
            g, bld_yr, m.gen_load_zone[g],
            m.BuildGen[g, bld_yr], m.BuildFlexibleEnergy[g, bld_yr],
            m.GenCapacity[g, bld_yr], m.FlexibleEnergyCapacity[g, bld_yr]
            ))
    """

    reporting.write_table(
        instance, instance.FLEXIBLE_GEN_TPS,
        output_file=os.path.join(outdir, "storage_dispatch.csv"),
        headings=("generation_project", "timepoint", "load_zone",
                  "FlexChargeMW", "DischargeMW", "flexStateOfCharge"),
        values=lambda m, g, t: (
            g, m.tp_timestamp[t], m.gen_load_zone[g],
            m.ChargeFlexible[g, t], m.DispatchGen[g, t],
            m.FlexStateOfCharge[g, t]
            ))
