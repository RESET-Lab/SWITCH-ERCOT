from __future__ import division
import os
from pyomo.environ import *
from switch_model.financials import capital_recovery_factor as crf

"""
TO-DO:
    - For SWITCH ERCOT:
        - Convert those tab files to csvs since they have the correct build restrictions
        - I am concerned that the build costs are the same for technologies across all load zones. Is that realistic and wouldnt that just incentivize SWITCH
        to build directly in a load zone and not rely on transmission?
            - Maybe there are other constraints I am not aware of that influence this too
    - add in reserves
    - Add in pre-existing infrastructure (grey hydrogen capacity?)
    - check cost functions (Use data sources from Javier's work to start for cost data, then go from there)
        - Perhaps allow costs to vary by period to account for learning rates and stuff? I think I remember SWITCH being able to do that?
        - The actual cost calculations in SWITCH make complete, but I will need to update and verify the inputs of course
        - Pipeline costs are currently calculated separately, double check this calculate and potentially combine into single cost function
    - PIPELINES:
        - reincorporate pipelines and start adding in pipeline specific information (currently just a clone of the transmission lines)
            - Blending limits for existing pipelines (DONE)
            - Compression
            - Options for retrofitting? AND option to build new pure hydrogen pipelines
        - determine best way to include existing pipeline infrastructure
            - I have a basic version now, take a look at the transmission_lines.csv file
            in Pedro's dataset to see what other information might be relevant
        - Remember to go back and correct the efficiencies
            - How should we be defining pipeline efficiency? It's not like it gets dissipated. Define in terms of pressure drop? Why does that matter though?
    - STORAGE:
        - THE BUG IS STILL HAPPENING
            - Temporary fix: Set the start and end of the period to have zero storage. No there is no storage happening anywhere
        - Possibly change period mass balance to allow period to end with more storage than it started?
            - If I were to do this I would then also have to apply the current TS storage constraint to the TS at the start of the period as well
        - Is 100% efficiency for storage reasonable?
            - Seems right for compressed gas
            - Might be nonzero for liquid storage; I seem to recall it had a nonnegligible boil-off
            - Also injection and withdrawal rates? Look at how SWITCH WECC defines power and energy capacities to get an understanding of the formulation
                - Not SWITCH WECC, but normal switch does not explicitly define the charge/discharge limit as a decision variable. They provide a ratio of the
                energy to power (which would correspond to duration) and used that to define the power limits
        - Need to implement costs for charge/discharge rates, since right now using short-term H2 storage, which I thought was not cost effective?
            - Yeah Bodal has investment costs for both power and energy, but right now we only have energy
            - Define as a decision variable, constraint the withdrawl and storage variables, and add to overall cost
        - Should I track the remaining storage when I go from one timeseries to the next? Would help with long duration storage?
            - YES
        - Ramping rates?
    - ELECTROLYSIS AND FUEL CELLS:
        - Minimum generation rates?
        - Add in SMR (start with no CCS, but add in CCS option eventually)
            - This one is trickier I think, since it will be an issue of fuel consumption as well. Will pull numbers from Bodal for now I think
            - SWITCH fuel consumption calculation: Fuel consumed by generator g in zone z at timepoint t = DispatchGen[g, tp]*gen_full_load_heat_rate[g]
                - So I can setup the build and dispatch decisions of SMR the same way I did the electrolyzers
                - But I need to factor in the costs and emission associated with the natural gas
                - AND I need to determine how to best implement the CCS part of it
            - Inputs needed:
                - Capital cost
                - Fixed O&M
                - Variable O&M
                - Lifetime
                - "Heat rate" (Ratio of incoming thermal energy from fuel (MMBTU) to the amount of hydrogen produced (kg))
                - Define fuel as natural gas specifically (input or internal definition?)
                - Carbon emissions
                - CCS enabled? (T/F)
                    - Or would I make this a decision variable perhaps?
        - Tie electrolysis production to energy produced by renewables
            - Sum of all the electrolysis production must be less than the sum of all the production from variable generation?
                - Done, but still need to figure out how to deal with stored renewables
            - Model was infeasible when I added the constraint to implement green hydrogen
                - Because of the capacity limits on wind and solar in these inputs, feasible when limits are removed
            - Currently debugging
                - Hydrogen demand is being met, but need to verify electricity demand still
                - Done, kind of. There are still some weird storage outcomes, but I think those will be fixed once I can track storage between timeseries
        - Constrain fuel cells so that it cannot consume more hydrogen than what is withdrawn from storage
            -This should already be covered by the mass balance, but will need to double check
                - Difficult to know for sure when the fuel cells are not being used ever in these tests really
            -Be careful to avoid double counting any hydrogen used.
                - In fact you may want to take the hydrogen used for fuel cells out of the mass balance
    - Load profile for hydrogen (transportation?)
        - Would make sense to include industry demand for hydrogen, since that is what most hydrogen is being used for in Texas currently
    - Account for geologic storage

"""


def define_dynamic_lists(m):
    """
    Zone_Power_Injections and Zone_Power_Withdrawals are lists of
    components that contribute to load-zone level power balance equations.
    sum(Zone_Power_Injections[z,t]) == sum(Zone_Power_Withdrawals[z,t])
        for all z,t
    Other modules may append to either list, as long as the components they
    add are indexed by [zone, timepoint] and have units of MW. Other modules
    often include Expressions to summarize decision variables on a zonal basis.
    """

    m.Hydrogen_Injections = []
    m.Hydrogen_Withdrawals = []

def define_dynamic_components(m):
    """
    Adds components to a Pyomo abstract model object to enforce the
    first law of thermodynamics at the level of load zone buses. Unless
    otherwise stated, all terms describing power are in units of MW and
    all terms describing energy are in units of MWh.

    Zone_Energy_Balance[load_zone, timepoint] is a constraint that mandates
    conservation of energy in every load zone and timepoint. This constraint
    sums the model components in the lists Zone_Power_Injections and
    Zone_Power_Withdrawals - each of which is indexed by (z, t) and
    has units of MW - and ensures they are equal. The term tp_duration_hrs
    is factored out of the equation for brevity.
    """

    m.Hnode_Mass_Balance = Constraint(
        m.HNODE_TPS,
        rule=lambda m, z, t: (
            sum(
                getattr(m, component)[z, t]
                for component in m.Hydrogen_Injections
            ) == sum(
                getattr(m, component)[z, t]
                for component in m.Hydrogen_Withdrawals)))

def define_components(m):

    """
    Defining technology details

    Units for hydrogen production and demand are kg/h
    Units for power generation are MW
    Units for storage capacity are MWh
    """

    m.HNODE_TPS = Set(
        dimen=2,
        initialize=lambda m: m.LOAD_ZONES * m.TIMEPOINTS)

    m.hydrogen_demand_kg = Param(
        m.HNODE_TPS)

    m.Hydrogen_Withdrawals.append('hydrogen_demand_kg')

    """
    # electrolyzer details
    """
    m.hydrogen_electrolyzer_capital_cost_per_kgh = Param(m.PERIODS)
    m.hydrogen_electrolyzer_fixed_cost_per_kgh_year = Param(m.PERIODS, default=0.0)
    m.hydrogen_electrolyzer_variable_cost_per_kg = Param(m.PERIODS, default=0.0)  # assumed to include any refurbishment needed
    m.hydrogen_electrolyzer_kg_per_mwh = Param() # assumed to deliver H2 at enough pressure for liquifier and daily buffering
    m.hydrogen_electrolyzer_life_years = Param()
    m.BuildElectrolyzerKgPerHour = Var(m.LOAD_ZONES, m.PERIODS, within=NonNegativeReals)
    m.ElectrolyzerCapacityKgPerHour = Expression(m.LOAD_ZONES, m.PERIODS, rule=lambda m, z, p:
        sum(m.BuildElectrolyzerKgPerHour[z, p_] for p_ in m.CURRENT_AND_PRIOR_PERIODS_FOR_PERIOD[p]))
        #if ((m.period_start[p] - m.period_start[p_]) < m.hydrogen_electrolyzer_life_years)))
    m.DispatchElectrolyzerKgPerHour = Var(m.LOAD_ZONES, m.TIMEPOINTS, within=NonNegativeReals) #by power sector TPS or H node TPS?
    m.Max_Dispatch_Electrolyzer = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        m.DispatchElectrolyzerKgPerHour[z, t] <= m.ElectrolyzerCapacityKgPerHour[z, m.tp_period[t]])

    #Let hydrogen production be a direct decision variable with DispatchElectrolyzer and then calculate the electricity requirements
    m.ConsumeElecMW = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        m.DispatchElectrolyzerKgPerHour[z, t] / m.hydrogen_electrolyzer_kg_per_mwh
    )

    #Define constraint so that the electricity consumed for electrolysis does not exceed the amount of generation from variable sources
    #Two ways we could define this:
        #1. I could define it at the load zone level (electrolysis < renewable gen in that load zone) if we only wanted to allow for local/onsite electrolysis
        #2. Define at a system level (total electrolysis < total renewable gen), which would allow the use of renewable gen from other load zones, but would be hard to track
        #that and would you even know if you were using renewables. I think I will tie it to "on-site" generation for now at the load zone level

    # Load zone level limitation
    #m.Electrolysis_Electricity_Consumption_Limit = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m,z,t:
    #    m.ConsumeElecMW[z,t] <= sum(m.DispatchGen[g, t] for (g, t) in m.GEN_TPS if (g in m.VARIABLE_GENS and g in m.GENS_IN_ZONE[z])))

    #System-wide constraint
    m.Electrolysis_Electricity_Consumption_Limit = Constraint(m.TIMEPOINTS, rule=lambda m, t:
        sum(m.ConsumeElecMW[z,t] for z in m.LOAD_ZONES) <=
        sum(m.DispatchGen[g, t] for (g, t) in m.GEN_TPS if g in m.VARIABLE_GENS))

    m.Hydrogen_Injections.append('DispatchElectrolyzerKgPerHour')
    m.Zone_Power_Withdrawals.append('ConsumeElecMW')

    """
    #SMR with CCS Details
    To-do:
    - Write expression for fuel consumption - Done
    - Decide on how to incorporate CCS - Done
        - Treat blue and grey hydrogen as separate build options for simplicity with different costs and emission intensities --> Done
        - Not sure where the CO2 emissions are used in Switch, but we may want to look into that later --> Calculated in dispatch.py
    - Minimum builds and ramping??
    - Note: Unlike current fuel tracking in Switch, this does not have the option for natural gas to be unavailable in a load zone. Consider adding later
    """

    m.hydrogen_smr_with_ccs_capital_cost_per_kgh = Param(m.PERIODS)
    m.hydrogen_smr_with_ccs_fixed_cost_per_kgh_year = Param(m.PERIODS, default=0.0)
    m.hydrogen_smr_with_ccs_variable_cost_per_kg = Param(m.PERIODS, default=0.0)  # assumed to include any refurbishment needed
    m.hydrogen_smr_with_ccs_kg_per_mmbtu = Param()
    m.hydrogen_smr_with_ccs_life_years = Param()
    m.hydrogen_smr_with_ccs_fuel_source = Param()
    m.hydrogen_smr_with_ccs_co2_per_h2 = Param()
    m.BuildSMRWithCCSKgPerHour = Var(m.LOAD_ZONES, m.PERIODS, within=NonNegativeReals)
    m.SMRWithCCSCapacityKgPerHour = Expression(m.LOAD_ZONES, m.PERIODS, rule=lambda m, z, p:
        sum(m.BuildSMRWithCCSKgPerHour[z, p_] for p_ in m.CURRENT_AND_PRIOR_PERIODS_FOR_PERIOD[p]))
        #if ((p - p_) < m.hydrogen_smr_with_ccs_life_years)))
    m.DispatchSMRWithCCSKgPerHour = Var(m.LOAD_ZONES, m.TIMEPOINTS, within=NonNegativeReals) #by power sector TPS or H node TPS?
    m.Max_Dispatch_SMR_With_CCS = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        m.DispatchSMRWithCCSKgPerHour[z, t] <= m.SMRWithCCSCapacityKgPerHour[z, m.tp_period[t]])

    m.hydrogen_smr_no_ccs_capital_cost_per_kgh = Param(m.PERIODS)
    m.hydrogen_smr_no_ccs_fixed_cost_per_kgh_year = Param(m.PERIODS, default=0.0)
    m.hydrogen_smr_no_ccs_variable_cost_per_kg = Param(m.PERIODS, default=0.0)  # assumed to include any refurbishment needed
    m.hydrogen_smr_no_ccs_kg_per_mmbtu = Param()
    m.hydrogen_smr_no_ccs_life_years = Param()
    m.hydrogen_smr_no_ccs_fuel_source = Param()
    m.hydrogen_smr_no_ccs_co2_per_h2 = Param()
    m.BuildSMRNoCCSKgPerHour = Var(m.LOAD_ZONES, m.PERIODS, within=NonNegativeReals)
    m.SMRNoCCSCapacityKgPerHour = Expression(m.LOAD_ZONES, m.PERIODS, rule=lambda m, z, p:
        sum(m.BuildSMRNoCCSKgPerHour[z, p_] for p_ in m.CURRENT_AND_PRIOR_PERIODS_FOR_PERIOD[p]))
        #if ((p - p_) < m.hydrogen_smr_no_ccs_life_years)))
    m.DispatchSMRNoCCSKgPerHour = Var(m.LOAD_ZONES, m.TIMEPOINTS, within=NonNegativeReals) #by power sector TPS or H node TPS?
    m.Max_Dispatch_SMR_No_CCS = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        m.DispatchSMRNoCCSKgPerHour[z, t] <= m.SMRNoCCSCapacityKgPerHour[z, m.tp_period[t]])

    m.SMRFuelUseRate = Var(m.LOAD_ZONES, m.TIMEPOINTS, within=NonNegativeReals)
    m.SMRFuelUseCalculation = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule = lambda m, z, t:
        m.SMRFuelUseRate[z, t] == m.DispatchSMRWithCCSKgPerHour[z, t]/m.hydrogen_smr_with_ccs_kg_per_mmbtu + m.DispatchSMRNoCCSKgPerHour[z, t]/m.hydrogen_smr_no_ccs_kg_per_mmbtu)

    m.SMRDispatchEmissions = Var(m.LOAD_ZONES, m.TIMEPOINTS, within=NonNegativeReals)
    m.SMREmissionsCalculation = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule = lambda m, z, t:
        m.SMRDispatchEmissions[z, t] == (m.DispatchSMRWithCCSKgPerHour[z, t]*m.hydrogen_smr_with_ccs_co2_per_h2 + m.DispatchSMRNoCCSKgPerHour[z,t]*m.hydrogen_smr_no_ccs_co2_per_h2)/1000) #divide by 1000 to convert to tonnes

    m.Hydrogen_Injections.append('DispatchSMRWithCCSKgPerHour')
    m.Hydrogen_Injections.append('DispatchSMRNoCCSKgPerHour')

    m.H2AnnualEmissions = Expression(m.PERIODS,
        rule=lambda m, period: sum(
            m.SMRDispatchEmissions[z, t] * m.tp_weight_in_year[t]
            for (z, t) in m.ZONE_TIMEPOINTS
            if m.tp_period[t] == period),
        doc="The system's annual emissions from hydrogen production, in metric tonnes of CO2 per year.")

    m.hydrogen_emission_cap = Param(m.PERIODS, within = NonNegativeReals, default = float('inf'))
    m.Enforce_H2_Carbon_Cap = Constraint(m.PERIODS, rule=lambda m, p:
        Constraint.Skip if m.hydrogen_emission_cap[p] == float('inf')
        else (m.H2AnnualEmissions[p]) <= m.hydrogen_emission_cap[p])

    """
    # storage tank details
    # to-do: Add in fixed cost
    """
    m.hydrogen_tank_capital_cost_per_kg = Param(m.PERIODS) #energy investment
    m.hydrogen_tank_capital_cost_per_kg_per_hour = Param(m.PERIODS) #power investment
    m.hydrogen_tank_minimum_size_kg = Param(default=0.0)
    m.hydrogen_tank_life_years = Param()
    m.BuildHydrogenTankKg = Var(m.LOAD_ZONES, m.PERIODS, within=NonNegativeReals) # in kg
    m.BuildHydrogenTankKgPower = Var(m.LOAD_ZONES, m.PERIODS, within=NonNegativeReals) # in kg
    m.HydrogenTankCapacityKg = Expression(m.LOAD_ZONES, m.PERIODS, rule=lambda m, z, p:
        sum(m.BuildHydrogenTankKg[z, p_] for p_ in m.CURRENT_AND_PRIOR_PERIODS_FOR_PERIOD[p]
        if ((p - p_) < m.hydrogen_tank_life_years)))
    m.HydrogenTankStorageRate = Expression(m.LOAD_ZONES, m.PERIODS, rule=lambda m, z, p:
        sum(m.BuildHydrogenTankKgPower[z, p_] for p_ in m.CURRENT_AND_PRIOR_PERIODS_FOR_PERIOD[p]))
        #if ((p - p_) < m.hydrogen_tank_life_years)))
    #m.StoreHydrogenKgPerHour = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, tp:
    #    m.tp_duration_hrs[tp] * m.CompressHydrogenKgPerHour[z, tp])
    m.StoreHydrogenKgPerHour = Var(m.LOAD_ZONES, m.TIMEPOINTS, within=NonNegativeReals)
    m.WithdrawHydrogenKgPerHour = Var(m.LOAD_ZONES, m.TIMEPOINTS, within=NonNegativeReals)
    m.HydrogenStorageStateKg = Var(m.LOAD_ZONES, m.TIMEPOINTS)#, within=NonNegativeReals, initialize=0)

    m.HydrogenNetStorageKg = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m,z,tp:
        (m.StoreHydrogenKgPerHour[z, tp] - m.WithdrawHydrogenKgPerHour[z, tp]) * m.tp_duration_hrs[tp])

    m.HydrogenStoreUpperLimit = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, tp:
        m.StoreHydrogenKgPerHour[z, tp] <= m.HydrogenTankStorageRate[z, m.tp_period[tp]])

    m.HydrogenWithdrawUpperLimit = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, tp:
        m.WithdrawHydrogenKgPerHour[z, tp] <= m.HydrogenTankStorageRate[z, m.tp_period[tp]])

    #m.HydrogenNetStorageLimit = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, tp:
    #    m.HydrogenNetStorageKg[z, tp] <= m.HydrogenTankStorageRate[z, m.tp_period[tp]])

    #Constraint to track storage state, but it breaks for periods with a single timepoint, so I added in something to account for that
    #That likely won't ever come up for real test cases, but it should be included so I can keep testing on the 3_zone_toy
    m.ts_previous = Param(
        m.TIMESERIES,
        within=m.TIMESERIES,
        initialize=lambda m, ts: m.TIMESERIES.prevw(ts))
    def Track_Storage_State_TP_rule(m, z, tp):
        currentTS = m.tp_ts[tp]
        currentPeriod = m.tp_period[tp]

        #for the rare case of a period with a single timepoint
        if m.tp_previous[tp] == tp:
            return m.HydrogenStorageStateKg[z, tp] == (m.StoreHydrogenKgPerHour[z, tp] - m.WithdrawHydrogenKgPerHour[z, tp]) * m.tp_duration_hrs[tp]

        else:
            return m.HydrogenStorageStateKg[z, tp] ==  m.HydrogenStorageStateKg[z, m.tp_previous[tp]] + \
            m.HydrogenNetStorageKg[z, tp]
            #(m.StoreHydrogenKgPerHour[z, tp] - m.WithdrawHydrogenKgPerHour[z, tp]) * m.tp_duration_hrs[tp]

    m.Track_Storage_State_TP = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=Track_Storage_State_TP_rule)

    def starting_tp_check(m, tp):
        currentTS = m.tp_ts[tp]
        currentPeriod = m.tp_period[tp]
        return (m.TPS_IN_TS[currentTS][1]==tp and m.TPS_IN_PERIOD[currentPeriod][1] != tp)

    m.IS_STARTING_TP = Set(initialize=m.TIMEPOINTS, filter=starting_tp_check)
    def Track_Storage_State_TS_rule(m, z, tp):
        currentTS = m.tp_ts[tp]
        previousTS = m.ts_previous[currentTS]
        finalPointOfSeries = m.TPS_IN_TS[previousTS][-1]
        return m.HydrogenStorageStateKg[z, tp] == m.HydrogenStorageStateKg[z, finalPointOfSeries] + \
        m.HydrogenNetStorageKg[z, tp]

    #At first glance this constraint appears to work, but the behavior of the model is MUCH different once this is in place, so I need
    #to look more closely before I can be sure
    m.Track_Storage_State_TS = Constraint(m.LOAD_ZONES, m.IS_STARTING_TP, rule=Track_Storage_State_TS_rule)


    #going to try a slightly different way to track annual storage. If I constrained it such that the first and last timepoint of a period
    #must have a storagestate of 0, I think that might also solve this weird bug?
    #m.Hydrogen_Conservation_of_Mass_Annual = Constraint(m.LOAD_ZONES, m.PERIODS, rule=lambda m, z, p:
    #    m.HydrogenStorageStateKg[z, m.TPS_IN_PERIOD[p][1]] - m.HydrogenStorageStateKg[z, m.TPS_IN_PERIOD[p][-1]] == 0
    #)

    m.Hydrogen_Conservation_of_Mass_Annual_Start = Constraint(m.LOAD_ZONES, m.PERIODS, rule=lambda m, z, p:
        m.HydrogenStorageStateKg[z, m.TPS_IN_PERIOD[p][1]] == 0
    )
    m.Hydrogen_Conservation_of_Mass_Annual_End = Constraint(m.LOAD_ZONES, m.PERIODS, rule=lambda m, z, p:
        m.HydrogenStorageStateKg[z, m.TPS_IN_PERIOD[p][-1]] == 0
    )

    """
    #I think this constraint is unnecessary
    #It is kind of necessary. Serves the same function as restricting StorageState to be nonnegativereals, but with some stuff that does not make sense mathematically
    """
    def Min_Storage_State_rule(m,z,tp):
        if m.tp_previous == tp:
            return m.HydrogenNetStorageKg[z,tp] >= 0
        else:
            return (m.HydrogenStorageStateKg[z, m.tp_previous[tp]] + m.HydrogenNetStorageKg[z, tp]) >=0
    m.Min_Storage_State = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=Min_Storage_State_rule)

    m.Max_Storage_State = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, tp:
        m.HydrogenStorageStateKg[z,tp] <= m.HydrogenTankCapacityKg[z, m.tp_period[tp]])

    #Seems odd to me that the original scaled it to a single year and not to the full period, but I will run with it for now
    #Scaling by year and scaling by period lead to the same outcome in this case
    #m.Hydrogen_Conservation_of_Mass_Annual = Constraint(m.LOAD_ZONES, m.PERIODS, rule=lambda m, z, p:
    #    sum(
    #        (m.StoreHydrogenKgPerHour[z, tp] - m.WithdrawHydrogenKgPerHour[z, tp])
    #            * m.tp_weight[tp]
    #            #* m.tp_weight_in_year[tp]
    #        for tp in m.TPS_IN_PERIOD[p]
    #    ) == 0
    #)

    #Also should definitely have some type of efficiency for these since round trip efficiency is VERY relevant for this
    #Round trip is relevant for fuel cells, not the storage itself really
    m.Hydrogen_Withdrawals.append('StoreHydrogenKgPerHour')
    m.Hydrogen_Injections.append('WithdrawHydrogenKgPerHour')

    """
    # fuel cell details
    """
    m.hydrogen_fuel_cell_capital_cost_per_mw = Param(m.PERIODS)
    m.hydrogen_fuel_cell_fixed_cost_per_mw_year = Param(m.PERIODS, default=0.0)
    m.hydrogen_fuel_cell_variable_cost_per_mwh = Param(m.PERIODS, default=0.0) # assumed to include any refurbishment needed
    m.hydrogen_fuel_cell_mwh_per_kg = Param()
    m.hydrogen_fuel_cell_life_years = Param()
    m.BuildFuelCellMW = Var(m.LOAD_ZONES, m.PERIODS, within=NonNegativeReals)
    m.FuelCellCapacityMW = Expression(m.LOAD_ZONES, m.PERIODS, rule=lambda m, z, p:
        sum(m.BuildFuelCellMW[z, p_] for p_ in m.CURRENT_AND_PRIOR_PERIODS_FOR_PERIOD[p]
        if ((p - p_) < m.hydrogen_fuel_cell_life_years)))
    m.DispatchFuelCellMW = Var(m.LOAD_ZONES, m.TIMEPOINTS, within=NonNegativeReals)
    m.Max_Dispatch_Fuel_Cell = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        m.DispatchFuelCellMW[z, t] <= m.FuelCellCapacityMW[z, m.tp_period[t]])
    m.ConsumeHydrogenKgPerHour = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        m.DispatchFuelCellMW[z, t] / m.hydrogen_fuel_cell_mwh_per_kg
    )
    #m.Hydrogen_For_Fuel_Cells_Limit = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m,z,t:
    #    m.ConsumeHydrogenKgPerHour[z,t] <= m.HydrogenStorageStateKg[z, m.tp_previous[t]])
    m.Hydrogen_Withdrawals.append('ConsumeHydrogenKgPerHour')
    m.Zone_Power_Injections.append('DispatchFuelCellMW')

    """
    # compressor details
    """
    m.hydrogen_compressor_capital_cost_per_kg_per_hour = Param(m.PERIODS)
    m.hydrogen_compressor_fixed_cost_per_kg_hour_year = Param(m.PERIODS, default=0.0)
    m.hydrogen_compressor_variable_cost_per_kg = Param(m.PERIODS, default=0.0)
    m.hydrogen_compressor_mwh_per_kg = Param()
    m.hydrogen_compressor_life_years = Param()

    m.BuildCompressorKgPerHour = Var(m.LOAD_ZONES, m.PERIODS, within=NonNegativeReals)  # capacity to build, measured in kg/hour of throughput
    m.CompressorCapacityKgPerHour = Expression(m.LOAD_ZONES, m.PERIODS, rule=lambda m, z, p:
        sum(m.BuildCompressorKgPerHour[z, p_] for p_ in m.CURRENT_AND_PRIOR_PERIODS_FOR_PERIOD[p]))
        #if ((p - p_) < m.hydrogen_compressor_life_years)))

    #Can expand these two equations to account for the compression used for other technologies like pipelines
    m.Max_Dispatch_Compressor = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        m.StoreHydrogenKgPerHour[z, t] <= m.CompressorCapacityKgPerHour[z, m.tp_period[t]])

    m.CompressHydrogenMW = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        m.StoreHydrogenKgPerHour[z, t] * m.hydrogen_compressor_mwh_per_kg
    )
    m.Zone_Power_Withdrawals.append('CompressHydrogenMW')

    """
    Cost calculations for hydrogen

    """
    m.HydrogenVariableCost = Expression(m.TIMEPOINTS, rule=lambda m, t:
        sum(
            m.DispatchElectrolyzerKgPerHour[z, t] * m.hydrogen_electrolyzer_variable_cost_per_kg[m.tp_period[t]]
            + m.DispatchSMRWithCCSKgPerHour[z, t] * m.hydrogen_smr_with_ccs_variable_cost_per_kg[m.tp_period[t]]
            + m.DispatchSMRNoCCSKgPerHour[z, t] * m.hydrogen_smr_no_ccs_variable_cost_per_kg[m.tp_period[t]]
            + m.SMRFuelUseRate[z, t] * m.fuel_cost[z, m.hydrogen_smr_with_ccs_fuel_source, m.tp_period[t]]
            + m.StoreHydrogenKgPerHour[z, t] * m.hydrogen_compressor_variable_cost_per_kg[m.tp_period[t]]
            + m.DispatchFuelCellMW[z, t] * m.hydrogen_fuel_cell_variable_cost_per_mwh[m.tp_period[t]]
            for z in m.LOAD_ZONES
        )
    )
    m.HydrogenFixedCostAnnual = Expression(m.PERIODS, rule=lambda m, p:
        sum(
            m.ElectrolyzerCapacityKgPerHour[z, p] * (
                m.hydrogen_electrolyzer_capital_cost_per_kgh[p] * crf(m.interest_rate, m.hydrogen_electrolyzer_life_years)
                + m.hydrogen_electrolyzer_fixed_cost_per_kgh_year[p]) * m.hydrogen_electrolyzer_kg_per_mwh
            + m.SMRWithCCSCapacityKgPerHour[z, p] * (
                m.hydrogen_smr_with_ccs_capital_cost_per_kgh[p] * crf(m.interest_rate, m.hydrogen_smr_with_ccs_life_years)
                + m.hydrogen_smr_with_ccs_fixed_cost_per_kgh_year[p])
            + m.SMRNoCCSCapacityKgPerHour[z, p] * (
                m.hydrogen_smr_no_ccs_capital_cost_per_kgh[p] * crf(m.interest_rate, m.hydrogen_smr_no_ccs_life_years)
                + m.hydrogen_smr_no_ccs_fixed_cost_per_kgh_year[p])
            + m.CompressorCapacityKgPerHour[z, p] * (
                m.hydrogen_compressor_capital_cost_per_kg_per_hour[p] * crf(m.interest_rate, m.hydrogen_compressor_life_years)
                + m.hydrogen_compressor_fixed_cost_per_kg_hour_year[p])
            + m.HydrogenTankCapacityKg[z, p] * (
                m.hydrogen_tank_capital_cost_per_kg[p] * crf(m.interest_rate, m.hydrogen_tank_life_years))
            + m.HydrogenTankStorageRate[z, p] * (
                m.hydrogen_tank_capital_cost_per_kg_per_hour[p] * crf(m.interest_rate, m.hydrogen_tank_life_years))
            + m.FuelCellCapacityMW[z, p] * (
                m.hydrogen_fuel_cell_capital_cost_per_mw[p] * crf(m.interest_rate, m.hydrogen_fuel_cell_life_years)
                + m.hydrogen_fuel_cell_fixed_cost_per_mw_year[p])
            for z in m.LOAD_ZONES
        )
    )
    m.Cost_Components_Per_TP.append('HydrogenVariableCost')
    m.Cost_Components_Per_Period.append('HydrogenFixedCostAnnual')


    """
    Code needed for H2 transport via pipelines
    Defining set of connectors, head and tail nodes, length, efficiency, capacity, and builds
    Copied over almost the entire transmission.build file since they are being modeled in the same way

    """

    m.PIPELINES = Set()
    m.pipes_hn1 = Param(m.PIPELINES, within=m.LOAD_ZONES)
    m.pipes_hn2 = Param(m.PIPELINES, within=m.LOAD_ZONES)
    m.pipes_length_km = Param(m.PIPELINES, within=NonNegativeReals)
    m.pipes_efficiency = Param(
        m.PIPELINES,
        within=PercentFraction)
    m.existing_pipes_cap = Param(
        m.PIPELINES,
        within=NonNegativeReals)
    m.blend_limit = Param(
        m.PIPELINES,
        within=PercentFraction
    )
    m.pipes_new_build_allowed = Param(
        m.PIPELINES, within=Boolean, default=True)
    m.PIPELINES_BLD_YRS = Set(
        dimen=2,
        initialize=m.PIPELINES * m.PERIODS,
        filter=lambda m, ptx, p: m.pipes_new_build_allowed[ptx])
    m.BuildPtx = Var(m.PIPELINES_BLD_YRS, within=NonNegativeReals)
    m.pipes_terrain_multiplier = Param(
        m.PIPELINES,
        within=NonNegativeReals,
        default=1)
    m.pipes_capital_cost_per_kg_km = Param(
        within=NonNegativeReals,
        default=1000)
    m.pipes_lifetime_yrs = Param(
        within=NonNegativeReals,
        default=20)
    m.pipes_fixed_om_fraction = Param(
        within=NonNegativeReals,
        default=0.03)
    m.pipes_cost_annual = Param(
        m.PIPELINES,
        within=NonNegativeReals,
        initialize=lambda m, ptx: (
            m.pipes_capital_cost_per_kg_km * m.pipes_terrain_multiplier[ptx] *
            m.pipes_length_km[ptx] * (crf(m.interest_rate, m.pipes_lifetime_yrs) +
            m.pipes_fixed_om_fraction)))
    m.PtxCapacityNameplate = Expression(
        m.PIPELINES, m.PERIODS,
        rule=lambda m, ptx, period: sum(
            m.BuildPtx[ptx, bld_yr]
            for bld_yr in m.PERIODS
            if bld_yr <= period and (ptx, bld_yr) in m.PIPELINES_BLD_YRS
        ) + m.existing_pipes_cap[ptx])
    m.pipes_derating_factor = Param(
        m.TRANSMISSION_LINES,
        within=PercentFraction,
        default=1)
    m.PtxCapacityNameplateAvailableForH2 = Expression(
        m.PIPELINES, m.PERIODS,
        rule=lambda m, ptx, period: (
            m.PtxCapacityNameplate[ptx, period])*m.blend_limit[ptx])# * m.pipes_derating_factor[ptx]))

    def init_DIRECTIONAL_PTX(m):
        ptx_dir = set()
        for ptx in m.PIPELINES:
            ptx_dir.add((m.pipes_hn1[ptx], m.pipes_hn2[ptx]))
            ptx_dir.add((m.pipes_hn2[ptx], m.pipes_hn1[ptx]))
        return ptx_dir
    m.DIRECTIONAL_PTX = Set(
        dimen=2,
        initialize=init_DIRECTIONAL_PTX)
    m.PTX_CONNECTIONS_TO_NODE = Set(
        m.LOAD_ZONES,
        initialize=lambda m, hn: set(
            h for h in m.LOAD_ZONES if (h,hn) in m.DIRECTIONAL_PTX))

    def init_pipes_d_line(m, node_from, node_to):
        for ptx in m.PIPELINES:
            if((m.pipes_hn1[ptx] == node_from and m.pipes_hn2[ptx] == node_to) or
               (m.pipes_hn2[ptx] == node_from and m.pipes_hn1[ptx] == node_to)):
                return ptx
    m.pipes_d_line = Param(
        m.DIRECTIONAL_PTX,
        within=m.PIPELINES,
        initialize=init_pipes_d_line)

    m.PIPES_TIMEPOINTS = Set(
        dimen=3,
        initialize=lambda m: m.DIRECTIONAL_PTX * m.TIMEPOINTS
    )
    m.DispatchPtx = Var(m.PIPES_TIMEPOINTS, within=NonNegativeReals)

    m.Maximum_DispatchPtx = Constraint(
        m.PIPES_TIMEPOINTS,
        rule=lambda m, node_from, node_to, tp: (
            m.DispatchPtx[node_from, node_to, tp] <=
            m.PtxCapacityNameplateAvailableForH2[m.pipes_d_line[node_from, node_to],
                                     m.tp_period[tp]]))

    m.PtxHydrogenSent = Expression(
        m.PIPES_TIMEPOINTS,
        rule=lambda m, node_from, node_to, tp: (
            m.DispatchPtx[node_from, node_to, tp]))
    m.PtxHydrogenReceived = Expression(
        m.PIPES_TIMEPOINTS,
        rule=lambda m, node_from, node_to, tp: (
            m.DispatchPtx[node_from, node_to, tp] *
            m.pipes_efficiency[m.pipes_d_line[node_from, node_to]]))

    def PTXPowerNet_calculation(m, h, tp):
        return (
            sum(m.PtxHydrogenReceived[node_from, h, tp]
                for node_from in m.PTX_CONNECTIONS_TO_NODE[h]) -
            sum(m.PtxHydrogenSent[h, node_to, tp]
                for node_to in m.PTX_CONNECTIONS_TO_NODE[h]))
    m.PTXPowerNet = Expression(
        m.LOAD_ZONES, m.TIMEPOINTS,
        rule=PTXPowerNet_calculation)

    #m.testReceived = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m,h,tp: sum(m.PtxHydrogenReceived[node_from, h, tp]
    #    for node_from in m.PTX_CONNECTIONS_TO_NODE[h]))
    #m.testSent = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m,h,tp: sum(m.PtxHydrogenSent[h, node_to, tp]
    #    for node_to in m.PTX_CONNECTIONS_TO_NODE[h]))

    # Register net transmission as contributing to zonal energy balance
    m.Hydrogen_Injections.append('PTXPowerNet')

    #m.Hydrogen_Injections.append('testReceived')
    #m.Hydrogen_Withdrawals.append('testSent')

    # An expression to summarize annual costs for the objective
    # function. Units should be total annual future costs in $base_year
    # real dollars. The objective function will convert these to
    # base_year Net Present Value in $base_year real dollars.
    m.PtxFixedCosts = Expression(
        m.PERIODS,
        rule=lambda m, p: sum(
            m.PtxCapacityNameplate[ptx, p] * m.pipes_cost_annual[ptx]
            for ptx in m.PIPELINES
        )
    )
    m.Cost_Components_Per_Period.append('PtxFixedCosts')


def load_inputs(m, switch_data, inputs_dir):

    """
    Need input files for: Hydrogen load zones, hydrogen time points, hydrogen projects (gen and storage), pipeline transport

    LOAD_ZONES.csv: LOAD_ZONES hydrogen_dbid existing_local_td coupled_load_zones, transmission stuff,
    need to add a column to load zone indicating what hydrogen node they are coupled to i think
    hydrogen_load.csv

    have to add a column to load_zones.csv to indicate the hydrogen node they are coupled to

    Do we need to define individual projects like generation_projects_info, or can we just provide general information and let the model decide how much to build?
        I think the latter but could be worth a broader discussion with Sergio
        Based on what Jeff Reeds was saying, electrolyzers (and presumably fuel cells?) are very compact and modular, so maybe it makes more sense to get the model flexibility
        in deciding how much is built instead of working off of projects.

    Do we want a local_td equivalent for hydrogen pipelines? Might be nice to include some representation of the distribution level of pipelines, even if it is simple

    Example syntax for loading input data below
    """
    #Can use this same file from original hydrogen module for technology information
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'hydrogen_generation_costs.csv'),
        optional=False,
        auto_select=True,
        index = m.PERIODS,
        param=(
            m.hydrogen_electrolyzer_capital_cost_per_kgh,
            m.hydrogen_electrolyzer_fixed_cost_per_kgh_year,
            m.hydrogen_electrolyzer_variable_cost_per_kg,
            m.hydrogen_fuel_cell_capital_cost_per_mw,
            m.hydrogen_fuel_cell_fixed_cost_per_mw_year,
            m.hydrogen_fuel_cell_variable_cost_per_mwh,
            m.hydrogen_compressor_capital_cost_per_kg_per_hour,
            m.hydrogen_compressor_fixed_cost_per_kg_hour_year,
            m.hydrogen_compressor_variable_cost_per_kg,
            m.hydrogen_tank_capital_cost_per_kg,
            m.hydrogen_tank_capital_cost_per_kg_per_hour,
            m.hydrogen_smr_with_ccs_capital_cost_per_kgh,
            m.hydrogen_smr_with_ccs_fixed_cost_per_kgh_year,
            m.hydrogen_smr_with_ccs_variable_cost_per_kg,
            m.hydrogen_smr_no_ccs_capital_cost_per_kgh,
            m.hydrogen_smr_no_ccs_fixed_cost_per_kgh_year,
            m.hydrogen_smr_no_ccs_variable_cost_per_kg,
        )
    )
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'hydrogen_emissions.csv'),
        optional=False,
        #optional_params=(m.hydrogen_emission_cap),
        auto_select=True,
        index=m.PERIODS,
        param=(m.hydrogen_emission_cap))

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'hydrogen_generation_params.csv'),
        optional=False,
        auto_select=True,
        param=(
            m.hydrogen_electrolyzer_kg_per_mwh,
            m.hydrogen_electrolyzer_life_years,
            m.hydrogen_fuel_cell_life_years,
            m.hydrogen_fuel_cell_mwh_per_kg,
            m.hydrogen_compressor_life_years,
            m.hydrogen_compressor_mwh_per_kg,
            m.hydrogen_tank_life_years,
            m.hydrogen_tank_minimum_size_kg,
            m.hydrogen_smr_with_ccs_kg_per_mmbtu,
            m.hydrogen_smr_with_ccs_life_years,
            m.hydrogen_smr_with_ccs_fuel_source,
            m.hydrogen_smr_with_ccs_co2_per_h2,
            m.hydrogen_smr_no_ccs_kg_per_mmbtu,
            m.hydrogen_smr_no_ccs_life_years,
            m.hydrogen_smr_no_ccs_fuel_source,
            m.hydrogen_smr_no_ccs_co2_per_h2
        )
    )
    #switch_data.load_aug(
    #    filename=os.path.join(inputs_dir, 'LOAD_ZONES.csv'),
    #    optional=False, auto_select=True,
    #    index=m.LOAD_ZONES,
    #    param=(
    #        m.coupled_load_zones,
    #        m.hydrogen_dbid
    #    )
    #)

    #switch_data.load_aug(
    #    filename=os.path.join(inputs_dir, 'load_zones.csv'),
    #    optional=False, auto_select=True,
    #    index=m.LOAD_ZONES,
    #    param=(
    #        m.coupled_LOAD_ZONES
    #    )
    #)

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'hydrogen_loads.csv'),
        optional=False, auto_select=True,
        param=(
            m.hydrogen_demand_kg,
        )
    )

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'pipelines.csv'),
        auto_select=True,
        index=m.PIPELINES,
        optional_params=(
            'pipes_dbid', 'pipes_derating_factor',
            'pipes_terrain_multiplier', 'pipes_new_build_allowed'
        ),
        param=(
            m.pipes_hn1, m.pipes_hn2,
            m.pipes_length_km, m.pipes_efficiency, m.existing_pipes_cap, m.blend_limit,
        )
    )
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'pipes_params.csv'),
        optional=True, auto_select=True,
        param=(
            m.pipes_capital_cost_per_kg_km, m.pipes_lifetime_yrs,
            m.pipes_fixed_om_fraction
        )
    )


def post_solve(instance, outdir):
    """
    Export storage build information to storage_builds.csv, and storage
    dispatch info to storage_dispatch.csv

    """

    import switch_model.reporting as reporting

    #Misc output table I can use for debugging
    #reporting.write_table(
    #    instance, instance.IS_STARTING_TP,
    #    output_file=os.path.join(outdir, "timetest.csv"),
    #    headings=("TP", 'TS'),
    #    values=lambda m, tp: (
    #        tp, m.tp_ts[tp]
    #    ))

    reporting.write_table(
        instance, instance.LOAD_ZONES*instance.PERIODS,
        output_file=os.path.join(outdir, "hydrogen_build.csv"),
        headings=("Load Zone", "period", "BuildElectrolyzerKgPerHour", "BuildSMRWithCCSKgPerHour",
         "BuildSMRNoCCSKgPerHour", "BuildFuelCellMW", "BuildHydrogenTankKg", "BuildHydrogenTankKgPerHour",
         "SMRWithCCSCapacityKgPerHour", "SMRNoCCSCapacityKgPerHour"),
        values=lambda m, z, p: (
            z, p, m.BuildElectrolyzerKgPerHour[z, p], m.BuildSMRWithCCSKgPerHour[z,p],
             m.BuildSMRNoCCSKgPerHour[z,p], m.BuildFuelCellMW[z, p], m.BuildHydrogenTankKg[z, p], m.HydrogenTankStorageRate[z, p],
             m.SMRWithCCSCapacityKgPerHour[z, p], m.SMRNoCCSCapacityKgPerHour[z, p]
            ))

    reporting.write_table(
        instance, instance.LOAD_ZONES*instance.PERIODS,
        output_file=os.path.join(outdir, "hydrogen_capacities.csv"),
        headings=("Load Zone", "period", "SMRWithCCSCapacityKgPerHour", "SMRNoCCSCapacityKgPerHour",
        "ElectrolyzerCapacityKgPerHour", "HydrogenTankCapacityKg", "HydrogenTankStorageRate",
        "FuelCellCapacityMW"),
        values=lambda m, z, p: (
            z, p, m.SMRWithCCSCapacityKgPerHour[z, p], m.SMRNoCCSCapacityKgPerHour[z, p],
            m.ElectrolyzerCapacityKgPerHour[z, p], m.HydrogenTankCapacityKg[z, p], m.HydrogenTankStorageRate[z, p],
            m.FuelCellCapacityMW[z, p]
            ))

    reporting.write_table(
        instance, instance.PERIODS,
        output_file=os.path.join(outdir, "h2_annual_emissions.csv"),
        headings=("period", "H2AnnualEmissions"),
        values=lambda m, p: (
            p, m.H2AnnualEmissions[p]
            ))

    reporting.write_table(
        instance, instance.HNODE_TPS,
        output_file=os.path.join(outdir, "hydrogen_dispatch.csv"),
        headings=("Load Zone", "timepoint", "StoreHydrogenKgPerHour", "WithdrawHydrogenKgPerHour",
        "DispatchElectrolyzerKgPerHour", "DispatchSMRWithCCSKgPerHour", "DispatchSMRNoCCSKgPerHour", "ConsumeElecMW",
        "ConsumeHydrogenKgPerHour","DispatchFuelCellMW", "HydrogenNetStorageKg", "HydrogenStorageStateKg", "NetTransmission"),
        values=lambda m, z, t: (
            z, m.tp_timestamp[t], m.StoreHydrogenKgPerHour[z, t], m.WithdrawHydrogenKgPerHour[z, t],
             m.DispatchElectrolyzerKgPerHour[z, t], m.DispatchSMRWithCCSKgPerHour[z,t], m.DispatchSMRNoCCSKgPerHour[z,t],
             m.ConsumeElecMW[z,t], m.ConsumeHydrogenKgPerHour[z,t], m.DispatchFuelCellMW[z, t],
             m.HydrogenNetStorageKg[z,t], m.HydrogenStorageStateKg[z,t], m.PTXPowerNet[z,t]
            ))
