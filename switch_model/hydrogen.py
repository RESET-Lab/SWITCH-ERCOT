from __future__ import division
import os
from pyomo.environ import *
from switch_model.financials import capital_recovery_factor as crf

"""
TO-DO:
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
            - Temporary fix: Set the start and end of the period to have zero storage. No there is no storage happening anywhere -> 2023 update: Go see if this is fixed
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
                - In fact you may want to take the hydrogen used for fuel cells out of the mass balance - why would I do that?
    - Account for geologic storage
    - Should I split these things into different modules? I think eventually that would be wise but for now we can consolidate

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

    m.HYDROGEN_GEN = Set()
    m.H2_GEN_PERIODS = Set(dimen=2)
    m.H2_GEN_BLD_YRS = Set(dimen=3, initialize=m.HYDROGEN_GEN*m.LOAD_ZONES*m.PERIODS)

    m.kg_per_unit = Param(m.HYDROGEN_GEN)
    m.life_years = Param(m.HYDROGEN_GEN)
    m.co2_per_kg = Param(m.HYDROGEN_GEN)
    m.h2_gen_fuel = Param(m.HYDROGEN_GEN)
    m.IS_ELECTROLYSIS = Param(m.HYDROGEN_GEN, within=Boolean)

    m.capital_cost_per_kgh = Param(m.H2_GEN_PERIODS)
    m.fixed_cost_per_kgh_year = Param(m.H2_GEN_PERIODS)
    m.variable_cost_per_kg = Param(m.H2_GEN_PERIODS)
    
    m.ELECTROLYSIS_GEN = Set(initialize=m.HYDROGEN_GEN, filter=lambda m, g:m.IS_ELECTROLYSIS[g])
    m.NG_GEN = Set(initialize=m.HYDROGEN_GEN, filter= lambda m, g: m.IS_ELECTROLYSIS[g]==False)

    m.BuildH2Gen = Var(m.H2_GEN_BLD_YRS, within=NonNegativeReals)

    m.H2GenCapacity = Expression(m.H2_GEN_BLD_YRS, rule=lambda m, g, z, p:
        sum(m.BuildH2Gen[g, z, p_] for p_ in m.CURRENT_AND_PRIOR_PERIODS_FOR_PERIOD[p]
        if ((m.period_start[p] - m.period_start[p_]) < m.life_years[g])))
    
    m.DispatchH2Gen = Var(m.HYDROGEN_GEN, m.LOAD_ZONES, m.TIMEPOINTS, within=NonNegativeReals)

    m.Max_H2_Dispatch = Constraint(m.HYDROGEN_GEN, m.LOAD_ZONES, m.TIMEPOINTS, rule=
        lambda m, g, z, t: m.DispatchH2Gen[g, z, t] <= m.H2GenCapacity[g, z, m.tp_period[t]])

    m.ConsumeElecMW = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        sum(m.DispatchH2Gen[g, z, t] / m.kg_per_unit[g] for g in m.ELECTROLYSIS_GEN)
    )

    #System-wide constraint
    
    #Hourly matching, with geographical matching and additionality
    m.Electrolysis_Electricity_Consumption_Limit = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        m.ConsumeElecMW[z,t] <= sum(m.DispatchGen[g, t] for g in m.VARIABLE_GENS_IN_ZONE[z] if (g,m.tp_period[t]) in m.NEW_GEN_BLD_YRS))
    

    #Annual matching
    #m.Electrolysis_Electricity_Consumption_Limit = Constraint(m.PERIODS, rule=lambda m, p:
    #    sum(m.ConsumeElecMW[z,t] for (z in m.LOAD_ZONES) and (t in m.TPS_IN_PERIOD[p])) <=
    #    sum(m.DispatchGen[g, t] for (g, t) in m.GEN_TPS if g in m.VARIABLE_GENS and t is in m.TPS_IN_PERIOD[p]))

    m.H2ZonalDispatch = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        sum(m.DispatchH2Gen[g, z , t] for g in m.HYDROGEN_GEN))

    m.Hydrogen_Injections.append('H2ZonalDispatch')
    m.Zone_Power_Withdrawals.append('ConsumeElecMW')

    m.H2FuelUseRate = Var(m.NG_GEN, m.LOAD_ZONES, m.TIMEPOINTS, within=NonNegativeReals)
    m.H2FuelUseCalculation = Constraint(m.NG_GEN, m.LOAD_ZONES, m.TIMEPOINTS, rule = lambda m, g, z, t:
        m.H2FuelUseRate[g, z, t] == m.DispatchH2Gen[g, z, t]/m.kg_per_unit[g])

    m.H2DispatchEmissions = Var(m.HYDROGEN_GEN, m.LOAD_ZONES, m.TIMEPOINTS, within=NonNegativeReals)
    m.H2EmissionsCalculation = Constraint(m.HYDROGEN_GEN, m.LOAD_ZONES, m.TIMEPOINTS, rule = lambda m, g, z, t:
        m.H2DispatchEmissions[g, z, t] == m.DispatchH2Gen[g, z, t]*m.co2_per_kg[g]) 
    
    m.H2AnnualEmissions = Expression(m.PERIODS,
        rule=lambda m, period: sum(
            sum(m.H2DispatchEmissions[g, z, t] for g in m.HYDROGEN_GEN) * m.tp_weight_in_year[t] / 1000 #Convert from kg to tonnes
            for (z, t) in m.ZONE_TIMEPOINTS
            if m.tp_period[t] == period),
        doc="The system's annual emissions from hydrogen production, in metric tonnes of CO2 per year.")

    m.System_Emissions.append('H2AnnualEmissions')

    m.hydrogen_emission_cap = Param(m.PERIODS, within = NonNegativeReals, default = float('inf'))
    m.Enforce_H2_Carbon_Cap = Constraint(m.PERIODS, rule=lambda m, p:
        Constraint.Skip if m.hydrogen_emission_cap[p] == float('inf')
        else (m.H2AnnualEmissions[p]) <= m.hydrogen_emission_cap[p])
    m.H2EmissionsCosts = Expression(m.PERIODS,
    rule=lambda m, p: 
        (m.H2AnnualEmissions[p])* m.carbon_cost_dollar_per_tco2[p],
    doc=("Enforces the carbon cap for hydrogen-related emissions."))
    m.Cost_Components_Per_Period.append('H2EmissionsCosts')

    """
    # fuel cell details
    """
    m.H2_CONVERTERS = Set()
    m.H2_CONV_PERIODS = Set(dimen=2)
    m.H2_CONV_BLD_YRS = Set(dimen=3, initialize=m.H2_CONVERTERS*m.LOAD_ZONES*m.PERIODS)
    m.hydrogen_conv_capital_cost_per_mw = Param(m.H2_CONVERTERS, m.PERIODS)
    m.hydrogen_conv_fixed_cost_per_mw_year = Param(m.H2_CONVERTERS, m.PERIODS, default=0.0)
    m.hydrogen_conv_variable_cost_per_mwh = Param(m.H2_CONVERTERS, m.PERIODS, default=0.0) # assumed to include any refurbishment needed
    m.hydrogen_conv_mwh_per_kg = Param(m.H2_CONVERTERS)
    m.hydrogen_conv_life_years = Param(m.H2_CONVERTERS)
    m.BuildConverterMW = Var(m.H2_CONV_BLD_YRS, within=NonNegativeReals)
    m.ConverterCapacityMW = Expression(m.H2_CONV_BLD_YRS, rule=lambda m, g, z, p:
        sum(m.BuildConverterMW[g, z, p_] for p_ in m.CURRENT_AND_PRIOR_PERIODS_FOR_PERIOD[p]
        if ((p - p_) < m.hydrogen_conv_life_years)))
    m.DispatchConverterMW = Var(m.H2_CONVERTERS, m.LOAD_ZONES, m.TIMEPOINTS, within=NonNegativeReals)
    m.Max_Dispatch_Converter = Constraint(m.H2_CONVERTERS, m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, g, z, t:
        m.DispatchConverterMW[g, z, t] <= m.ConverterCapacityMW[g, z, m.tp_period[t]])
    m.ConsumeHydrogenKgPerHour = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        sum(m.DispatchConverterMW[g, z, t] / m.hydrogen_conv_mwh_per_kg[g] for g in m.H2_CONVERTERS)
    )
    #m.Hydrogen_For_Fuel_Cells_Limit = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m,z,t:
    #    m.ConsumeHydrogenKgPerHour[z,t] <= m.HydrogenStorageStateKg[z, m.tp_previous[t]])
    m.Hydrogen_Withdrawals.append('ConsumeHydrogenKgPerHour')
    m.Zone_Power_Injections.append('DispatchConverterMW')


    """
    # storage tank details
    """
    m.H2_STORAGE_PROJECTS = Set()
    m.H2_STORAGE_BUILD_YRS = Set(
        dimen=2,
        initialize = m.H2_STORAGE_PROJECTS * m.PERIODS
        )
    m.H2_STORAGE_BUILD_YRS_ZONES = Set(
        dimen=3,
        initialize = m.H2_STORAGE_PROJECTS * m.LOAD_ZONES * m.PERIODS
        )
    m.H2_STORAGE_TPS = Set(
        dimen=3,
        initialize = m.H2_STORAGE_PROJECTS * m.LOAD_ZONES * m.TIMEPOINTS
        )

    m.h2_storage_capital_cost_per_kg = Param(m.H2_STORAGE_BUILD_YRS) #energy investment
    m.h2_storage_capital_cost_per_kg_per_hour = Param(m.H2_STORAGE_BUILD_YRS) #power investment
    m.h2_storage_fixed_om_per_kg = Param(m.H2_STORAGE_BUILD_YRS) #energy investment
    m.h2_storage_fixed_om_per_kg_per_hour = Param(m.H2_STORAGE_BUILD_YRS) #power investment

    m.h2_storage_minimum_size_kg = Param(m.H2_STORAGE_PROJECTS, within=NonNegativeReals, default=0.0)
    m.h2_storage_life_years = Param(m.H2_STORAGE_PROJECTS, default=20)
    m.BuildH2StorageKg = Var(m.H2_STORAGE_BUILD_YRS_ZONES, within=NonNegativeReals) # in kg
    m.BuildH2StorageKgPower = Var(m.H2_STORAGE_BUILD_YRS_ZONES, within=NonNegativeReals) # in kg
    m.H2StorageCapacityKg = Expression(m.H2_STORAGE_BUILD_YRS_ZONES, rule=lambda m, s, z, p:
        sum(m.BuildH2StorageKg[s, z, p_] for p_ in m.CURRENT_AND_PRIOR_PERIODS_FOR_PERIOD[p]
        if ((p - p_) < m.h2_storage_life_years[s])))
    m.H2StorageCapacityKgPerHour = Expression(m.H2_STORAGE_BUILD_YRS_ZONES, rule=lambda m, s, z, p:
        sum(m.BuildH2StorageKgPower[s, z, p_] for p_ in m.CURRENT_AND_PRIOR_PERIODS_FOR_PERIOD[p]
        if ((p - p_) < m.h2_storage_life_years[s])))
    #m.StoreHydrogenKgPerHour = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, tp:
    #    m.tp_duration_hrs[tp] * m.CompressHydrogenKgPerHour[z, tp])

    # Defining minimum build capacities
    m.H2BuildMinStorageCap = Var(
        m.H2_STORAGE_BUILD_YRS_ZONES,
        within=Binary)
    m.Enforce_Min_Build_Lower_H2 = Constraint(
        m.H2_STORAGE_BUILD_YRS_ZONES,
        rule=lambda m, g, z, p: (
            m.H2BuildMinStorageCap[g, z, p] * m.h2_storage_minimum_size_kg[g]
            <= m.BuildH2StorageKg[g, z, p]))


    m.StoreHydrogenKgPerHour = Var(m.H2_STORAGE_TPS, within=NonNegativeReals)
    m.WithdrawHydrogenKgPerHour = Var(m.H2_STORAGE_TPS, within=NonNegativeReals)
    m.HydrogenStorageStateKg = Var(m.H2_STORAGE_TPS)#, within=NonNegativeReals, initialize=0)

    m.HydrogenNetStorageKg = Expression(m.H2_STORAGE_TPS, rule=lambda m, s, z, tp:
        (m.StoreHydrogenKgPerHour[s, z, tp] - m.WithdrawHydrogenKgPerHour[s, z, tp]) * m.tp_duration_hrs[tp])

    m.HydrogenStoreUpperLimit = Constraint(m.H2_STORAGE_TPS, rule=lambda m, s, z, tp:
        m.StoreHydrogenKgPerHour[s, z, tp] <= m.H2StorageCapacityKgPerHour[s, z, m.tp_period[tp]])

    m.HydrogenWithdrawUpperLimit = Constraint(m.H2_STORAGE_TPS, rule=lambda m, s, z, tp:
        m.WithdrawHydrogenKgPerHour[s, z, tp] <= m.H2StorageCapacityKgPerHour[s, z, m.tp_period[tp]])

    #m.HydrogenNetStorageLimit = Constraint(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, tp:
    #    m.HydrogenNetStorageKg[z, tp] <= m.H2StorageCapacityKgPerHour[z, m.tp_period[tp]])

    
    #Constraint to track storage state, but it breaks for periods with a single timepoint, so I added in something to account for that
    #That likely won't ever come up for real test cases, but it should be included so I can keep testing on the 3_zone_toy
    m.ts_previous = Param(
        m.TIMESERIES,
        within=m.TIMESERIES,
        initialize=lambda m, ts: m.TIMESERIES.prevw(ts))
    def Track_Storage_State_TP_rule(m, s, z, tp):

        #for the rare case of a period with a single timepoint
        if m.tp_previous[tp] == tp:
            return m.HydrogenStorageStateKg[s, z, tp] == (m.StoreHydrogenKgPerHour[s, z, tp] - m.WithdrawHydrogenKgPerHour[s, z, tp]) * m.tp_duration_hrs[tp]

        else:
            return m.HydrogenStorageStateKg[s, z, tp] ==  m.HydrogenStorageStateKg[s, z, m.tp_previous[tp]] + \
            m.HydrogenNetStorageKg[s, z, tp]
            #(m.StoreHydrogenKgPerHour[z, tp] - m.WithdrawHydrogenKgPerHour[z, tp]) * m.tp_duration_hrs[tp]

    m.Track_Storage_State_TP = Constraint(m.H2_STORAGE_TPS, rule=Track_Storage_State_TP_rule)

    def starting_tp_check(m, tp):
        currentTS = m.tp_ts[tp]
        currentPeriod = m.tp_period[tp]
        return (m.TPS_IN_TS[currentTS][1]==tp and m.TPS_IN_PERIOD[currentPeriod][1] != tp)

    m.IS_STARTING_TP = Set(initialize=m.TIMEPOINTS, filter=starting_tp_check)
    def Track_Storage_State_TS_rule(m, s, z, tp):
        currentTS = m.tp_ts[tp]
        previousTS = m.ts_previous[currentTS]
        finalPointOfSeries = m.TPS_IN_TS[previousTS][-1]
        return m.HydrogenStorageStateKg[s, z, tp] == m.HydrogenStorageStateKg[s, z, finalPointOfSeries] + \
        m.HydrogenNetStorageKg[s, z, tp]

    #At first glance this constraint appears to work, but the behavior of the model is MUCH different once this is in place, so I need
    #to look more closely before I can be sure. Did I check this? I think so? WHo knows anymore
    m.Track_Storage_State_TS = Constraint(m.H2_STORAGE_PROJECTS, m.LOAD_ZONES, m.IS_STARTING_TP, rule=Track_Storage_State_TS_rule)


    #going to try a slightly different way to track annual storage. If I constrained it such that the first and last timepoint of a period
    #must have a storagestate of 0, I think that might also solve this weird bug?
    #m.Hydrogen_Conservation_of_Mass_Annual = Constraint(m.LOAD_ZONES, m.PERIODS, rule=lambda m, z, p:
    #    m.HydrogenStorageStateKg[z, m.TPS_IN_PERIOD[p][1]] - m.HydrogenStorageStateKg[z, m.TPS_IN_PERIOD[p][-1]] == 0
    #)

    m.Hydrogen_Conservation_of_Mass_Annual_Start = Constraint(m.H2_STORAGE_BUILD_YRS_ZONES, rule=lambda m, s, z, p:
        m.HydrogenStorageStateKg[s, z, m.TPS_IN_PERIOD[p][1]] == 0
    )
    m.Hydrogen_Conservation_of_Mass_Annual_End = Constraint(m.H2_STORAGE_BUILD_YRS_ZONES, rule=lambda m, s, z, p:
        m.HydrogenStorageStateKg[s, z, m.TPS_IN_PERIOD[p][-1]] == 0
    )

    def Min_Storage_State_rule(m, s, z, tp):
        if m.tp_previous[tp] == tp:
            return m.HydrogenNetStorageKg[s, z,tp] >= 0
        else:
            return (m.HydrogenStorageStateKg[s, z, m.tp_previous[tp]] + m.HydrogenNetStorageKg[s, z, tp]) >=0
    m.Min_Storage_State = Constraint(m.H2_STORAGE_TPS, rule=Min_Storage_State_rule)

    m.Max_Storage_State = Constraint(m.H2_STORAGE_TPS, rule=lambda m, s, z, tp:
        m.HydrogenStorageStateKg[s, z, tp] <= m.H2StorageCapacityKg[s, z, m.tp_period[tp]])

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
    m.StoreHydrogenKgPerHourZonal = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        sum(m.StoreHydrogenKgPerHour[s, z, t] for s in m.H2_STORAGE_PROJECTS))
    m.WithdrawHydrogenKgPerHourZonal = Expression(m.LOAD_ZONES, m.TIMEPOINTS, rule=lambda m, z, t:
        sum(m.WithdrawHydrogenKgPerHour[s, z, t] for s in m.H2_STORAGE_PROJECTS))

    m.Hydrogen_Withdrawals.append('StoreHydrogenKgPerHourZonal')
    m.Hydrogen_Injections.append('WithdrawHydrogenKgPerHourZonal')





    """
    Code needed for H2 transport via pipelines
    Defining set of connectors, head and tail nodes, length, efficiency, capacity, and builds
    Copied over almost the entire transmission.build file since they are being modeled in the same way

    """

    m.TRANSPORT_PROJECTS = Set()
    m.TRANSPORT_ROUTES = Set()
    m.pipes_hn1 = Param(m.TRANSPORT_ROUTES, within=m.LOAD_ZONES)
    m.pipes_hn2 = Param(m.TRANSPORT_ROUTES, within=m.LOAD_ZONES)
    m.pipes_length_km = Param(m.TRANSPORT_ROUTES, within=NonNegativeReals)
    m.pipes_efficiency = Param(
        m.TRANSPORT_PROJECTS,
        within=PercentFraction)
    m.existing_pipes_cap = Param(
        m.TRANSPORT_ROUTES,
        within=NonNegativeReals)

    m.blend_limit = Param(
        m.TRANSPORT_PROJECTS,
        within=PercentFraction
    )
    m.is_pre_existing = Param(
        m.TRANSPORT_PROJECTS,
        within=Boolean
    )

    m.pipes_new_build_allowed = Param(
        m.TRANSPORT_ROUTES, within=Boolean, default=True)
    m.PIPELINES_BLD_YRS = Set(
        dimen=3,
        initialize=m.TRANSPORT_ROUTES * m.TRANSPORT_PROJECTS * m.PERIODS,
        filter=lambda m, r, p, period: m.pipes_new_build_allowed[r])
    m.PROJECT_YRS = Set(
        dimen=2,
        initialize=m.TRANSPORT_PROJECTS * m.PERIODS)
    m.BuildPtx = Var(m.PIPELINES_BLD_YRS, within=NonNegativeReals)
    m.pipes_terrain_multiplier = Param(
        m.TRANSPORT_ROUTES,
        within=NonNegativeReals,
        default=1)
    m.pipes_capital_cost_per_kg_km = Param(
        m.TRANSPORT_PROJECTS,
        within=NonNegativeReals,
        default=1000)
    m.pipes_lifetime_yrs = Param(
        m.TRANSPORT_PROJECTS,
        within=NonNegativeReals,
        default=20)
    m.pipes_fixed_om_fraction = Param(
        m.TRANSPORT_PROJECTS,
        within=NonNegativeReals,
        default=0.03)
    m.pipes_cost_annual = Param(
        m.TRANSPORT_ROUTES, m.TRANSPORT_PROJECTS,
        within=NonNegativeReals,
        initialize=lambda m, r, p: (
            m.pipes_capital_cost_per_kg_km[p] * m.pipes_terrain_multiplier[r] *
            m.pipes_length_km[r] * (crf(m.interest_rate, m.pipes_lifetime_yrs[p]) +
            m.pipes_fixed_om_fraction[p])))

    def PtxCapacityNameplate_rule(m, r, p, period):

        cap = sum(m.BuildPtx[r, p, bld_yr]
            for bld_yr in m.PERIODS
            if bld_yr <= period and (r, p, bld_yr) in m.PIPELINES_BLD_YRS)

        # Returning to this code 2 years later, is this double counting pre-existing builds??
        # No I don't think so since the pre existing capacity would not be counted in the BuildPtx variable
        if m.is_pre_existing[p]:
            cap += m.existing_pipes_cap[r]

        return cap

    m.PtxCapacityNameplate = Expression(
        m.PIPELINES_BLD_YRS,
        rule=PtxCapacityNameplate_rule)

    m.pipes_derating_factor = Param(
        m.TRANSPORT_ROUTES,
        within=PercentFraction,
        default=1)

    m.PtxCapacityNameplateAvailableForH2 = Expression(
        m.PIPELINES_BLD_YRS,
        rule=lambda m, r, p, period:
            m.PtxCapacityNameplate[r, p, period] * m.blend_limit[p]
        )

    def init_DIRECTIONAL_PTX(m):
        ptx_dir = set()
        for ptx in m.TRANSPORT_ROUTES:
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
        for ptx in m.TRANSPORT_ROUTES:
            if((m.pipes_hn1[ptx] == node_from and m.pipes_hn2[ptx] == node_to) or
               (m.pipes_hn2[ptx] == node_from and m.pipes_hn1[ptx] == node_to)):
                return ptx
    m.pipes_d_line = Param(
        m.DIRECTIONAL_PTX,
        within=m.TRANSPORT_ROUTES,
        initialize=init_pipes_d_line)

    m.PIPES_TIMEPOINTS = Set(
        dimen=4,
        initialize=lambda m: m.DIRECTIONAL_PTX * m.TRANSPORT_PROJECTS * m.TIMEPOINTS
    )
    m.DispatchPtx = Var(m.PIPES_TIMEPOINTS, within=NonNegativeReals)

    m.Maximum_DispatchPtx = Constraint(
        m.PIPES_TIMEPOINTS,
        rule=lambda m, node_from, node_to, p, tp: (
            m.DispatchPtx[node_from, node_to, p, tp] <=
            m.PtxCapacityNameplateAvailableForH2[m.pipes_d_line[node_from, node_to], p,
                                     m.tp_period[tp]]))

    m.PtxHydrogenSent = Expression(
        m.PIPES_TIMEPOINTS,
        rule=lambda m, node_from, node_to, p, tp: (
            m.DispatchPtx[node_from, node_to, p, tp]))

    m.PtxHydrogenReceived = Expression(
        m.PIPES_TIMEPOINTS,
        rule=lambda m, node_from, node_to, p, tp: (
            m.DispatchPtx[node_from, node_to, p, tp] *
            m.pipes_efficiency[p]))
            #m.PtxEffectiveEfficiency[m.pipes_d_line[node_from, node_to], m.tp_period[tp]]))
            #0.98))

    def PTXPowerNet_calculation(m, h, tp):
        return (sum(
                    sum(m.PtxHydrogenReceived[node_from, h, p, tp]
                        for node_from in m.PTX_CONNECTIONS_TO_NODE[h]) -
                    sum(m.PtxHydrogenSent[h, node_to, p, tp]
                        for node_to in m.PTX_CONNECTIONS_TO_NODE[h])
                for p in m.TRANSPORT_PROJECTS))

    m.PTXPowerNet = Expression(
        m.LOAD_ZONES, m.TIMEPOINTS,
        rule=PTXPowerNet_calculation)


    """
    Hydrogen pipeline line packing pseudocode

    1. Define new storage projects called pipeline
    2. Add optional input called IS_TX or IS_PIPELINE. Will default to False, but set to true for this project
    3. Set the cost of this storage project to 0, since the development of this project is already covered in transmission costs (unless there are other costs to consider?)
    4. Within the module, constrain the m.BuildH2StorageKg variable for all storage projects where IS_TX = True
    5. If true, force m.BuildH2StorageKg to be equal to half of the total pipeline capacity connected to this load zone, scaled based on allowable pressures and things like that
        a. Assuming 50-50 split between each load zone seems fine to me


    Look to see if there are any additional costs associated with linepacking and look to figure out a way to estimate storage capacity for a pipeline
    """

    m.IS_PX = Param(m.H2_STORAGE_PROJECTS, default=False, within=Boolean)

    m.LINEPACKING_PROJECTS = Set(initialize = m.H2_STORAGE_PROJECTS,
                            filter = lambda m, g : m.IS_PX[g])
    
    m.LINEPACKING_YRS_ZONES = Set(
        dimen=3,
        initialize = m.LINEPACKING_PROJECTS * m.LOAD_ZONES * m.PERIODS
    )

    m.LinePackingCapacity = Expression(m.LINEPACKING_YRS_ZONES, rule = lambda m,g,z,p:
        sum(m.PtxCapacityNameplateAvailableForH2[r, pr, p]*0.5 
        for (r, pr, per) in m.PIPELINES_BLD_YRS 
        if (per==p) and ((m.pipes_hn1[r]==z) or (m.pipes_hn2[r]==z)))
    )

    m.LinePackingConstraint = Constraint(m.LINEPACKING_YRS_ZONES, rule=lambda m,g,z,p:
        m.H2StorageCapacityKg[g,z,p] <= m.LinePackingCapacity[g,z,p])

    # Register net transmission as contributing to zonal energy balance
    m.Hydrogen_Injections.append('PTXPowerNet')

    m.PtxFixedCosts = Expression(
        m.PERIODS,
        rule=lambda m, period: sum(
            m.PtxCapacityNameplate[r, p, period] * m.pipes_cost_annual[r, p]
            for (r, p, per) in m.PIPELINES_BLD_YRS if per==period
        )
    )
    m.Cost_Components_Per_Period.append('PtxFixedCosts')

    # An expression to summarize annual costs for the objective
    # function. Units should be total annual future costs in $base_year
    # real dollars. The objective function will convert these to
    # base_year Net Present Value in $base_year real dollars.

    """
    Cost calculations for hydrogen infrastructure (not pipelines)
    """
    m.HydrogenVariableCost = Expression(m.TIMEPOINTS, rule=lambda m, t:
        sum(
            sum(m.DispatchH2Gen[g,z,t]*m.variable_cost_per_kg[g, m.tp_period[t]] for g in m.HYDROGEN_GEN)
            + sum(m.H2FuelUseRate[g, z, t] * m.fuel_cost[z, m.h2_gen_fuel[g], m.tp_period[t]] for g in m.NG_GEN)
            + m.DispatchFuelCellMW[z, t] * m.hydrogen_fuel_cell_variable_cost_per_mwh[m.tp_period[t]]
            for z in m.LOAD_ZONES
        )
    )

    def gen_build_can_operate_in_period(m, g, build_year, period):
        if build_year in m.PERIODS:
            online = m.period_start[build_year]
        else:
            online = build_year
        retirement = online + m.life_years[g]
        return (
            online <= m.period_start[period] < retirement
    )

    def conv_build_can_operate_in_period(m, g, build_year, period):
        if build_year in m.PERIODS:
            online = m.period_start[build_year]
        else:
            online = build_year
        retirement = online + m.hydrogen_conv_life_years[g]
        return (
            online <= m.period_start[period] < retirement
    )
    # This is probably more correct, but is a different behavior
    # mid_period = m.period_start[period] + 0.5 * m.period_length_years[period]
    # return online <= m.period_start[period] and mid_period <= retirement


    m.H2_BLD_YRS_FOR_GEN_PERIOD = Set(
    m.HYDROGEN_GEN, m.PERIODS,
    initialize=lambda m, g, period: set(
        bld_yr for (gen, bld_yr) in m.H2_GEN_PERIODS
        if gen == g and
            gen_build_can_operate_in_period(m, g, bld_yr, period)))
    
    m.H2_CONV_BLD_YRS_FOR_GEN_PERIOD = Set(
    m.H2_CONVERTERS, m.PERIODS,
    initialize=lambda m, g, period: set(
        bld_yr for (gen, bld_yr) in m.H2_CONV_PERIODS
        if gen == g and
            conv_build_can_operate_in_period(m, g, bld_yr, period)))

    #Review this later, because this might have to change from using capacities to using Build variables for everything -> Done
    m.HydrogenFixedCostAnnual = Expression(m.PERIODS, rule=lambda m, p:
        sum(
            sum(m.BuildH2Gen[g, z, p_] * m.capital_cost_per_kgh[g, p_] * crf(m.interest_rate, m.life_years[g]) for (g,p_) in m.H2_BLD_YRS_FOR_GEN_PERIOD)
            + sum(m.BuildH2Gen[g, z, p_] * m.fixed_cost_per_kgh_year[g, p_] for (g,p_) in m.H2_BLD_YRS_FOR_GEN_PERIOD)
            + sum(m.BuildH2StorageKg[s, z, p] * (
                m.h2_storage_capital_cost_per_kg[s, p] * crf(m.interest_rate, m.h2_storage_life_years[s]))
            + m.BuildH2StorageKgPower[s, z, p] * (
                m.h2_storage_capital_cost_per_kg_per_hour[s, p] * crf(m.interest_rate, m.h2_storage_life_years[s])) for s in m.H2_STORAGE_PROJECTS)
            + sum(m.BuildConverterMW[g, z, p] * (
                m.hydrogen_conv_capital_cost_per_mw[g, p] * crf(m.interest_rate, m.hydrogen_conv_life_years[g])
                + m.hydrogen_conv_fixed_cost_per_mw_year[g, p]) for g in m.H2_CONV_BLD_YRS_FOR_GEN_PERIOD)
            for z in m.LOAD_ZONES
        )
    )
    m.Cost_Components_Per_TP.append('HydrogenVariableCost')
    m.Cost_Components_Per_Period.append('HydrogenFixedCostAnnual')
    
    m.h2_credit_years = Set(
        dimen=2
    )
    m.h2_ptc_value = Param(
        m.h2_credit_years,
        default=0,
        domain=NonNegativeReals#,
        #doc="Production Tax Credit (PTC) for given technology by period. Data in $/MWh",
    )
    m.h2_itc_value = Param(
        m.h2_credit_years,
        default=0,
        domain=NonNegativeReals#,
        #doc="Investment Tax Credit (ITC) for given technology by period. Data in $/MW",
    )
    m.h2_carbon_capture_credit = Param(
        m.h2_credit_years,
        default=0,
        domain=NonNegativeReals
    )


    # Create a set that has build capacity constrained by year (both caps of the PTC).
    # The two caps of the PTC are that the generator must be built prior to 2035 in
    # order to recieve credit, and that all generators will no longer recieve the PTC
    # starting in 2040.
    m.h2_ptc_eligible_yrs = Set(
        m.HYDROGEN_GEN,
        m.PERIODS,
        ordered=False,
        initialize=lambda m, g, period: set(
            bld_yr
            for bld_yr in m.H2_BLD_YRS_FOR_GEN_PERIOD[g, period]
            if 2025 <= bld_yr <= 2035 and period <= 2040 #projects built in 2030 would be able to receive the credit in 2040, I think? Check this later. Yes!
        ),
    )
    # Calculate the total eligible PTC capacity per period
    m.H2_PTC_Capacity = Expression(
        m.HYDROGEN_GEN,
        m.PERIODS,
        rule=lambda m, g, period: sum(
            sum(m.BuildH2Gen[g, z, bld_yr] for z in m.LOAD_ZONES) for bld_yr in m.h2_ptc_eligible_yrs[g, period]
        ),
    )

    # Same as PTC_Capacity but per timepoint
    m.H2_PTC_CapacityInTP = Expression(
        m.HYDROGEN_GEN, m.TIMEPOINTS, rule=lambda m, g, t: m.H2_PTC_Capacity[g, m.tp_period[t]]
    )

    # Create PTC variable that will either return the PTC Capacity or the DispatchGen
    # whichever is minimum.
    m.H2PTC = Var(m.HYDROGEN_GEN, m.TIMEPOINTS, domain=NonNegativeReals)

    m.H2PTC_lower_bound = Constraint(
        m.HYDROGEN_GEN, m.TIMEPOINTS, rule=lambda m, g, t: m.H2PTC[g, t] <= m.H2_PTC_CapacityInTP[g, t]
    )
    m.H2PTC_upper_bound = Constraint(
        m.HYDROGEN_GEN, m.TIMEPOINTS, rule=lambda m, g, t: m.H2PTC[g, t] <= sum(m.DispatchH2Gen[g, z, t] for z in m.LOAD_ZONES)
    )

    # Calculate PTC
    m.H2_PTC_per_tp = Expression(
        m.TIMEPOINTS,
        rule=lambda m, t: sum(
            -m.H2PTC[g, t] * m.h2_ptc_value[m.tp_period[t], g]
            for g in m.HYDROGEN_GEN
            if g in set([item[1] for item in m.h2_credit_years.data()])
            and m.tp_period[t] <= 2045
        ),
    )
    m.Cost_Components_Per_TP.append("H2_PTC_per_tp")

    # Calculate Carbon Capture tax credit
    m.h2_ccs_eligible_yrs = Set(
        m.HYDROGEN_GEN,
        m.PERIODS,
        ordered=False,
        initialize=lambda m, g, period: set(
            bld_yr
            for bld_yr in m.H2_BLD_YRS_FOR_GEN_PERIOD[g, period]
            if 2025 <= bld_yr < 2035 and (period-bld_yr) < 12 
        ),
    )
    # Calculate the total eligible PTC capacity per period
    m.H2_CCS_Capacity = Expression(
        m.HYDROGEN_GEN,
        m.PERIODS,
        rule=lambda m, g, period: sum(
            sum(m.BuildH2Gen[g, z, bld_yr] for z in m.LOAD_ZONES) for bld_yr in m.h2_ccs_eligible_yrs[g, period]
        ),
    )

    # Same as PTC_Capacity but per timepoint
    m.H2_CCS_CapacityInTP = Expression(
        m.HYDROGEN_GEN, m.TIMEPOINTS, rule=lambda m, g, t: m.H2_CCS_Capacity[g, m.tp_period[t]]
    )

    # Create PTC variable that will either return the PTC Capacity or the DispatchGen
    # whichever is minimum.
    m.H2_CCS_credit = Var(m.HYDROGEN_GEN, m.TIMEPOINTS, domain=NonNegativeReals)

    m.H2_CCS_credit_lower_bound = Constraint(
        m.HYDROGEN_GEN, m.TIMEPOINTS, rule=lambda m, g, t: m.H2_CCS_credit[g, t] <= m.H2_CCS_CapacityInTP[g, t]
    )
    m.H2_CCS_credit_upper_bound = Constraint(
        m.HYDROGEN_GEN, m.TIMEPOINTS, rule=lambda m, g, t: m.H2_CCS_credit[g, t] <= sum(m.DispatchH2Gen[g, z, t] for z in m.LOAD_ZONES)
    )

    # Calculate 45Q PTC
    #Note that this calculation uses an estimated equivalent $/kg value for the 45Q credit based on potential capture rates for SMR, not a precise calculation based
    #on captured carbon. May update that later.
    m.H2_CCS_credit_per_tp = Expression(
        m.TIMEPOINTS,
        rule=lambda m, t: sum(
            -m.H2_CCS_credit[g, t] * m.h2_carbon_capture_credit[m.tp_period[t], g]
            for g in m.HYDROGEN_GEN
            if g in set([item[1] for item in m.h2_credit_years.data()])
            and m.tp_period[t] < 2045
        ),
    )
    m.Cost_Components_Per_TP.append("H2_CCS_credit_per_tp")

def load_inputs(m, switch_data, inputs_dir):

    """
    Need input files for: Hydrogen load zones, hydrogen time points, hydrogen projects (gen and storage), pipeline transport

    Do we want a local_td equivalent for hydrogen pipelines? Might be nice to include some representation of the distribution level of pipelines, even if it is simple

    Example syntax for loading input data below
    """
    #Can use this same file from original hydrogen module for technology information
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'hydrogen_generation_costs.csv'),
        optional=False,
        auto_select=True,
        index = m.H2_GEN_PERIODS,
        param=(
            m.capital_cost_per_kgh,
            m.fixed_cost_per_kgh_year,
            m.variable_cost_per_kg,
        )
    )

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'hydrogen_projects_costs.csv'),
        optional=False,
        auto_select=True,
        index = m.H2_CONV_PERIODS,
        param=(
            m.hydrogen_conv_capital_cost_per_mw,
            m.hydrogen_conv_fixed_cost_per_mw_year,
            m.hydrogen_conv_variable_cost_per_mwh,
        )
    )
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'hydrogen_generation_params.csv'),
        optional=False,
        auto_select=True,
        index=m.HYDROGEN_GEN,
        param=(
            m.kg_per_unit,
            m.life_years,
            m.co2_per_kg,
            m.h2_gen_fuel,
            m.IS_ELECTROLYSIS,
        )
    )
    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'hydrogen_projects_params.csv'),
        optional=False,
        auto_select=True,
        index=m.H2_CONVERTERS,
        param=(
            m.hydrogen_conv_life_years,
            m.hydrogen_conv_mwh_per_kg,
        )
    )

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'hydrogen_emissions.csv'),
        optional=False,
        #optional_params=(m.hydrogen_emission_cap),
        auto_select=True,
        index=m.PERIODS,
        param=(m.hydrogen_emission_cap)
    )

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'hydrogen_storage_costs.csv'),
        optional=False,
        auto_select=True, index=m.H2_STORAGE_BUILD_YRS,
        param=(
            m.h2_storage_capital_cost_per_kg,
            m.h2_storage_capital_cost_per_kg_per_hour,
            m.h2_storage_fixed_om_per_kg,
            m.h2_storage_fixed_om_per_kg_per_hour
        )
    )

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'hydrogen_storage_params.csv'),
        optional=False,
        auto_select=True, index=m.H2_STORAGE_PROJECTS,
        param=(
            m.h2_storage_life_years,
            m.h2_storage_minimum_size_kg,
            m.IS_PX
        )
    )

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
        index=m.TRANSPORT_ROUTES,
        optional_params=(
            'pipes_dbid', 'pipes_derating_factor',
            'pipes_terrain_multiplier', 'pipes_new_build_allowed'
        ),
        param=(
            m.pipes_hn1, m.pipes_hn2,
            m.pipes_length_km, m.existing_pipes_cap
        )
    )

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'pipes_params.csv'),
        optional=True, auto_select=True,
        index=m.TRANSPORT_PROJECTS,
        param=(
            m.pipes_capital_cost_per_kg_km, m.pipes_lifetime_yrs,
            m.pipes_fixed_om_fraction, m.blend_limit, m.pipes_efficiency, m.is_pre_existing
        )
    )

    switch_data.load_aug(
        filename=os.path.join(inputs_dir, 'h2_tax_credits.csv'),
        autoselect=True,
        index=m.h2_credit_years,
        param=(m.h2_ptc_value, m.h2_itc_value, m.h2_carbon_capture_credit))


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
        instance, instance.H2_STORAGE_PROJECTS*instance.LOAD_ZONES*instance.PERIODS,
        output_file=os.path.join(outdir, "hydrogen_storage_build.csv"),
        headings=("H2_STORAGE_PROJECTS","Load Zone", "period", "BuildH2StorageKg", "BuildH2StorageKgPower",
         "H2StorageCapacityKg", "H2StorageCapacityKgPerHour"),
        values=lambda m, s, z, p: (
            s, z, p, m.BuildH2StorageKg[s, z, p], m.BuildH2StorageKgPower[s, z, p],
            m.H2StorageCapacityKg[s, z, p], m.H2StorageCapacityKgPerHour[s, z, p]
            ))
    reporting.write_table(
        instance, instance.H2_STORAGE_PROJECTS*instance.LOAD_ZONES*instance.TIMEPOINTS,
        output_file=os.path.join(outdir, "hydrogen_storage_dispatch.csv"),
        headings=("H2_STORAGE_PROJECTS","Load Zone", "timepoint", "StoreHydrogenKgPerHour", "WithdrawHydrogenKgPerHour",
         "HydrogenStorageStateKg", "HydrogenNetStorageKg"),
        values=lambda m, s, z, tp: (
            s, z, tp, m.StoreHydrogenKgPerHour[s, z, tp], m.WithdrawHydrogenKgPerHour[s, z, tp],
            m.HydrogenStorageStateKg[s, z, tp], m.HydrogenNetStorageKg[s, z, tp]
            ))

    reporting.write_table(
        instance, instance.HYDROGEN_GEN*instance.LOAD_ZONES*instance.PERIODS,
        output_file=os.path.join(outdir, "hydrogen_cap.csv"),
        headings=("HYDROGEN_GEN","Load Zone", "period", "BuildH2Gen", "H2GenCapacity"),
        values=lambda m, g, z, p: (
            g, z, p, m.BuildH2Gen[g, z, p], m.H2GenCapacity[g, z, p]
            ))
    reporting.write_table(
        instance, instance.H2_CONVERTERS*instance.LOAD_ZONES*instance.PERIODS,
        output_file=os.path.join(outdir, "hydrogen_cap_converters.csv"),
        headings=("H2_CONVERTERS","Load Zone", "period", "ConverterCapacityMW"),
        values=lambda m, g, z, p: (
            g, z, p, m.ConverterCapacityMW[g, z, p]
            ))

    reporting.write_table(
        instance, instance.PERIODS,
        output_file=os.path.join(outdir, "h2_annual_emissions.csv"),
        headings=("period", "H2AnnualEmissions"),
        values=lambda m, p: (
            p, m.H2AnnualEmissions[p]
            ))

    reporting.write_table(
        instance, instance.HYDROGEN_GEN*instance.LOAD_ZONES*instance.TIMEPOINTS,
        output_file=os.path.join(outdir, "hydrogen_dispatch.csv"),
        headings=("HYDROGEN_GEN", "Load Zone", "timepoint", "DispatchH2Gen"),
        values=lambda m, g, z, t: (
            g, z, m.tp_timestamp[t],
             m.DispatchH2Gen[g, z, t]
            ))
    reporting.write_table(
        instance, instance.H2_CONVERTERS*instance.LOAD_ZONES*instance.TIMEPOINTS,
        output_file=os.path.join(outdir, "converter_dispatch.csv"),
        headings=("H2_CONVERTERS", "Load Zone", "timepoint", "DispatchConverterMW"),
        values=lambda m, g, z, t: (
            g, z, m.tp_timestamp[t],
             m.DispatchConverterMW[g, z, t]
            ))
    reporting.write_table(
        instance, instance.HNODE_TPS,
        output_file=os.path.join(outdir, "hydrogen_dispatch_other.csv"),
        headings=("Load Zone", "timepoint",
        "ConsumeElecMW", "ConsumeHydrogenKgPerHour", "NetTransmission"),
        values=lambda m, z, t: (
            z, m.tp_timestamp[t],
             m.ConsumeElecMW[z,t], m.ConsumeHydrogenKgPerHour[z,t], m.PTXPowerNet[z,t]
            ))

    reporting.write_table(
        instance, instance.PIPELINES_BLD_YRS,
        output_file=os.path.join(outdir, "PtxNameplateCapacity.csv"),
        headings=("TRANSPORT_ROUTE", "TRANSPORT_PROJECT", "period", "PtxNameplateCapacity", "PtxNameplateCapacityForH2"),
        values=lambda m, r, p, period: (
            r, p, period, m.PtxCapacityNameplate[r, p, period], m.PtxCapacityNameplateAvailableForH2[r, p, period]
            ))
