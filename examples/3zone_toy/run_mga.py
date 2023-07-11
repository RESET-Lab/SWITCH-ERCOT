import os
import pandas as pd
import numpy as np

"""
To-do list:
- Update the method for scaling to a period so that it pulls the corrects values from the actual csv files instead of assuming 10
- Clean up the I/O and see if we can streamline any of that
    - Also make sure the I/O works fine on a Linux system (it should?)
- Add an input to define the MGA method used
- Double check the math for the normalized method to make sure the MGA coefficient calculations are correct
"""

#MGA parameters
iterations = 2         #Number of additional iterations you want to run following the cost-optimal solution
slack = .10             #The slack cost, which will the upper bound of the MGA decision space
method = "integer"
"""
#to-do
Add in stopping criteria
"""
#Run the cost optimal solution
os.system("switch solve --verbose --sorted-output --outputs-dir optimal_outputs")

#Maybe change I/O syntax to match SWITCH syntax in reporting modules
#Append the mpa.py module to modules.txt file for the MGA iterations
file = open("inputs/modules.txt", "a")
file.write("\nswitch_model.mga")
file.close()

#Pull the optimal cost from the outputs of the first model to provide to next iterations
file = open("optimal_outputs/total_cost.txt", "r")
cost = file.readlines() #pull cost from text file
cost = cost[0][:-2] #remove \n from the string, leaving just the optimal cost
file.close()


#Send the slack and the optimal cost from first run to an input file for the MGA module
mga_input = pd.DataFrame({"mga_slack" : [slack], "opt_cost" : [cost]})
mga_input.to_csv("inputs/mga.csv", index=False)

#Read in periods.csv to get the number of years per period to scale the dispatch of typical year up to the dispatch of the full period
#period = pd.read_csv("inputs/periods.csv", index_col="INVESTMENT_PERIOD")
#period = period.rename_axis(index='period')
#Add new column with number of years in a period
#period['period_length'] = period.apply(
#    lambda p : p['period_end']-p['period_start'],
#    axis=1)
#print(period)

#For loop to run the MGA iterations
for i in range(iterations):

    #this will probably be different for the tutorial case than the real runs since the syntax is different
    if i == 0:

        distribution_loss = 0.053
        transmissionLines = pd.read_csv("inputs/transmission_lines.csv")
        genProjectsInfo = pd.read_csv("inputs/generation_projects_info.csv")
        dispatchGen = pd.read_csv("optimal_outputs/DispatchGen.csv")
        dispatchTx = pd.read_csv("optimal_outputs/DispatchTx.csv")
        timeSeries = pd.read_csv("inputs/timeseries.csv")
        timePoints = pd.read_csv("inputs/timepoints.csv")
        zones = pd.read_csv("inputs/load_zones.csv")

        #combine timescale information
        timeScales = timePoints.merge(timeSeries, left_on="timeseries", right_on="TIMESERIES")
        timeScales = timeScales.drop("TIMESERIES", axis=1)

        #merge genProjectsInfo and dispatchGen to get technology type and dispatch together
        intermediate = genProjectsInfo.merge(dispatchGen, left_on="GENERATION_PROJECT", right_on="GEN_TPS_1")
        dispatchTechs = intermediate.merge(timeScales, left_on="GEN_TPS_2", right_on="timepoint_id")
        dispatchTechs = dispatchTechs.drop(["gen_connect_cost_per_mw", "gen_capacity_limit_mw", "gen_full_load_heat_rate", "gen_variable_om", "gen_max_age", "gen_min_build_capacity", "gen_scheduled_outage_rate", "gen_forced_outage_rate",
                                            "gen_is_variable", "gen_is_baseload", "gen_is_cogen", "gen_unit_size", "GEN_TPS_1", "GEN_TPS_2", "timeseries", "gen_ccs_capture_efficiency", "gen_storage_efficiency", "gen_store_to_release_ratio"], axis=1)

        if "gen_ccs_energy_load" in dispatchTechs.columns:
            dispatchTechs["gen_ccs_energy_load"] = pd.to_numeric(dispatchTechs["gen_ccs_energy_load"].replace(".", 0))
        else:
            dispatchTechs["gen_ccs_energy_load"] = 0

        dispatchTechs["Total Dispatch (GWh)"] = dispatchTechs["DispatchGen"]*dispatchTechs["ts_duration_of_tp"]*dispatchTechs["ts_scale_to_period"]/1000
        dispatchTechs["Dispatch after CCS (GWh)"] = (dispatchTechs["DispatchGen"]-dispatchTechs["gen_ccs_energy_load"]*dispatchTechs["DispatchGen"])*dispatchTechs["ts_duration_of_tp"]*dispatchTechs["ts_scale_to_period"]/1000
        dispatchTechs.loc[dispatchTechs["Dispatch after CCS (GWh)"] < 0, "Dispatch after CCS (GWh)"] = 0

        tx_merged_1 = dispatchTx.merge(transmissionLines, left_on=["TRANS_TIMEPOINTS_1", "TRANS_TIMEPOINTS_2"], right_on=["trans_lz1", "trans_lz2"])
        tx_merged_2 = dispatchTx.merge(transmissionLines, left_on=["TRANS_TIMEPOINTS_2", "TRANS_TIMEPOINTS_1"], right_on=["trans_lz1", "trans_lz2"])
        tx_merged = tx_merged_1.append(tx_merged_2)

        tx_scaled = tx_merged.merge(timeScales, left_on="TRANS_TIMEPOINTS_3", right_on="timepoint_id")
        tx_scaled = tx_scaled.drop(["trans_lz1", "trans_lz2", "trans_length_km", "timepoint_id", "TRANSMISSION_LINE", "existing_trans_cap", "timeseries", "ts_num_tps"], axis=1)
        tx_scaled["tx_losses_GWh"] = tx_scaled["DispatchTx"]*(1-tx_scaled["trans_efficiency"])*tx_scaled["ts_duration_of_tp"]*tx_scaled["ts_scale_to_period"]/1000

        #Needed to calculate the total transmission losses associated with each load zone and timepoint
        tx_losses_zone_tp = pd.DataFrame()
        for index, row in zones.iterrows():
            for index2, row2 in timeScales.iterrows():
                #print(tx_scaled.loc[tx_scaled["TRANS_TIMEPOINTS_1"]==row["LOAD_ZONE"], ["TRANS_TIMEPOINTS_3","tx_losses_GWh"]])
                loss = tx_scaled.loc[(tx_scaled["TRANS_TIMEPOINTS_1"]==row["LOAD_ZONE"]) & (tx_scaled["TRANS_TIMEPOINTS_3"]==row2["timepoint_id"]), "tx_losses_GWh"].sum()
                zone_dispatch = dispatchTechs.loc[(dispatchTechs["gen_load_zone"]==row["LOAD_ZONE"]) & (dispatchTechs["timepoint_id"]==row2["timepoint_id"]), "Dispatch after CCS (GWh)"].sum()
                test = pd.Series([row["LOAD_ZONE"], row2["timepoint_id"], loss, zone_dispatch])
                test=pd.DataFrame(test)
                tx_losses_zone_tp=pd.concat([tx_losses_zone_tp, test.transpose()], axis=0)

        tx_losses_zone_tp = tx_losses_zone_tp.rename(columns={0:"Load Zone", 1:"Timepoint", 2:"TX Losses", 3:"Total Zonal Dispatch after CCS"})

        dispatchTechsTx = dispatchTechs.merge(tx_losses_zone_tp, left_on=["gen_load_zone", "timepoint_id"], right_on=["Load Zone", "Timepoint"])
        dispatchTechsTx["Proportion Dispatch"] = dispatchTechsTx["Dispatch after CCS (GWh)"] / dispatchTechsTx["Total Zonal Dispatch after CCS"]
        dispatchTechsTx["Dispatch after T&D Losses (GWh)"] = (1-distribution_loss)*(dispatchTechsTx["Dispatch after CCS (GWh)"] - dispatchTechsTx["Proportion Dispatch"]*dispatchTechsTx["TX Losses"])
        totalDispatchByTech = dispatchTechsTx.loc[:,["gen_tech", "Dispatch after T&D Losses (GWh)"]].groupby("gen_tech").sum()

        coefficients = np.zeros(len(totalDispatchByTech))
        totalDispatchByTech['hsj_coefficients'] = coefficients

        #also send this over to the outputs for convenience
        totalDispatchByTech.to_csv("optimal_outputs/dispatch_coeff.csv")

    else:

        #post_solve in mga.py already sorted through the data so we just need to read it back in
        prev_output = output_folder + "/dispatch_coeff.csv"
        totalDispatchByTech = pd.read_csv(prev_output, index_col='TECH_TYPE')


    #move that data over to the input directory for next iteration

    totalDispatchByTech.to_csv("inputs/prev_run_gen.csv")
    print(totalDispatchByTech)

    #make new filename for new outputs
    output_folder = "mga_" + method + "_" + str(i+1) + "_slack" + str(slack*100)

    #Running new model for each iteration, now with the MGA module and sending outputs to unique folders
    os.system("switch solve --verbose --sorted-output --outputs-dir " + output_folder)


#Move this to the front later and put it into an if statement to check if mga.py is in the list before running the optimal model

#Removes the mga.py module from modules.txt for future runs
file = open("inputs/modules.txt", "r")
lines = file.readlines()
lines = lines[:-1]
file.close()

file = open("inputs/modules.txt", "w")
for mod in lines:
    file.write(mod)
file.close()
