import os
import pandas as pd
import numpy as np

"""
To-do list:
- Update the method for scaling to a period so that it pulls the corrects values from the actual csv files instead of assuming 10
- Clean up the I/O and see if we can streamline any of that
    - Also make sure the I/O works fine on a Linux system (it should?)
"""

#MGA parameters
iterations = 15          #Number of additional iterations you want to run following the cost-optimal solution
slack = .1                #The slack cost, which will the upper bound of the MGA decision space
"""
#to-do
Add in stopping criteria
"""
#Run the cost optimal solution
os.system("switch solve --verbose --sorted-output --outputs-dir optimal_outputs")

#Maybe change I/O syntax to match SWITCH syntax in reporting modules
#Append the mpa.py module to modules.txt file for the MGA iterations
file = open("inputs/modules.txt", "a")
file.write("\nmga")
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
        prev_output = "optimal_outputs/dispatch_annual_summary.csv"

        #breaking down the data from initial run into the components necessary for mga
        gen_tech = pd.read_csv(prev_output)
        gen_tech=gen_tech.rename(columns={'gen_tech' : 'TECH_TYPE', 'Energy_GWh_typical_yr' : 'prev_gen_tech'})

        #Scaling dispatch by number of years in period
        #I definitely don't need this since the tp_weight is already determined by the length of the period!
        #No, I do need this since I using the Dispatch for a typical year for the calculation, not the DispatchGen
        def dispatch_scaling(g):
            return g['prev_gen_tech'] * 10
        gen_tech['prev_gen_tech'] = gen_tech.apply(
            dispatch_scaling,#* p[ '2020' ][ 'period_length' ],
            axis=1,
            )
        #print(gen_tech)
        gen_tech = gen_tech.groupby("TECH_TYPE").sum()
        gen_tech = gen_tech.drop(["period", "VariableCost_per_yr", "DispatchEmissions_tCO2_per_typical_yr"], axis=1)

        coefficients = np.zeros(len(gen_tech))
        gen_tech['hsj_coefficients'] = coefficients

        #also send this over to the outputs for convenience
        gen_tech.to_csv("optimal_outputs/dispatch_coeff.csv")

    else:

        #post_solve in mga.py already sorted through the data so we just need to read it back in
        prev_output = output_folder + "/dispatch_coeff.csv"
        gen_tech = pd.read_csv(prev_output, index_col='TECH_TYPE')


    #move that data over to the input directory for next iteration

    gen_tech.to_csv("inputs/prev_run_gen.csv")
    print(gen_tech)

    #make new filename for new outputs
    output_folder = "mga_outputs_" + str(i+1)

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
