'''
Description: In childhood survival project, the script is used to generate a scratch table with any intermediate metrics to calculate 
            final results and final results table with excess deaths and total prevelance for a location/cause/year. The scratch table
            will be stored in /ihme/cancer/childhood_survival/ped_survival_scratch/cause_id/{location_id}_{year_id}.csv. The final 
            results table will be stored in /ihme/cancer/childhood_survival/results_table/cause_id/{location_id}_{year_id}.csv.

How to use: Run the following command on cluster in can_dev_env: python get_excess_deaths.py cause_id location_id year_id
Contributors: Weijia Fu
Date: 9/15/2021
'''
import numpy as np
import pandas as pd
import os
from sys import argv

from db_queries import get_outputs
from db_queries import get_envelope
from db_queries import get_ids
from db_queries import get_age_weights

def fill_in_year(table):
    '''Fill in year_id, year since diagnosis (year_since_dx) and follow up year( f/u year) '''
    table['year_id'] = year_id
    table['year_dx'] = table['year_id'] - table['year_since_dx']
    table['f/u_year'] = table['year_id'] + table['year_since_dx'] - 10
    return(table)

def fill_in_acmr(table, year_id, location_id):
    '''Fill in all cause mortality rate (ACMR) pulled from GBD mortality envelope'''
    acmr = get_envelope(year_id = year_id, sex_id= [1,2], gbd_round_id = 7, with_hiv=1, rates=1, location_id = location_id, age_group_id = 'all', decomp_step='step3')
    acmr = acmr[['age_group_id', 'year_id', 'sex_id', 'mean']]
    table = table.merge(acmr, on = ['age_group_id', 'year_id', 'sex_id'], how = 'left')
    table = table.rename({'mean': 'ACMR'}, axis = 1)
    table.loc[table['months_alive'] == 0, 'ACMR'] = 0
    return(table)

def fill_in_smr(table):
    '''Fill in standardized mortality ratio (SMR)'''
    conditions = [(table['year_since_dx'] < 40), (table['year_since_dx'] >= 40) & (table['year_since_dx'] < 60), table['year_since_dx'] >= 60]
    #conditions = [(table['age_group_id'] < 13), (table['age_group_id'] >= 13) & (table['age_group_id'] < 17), table['age_group_id'] >= 17]
    values = [5, 3, 2]
    #values = [5,5,5]
    #values = [5, 2, 2]
    table['SMR'] = np.select(conditions, values)
    return(table)

def get_cause_name(cause_id):
    '''Get cause name given cause_id from get_ids function'''
    cause_meta = get_ids('cause')
    acause = cause_meta[cause_meta['cause_id'] == cause_id]['acause'].to_string(index = False)
    return(acause)

def get_survival(acause, location_id, year_id):
    '''Pull absolute survival from GBD survival estimates compiled by compile_surv.py'''
    surv_year = pd.read_csv("/ihme/cancer/childhood_survival/abs_survival/{}_ped_abs_surv.csv".format(acause))
    surv_year = surv_year[['age_group_id', 'year_id', 'sex_id', 'location_id', 'surv_abs_mean']]
    surv_year_loc = surv_year.loc[(surv_year['location_id'] == location_id) & (surv_year['year_id'] == year_id)]
    surv_year_loc = surv_year_loc.rename({'surv_abs_mean': 'abs_surv'}, axis = 1)
    # generate survival under 5 
    surv_year_loc['age_group_id'].replace({2:1}, inplace = True)
    surv_year_loc = surv_year_loc.drop('year_id', axis = 1)
    return(surv_year_loc)

def generate_case_begin(cause_id, location_id, year_id, table):
    '''Generate number of cases after 10 years of diagnosis (case_beginning) based on GBD absolute survival and 
    incidence estimates.
    '''
    # Note: The test code used GBD 2019 incidence. Can change to gbd_round_id 7 when GBD 2020 results are finalized.
    can_inc_dx = get_outputs('cause', cause_id = cause_id, gbd_round_id = 6, decomp_step = 'step4', version='latest', measure_id=6, 
                            sex_id= [1,2], year_id = year_id - 10, location_id = locaiton_id, age_group_id = [1, 6, 7, 8], metric_id = 1)
    can_inc_dx = can_inc_dx[['age_group_id', 'sex_id', 'year_id', 'val']]
    can_inc_dx = can_inc_dx.rename({'val': 'inc'}, axis = 1)
    acause = get_cause_name(cause_id)
    surv_year_loc = get_survival(acause, location_id, year_id)  
    inc_table = can_inc_dx.merge(surv_year_loc, on = ['age_group_id', 'sex_id'])
    inc_table['case_beginning'] = inc_table['inc'] * inc_table['abs_surv']
    table.loc[(table['sex_id'] == 1) & (table['age_dx'] == ' 0-4') & (table['year_since_dx'] == 10), 'case_beginning'] = float(inc_table.loc[(inc_table['sex_id'] == 1) & (inc_table['age_group_id'] == 1)]['case_beginning'])
    table.loc[(table['sex_id'] == 1) & (table['age_dx'] == ' 5-9') & (table['year_since_dx'] == 10), 'case_beginning'] = float(inc_table.loc[(inc_table['sex_id'] == 1) & (inc_table['age_group_id'] == 6)]['case_beginning'])
    table.loc[(table['sex_id'] == 1) & (table['age_dx'] == ' 10-14') & (table['year_since_dx'] == 10), 'case_beginning'] = float(inc_table.loc[(inc_table['sex_id'] == 1) & (inc_table['age_group_id'] == 7)]['case_beginning'])
    table.loc[(table['sex_id'] == 1) & (table['age_dx'] == ' 15-19') & (table['year_since_dx'] == 10), 'case_beginning'] = float(inc_table.loc[(inc_table['sex_id'] == 1) & (inc_table['age_group_id'] == 8)]['case_beginning'])
    table.loc[(table['sex_id'] == 2) & (table['age_dx'] == ' 0-4') & (table['year_since_dx'] == 10), 'case_beginning'] = float(inc_table.loc[(inc_table['sex_id'] == 2) & (inc_table['age_group_id'] == 1)]['case_beginning'])
    table.loc[(table['sex_id'] == 2) & (table['age_dx'] == ' 5-9') & (table['year_since_dx'] == 10), 'case_beginning'] = float(inc_table.loc[(inc_table['sex_id'] == 2) & (inc_table['age_group_id'] == 6)]['case_beginning'])
    table.loc[(table['sex_id'] == 2) & (table['age_dx'] == ' 10-14') & (table['year_since_dx'] == 10), 'case_beginning'] = float(inc_table.loc[(inc_table['sex_id'] == 2) & (inc_table['age_group_id'] == 7)]['case_beginning'])
    table.loc[(table['sex_id'] == 2) & (table['age_dx'] == ' 15-19') & (table['year_since_dx'] == 10), 'case_beginning'] = float(inc_table.loc[(inc_table['sex_id'] == 2) & (inc_table['age_group_id'] == 8)]['case_beginning'])
    return(table)

def fill_in_death_case(table_sub):
    '''Fill in deaths and case_begining columns for one sub-table'''
    table_sub['deaths'] = np.nan
    index = table_sub.index
    for i in range(index[0], index[-1] + 1):
        if (i < index[-1]) :
            table_sub.loc[i,'deaths'] = table_sub.loc[i,'case_beginning'] * (table_sub.loc[i,'ACMR'] * (table_sub.loc[i, 'SMR'] - 1))
            table_sub.loc[i+1, 'case_beginning'] = table_sub.loc[i, 'case_beginning'] - table_sub.loc[i, 'deaths']
        else:
            table_sub.loc[i,'deaths'] = table_sub.loc[i,'case_beginning'] * (table_sub.loc[i,'ACMR'] * (table_sub.loc[i, 'SMR'] - 1))
    return(table_sub)

def concat_table_sub(table):
    '''Fill in deaths and case_beginning columns for the whole table '''
    outputs = []
    for age in table['age_dx'].unique():
        for sex in table['sex_id'].unique():
            table_sub = table[(table['sex_id'] == sex) & (table['age_dx'] == age)]
            table_sub_filled = fill_in_death_case(table_sub)
            outputs.append(table_sub_filled)
    outputs = pd.concat(outputs)
    return(outputs)

def fill_in_others(table):
    '''Fill in other columns '''
    table['case_remaining'] = table['case_beginning'] - table['deaths']
    table['person_months'] = table['deaths'] * table['months_alive']
    table['excess_deaths'] = table['deaths'] - table['case_beginning'] * table['ACMR'] 
    return(table)

def get_results_table(table, year_id):
    '''Generate results table including excess deaths and total prevalence'''
    results_table = pd.DataFrame(columns=['age_dx', 'sex_id', 'excess_deaths', 'prev_deaths_months', 'prev_alive_months', 'total_prev_months'])
    for age in table['age_dx'].unique():
        for sex in table['sex_id'].unique():
            table_sub = table[(table['sex_id'] == sex) & (table['age_dx'] == age)]
            deaths = table_sub['deaths'].sum()
            last_index = table_sub.index[-1]
            alive = table.loc[last_index, 'case_remaining']
            excess_deaths = table_sub['excess_deaths'].sum()
            prev_deaths_months = table_sub['person_months'].sum()
            prev_alive_months = alive * (table.loc[last_index, 'months_alive'] + 6)
            total_prev = prev_alive_months + prev_deaths_months
            df = {'age_dx': age, 'sex_id': sex, 'deaths': deaths, 'alive': alive, 'excess_deaths': excess_deaths, 'prev_deaths_months': prev_deaths_months, 'prev_alive_months': prev_alive_months,
                  'total_prev_months': total_prev, 'prev_deaths_year': prev_deaths_months/12, 'prev_alive_year': prev_alive_months/12, 'total_prev_year': total_prev/12}
            results_table = results_table.append(df, ignore_index= True)
    results_table['year_id'] = year_id
    return(results_table)

def get_results_all_ages(results_table, year_id):
    '''Generate final results for all ages '''
    results_table = results_table.groupby('sex_id', as_index=False).sum()
    results_table['age_dx'] = 'all ages'
    results_table['year_id'] = year_id
    return(results_table)

if __name__ == "__main__": 
    cause_id = int(argv[1])
    locaiton_id = int(argv[2])
    year_id = int(argv[3])
    table = pd.read_csv("/snfs1/temp/fuwj0528/childhood_survival/table1.csv")
    table = fill_in_year(table)
    table = fill_in_acmr(table, year_id, locaiton_id)
    table = fill_in_smr(table)
    table = generate_case_begin(cause_id, locaiton_id, year_id, table)
    table = concat_table_sub(table)
    table = fill_in_others(table)
    path_scratch = '/ihme/cancer/childhood_survival/ped_survival_scratch/{}/'.format(cause_id)
    if not os.path.exists(path_scratch):
        os.mkdir(path_scratch)
    table.to_csv('{}{}_{}.csv'.format(path_scratch, locaiton_id, year_id), index = False)
    results_table = get_results_table(table, year_id)
    results_table_all_ages = get_results_all_ages(results_table, year_id)
    results_table = pd.concat([results_table, results_table_all_ages])
    results_table['location_id'] = locaiton_id
    path_results = '/ihme/cancer/childhood_survival/results_table/{}/'.format(cause_id)
    if not os.path.exists(path_results):
        os.mkdir(path_results)
    results_table.to_csv('{}{}_{}.csv'.format(path_results, locaiton_id, year_id), index = False)
