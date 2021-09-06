import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#Distribute MIC day data to strain folder
def MIC_data_distribution(main_path, branch_path):
    data_path = main_path + '/raw_data' + branch_path + '/day_data'
    Day_list = os.listdir(data_path)
    if '.DS_Store' in Day_list:
        Day_list.remove('.DS_Store')
    
    for day in np.arange(len(Day_list)):
        Day_list[day] = Day_list[day].rstrip('.csv')

    for day in np.arange(len(Day_list)):
        MIC_data = pd.read_csv(data_path + '/' + Day_list[day] + '.csv', encoding = 'shift-jis')
        strain_list = list(MIC_data.columns)
        strain_list.remove('Cp_con')
        for strain in np.arange(len(strain_list)):
            MIC_day_data = MIC_data.loc[:,['Cp_con',strain_list[strain]]].dropna(subset = [strain_list[strain]])
            save_path = main_path + '/raw_data' + branch_path + '/' + strain_list[strain] + '/' + Day_list[day] + '.csv'
            MIC_day_data.to_csv(save_path, index = False)
            
def MIC_figure(main_path, branch_path):
    folder_path = main_path + '/raw_data' + branch_path
    save_path = main_path + '/figure' + branch_path

    Folder_list=os.listdir(folder_path)
    Folder_list.remove('day_data')
    if '.DS_Store' in Folder_list:
        Folder_list.remove('.DS_Store')

    plt.rcParams['font.size'] = 16
    min_func = lambda x: max(x, 0.01)


    for strain in Folder_list:
        data_path = folder_path + '/' + strain
        Day_list = os.listdir(data_path)
        strain_total = pd.DataFrame({'Cp_con':[], strain:[]})
        for day in Day_list:
            day_data = pd.read_csv(data_path + '/' + day)
            strain_total = pd.concat([strain_total, day_data])
        strain_average = strain_total.groupby('Cp_con', as_index = False).mean()
        strain_total[strain] = strain_total[strain].map(min_func)

        fig, ax = plt.subplots(figsize = (12,8))
        ax.set_xlabel('Cp concentration (µg/mL)')
        ax.set_ylabel('OD$_{600}$')
        ax.set_xlim(0.1,1000)
        ax.set_ylim(0.01,2)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.plot(strain_average['Cp_con'], strain_average[strain], color = 'black', label = strain)
        ax.legend()
        ax.plot(strain_total['Cp_con'], strain_total[strain], marker = 'o', lw = 0, markersize = 3, color = 'grey')
        ax.hlines(0.01,0.1,1000, linestyle = '--')
        plt.savefig(save_path + '/' + strain + '.pdf', transparent = True, bbox_inches = 'tight')
        plt.close()