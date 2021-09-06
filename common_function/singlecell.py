import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import random
import glob
import os, sys
import math

#Global variables
fluorescence1 = 'TxRed'#'TxRed' or 'YFP'
fluorescence2 = 'YFP'
ratio = 'fluorescence_ratio'
dir_predel = '/Users/yuta_wakamotolab/Desktop/yuta_Wakamotolab/FpLrT19_20/FpLrT19d_Result/191218'
a1 = 'Mean_' + fluorescence1
b1 = 'Background_' + fluorescence1
a2 = 'Mean_' + fluorescence2
b2 = 'Background_' + fluorescence2
rate = 'elongation_rate'
plt.rcParams["font.size"] = 16


#Integrate result data
#revise and leave only nesessary data(Slice, Area, TxRed, YFP)
#210129 revise
def Single_cell_data_integrator(path):
    pre_list = ['YK0085','YK0138','YK0139']
    path_list = glob.glob(os.path.join(path, '**/ph'), recursive = True)
    for i in np.arange(len(path_list)):
        path_list[i] = path_list[i].rstrip('/ph')
        folder_list = os.listdir(path_list[i])
        folder_list.remove('ph')
        if '.DS_Store' in folder_list:
            folder_list.remove('.DS_Store')
        if 'integrated' in folder_list:
            folder_list.remove('integrated')
        else:
            os.mkdir(path_list[i] + '/integrated')


        Number = os.listdir(path_list[i] + '/ph')
        if '.DS_Store' in Number:
            Number.remove('.DS_Store')

        for j in np.arange(len(Number)):
            Number[j] = Number[j].lstrip('Results')
            Number[j] = Number[j].rstrip('.csv')

        for number in Number:
            Result_ph = pd.read_csv(path_list[i] + '/ph/Results'+str(number).zfill(4)+'.csv', encoding = 'shift_jis')
            Result_ph = Result_ph.rename(columns={'Mean':'Mean_ph','Background':'Background_ph'})
            
            for folder in folder_list:
                if folder == 'TxRed':
                    result_path = path_list[i] + '/' + folder + '/Results_Tx'
                elif folder == 'YFP':
                    result_path = path_list[i] + '/' + folder + '/Results_YFP'
                Result = pd.read_csv(result_path+str(number).zfill(4)+'.csv', encoding = 'shift_jis')

#                 Result_ph['Mean_' + folder] = Result['Mean']
#                 Result_ph['Background_' + folder] = Result['Background']
                Result_ph[folder] = Result['Mean'] - Result['Background']
            if any(x in path_list[i] for x in pre_list):
                number = number[0:2] + '3' + number[3]
            Result_ph['id_number'] = number[2]
            Result_ph['file'] = str(number).zfill(4)
            column_list = ['Slice', 'Area', 'file', 'id_number'] + folder_list
            Result_nessesary = Result_ph.loc[:, column_list]
            Result_nessesary.to_csv(path_list[i] + '/integrated/Results'+str(number).zfill(4)+'.csv', index = False)
            
##Single-cell analysis
def Read_results_data(Number, dr = './'):
    Result_all = pd.concat([pd.read_csv(dr+'/integrated/Results'+str(i).zfill(4)+'.csv', encoding = 'shift_jis') for i in Number])
        
    return Result_all.reset_index(drop = True)

def Calculate_normalize_term(Results, dr = './', BL_time = 0):
    if all([x in Results.columns.tolist() for x in ['TxRed', 'YFP']]):
        norm_Tx = np.mean(np.array(Results[Results['Slice'] < BL_time - 4][fluorescence1]))
        norm_YFP = np.mean(np.array(Results[Results['Slice'] < BL_time - 4][fluorescence2]))
        norm = norm_Tx/norm_YFP
    else:
        norm = 1
    
    return norm

def Calculate_single_normalize_term(Results, fluorescence1 = 'TxRed', dr = './', BL_time = 0):
    if all([x in Results.columns.tolist() for x in ['TxRed', 'YFP']]):
        norm = np.mean(np.array(Results[Results['Slice'] < BL_time - 4][fluorescence1]))
    else:
        norm = 1
    
    return norm


def All_Result_time_elongation_fluorescence_ratio(Result_all, id_number, directory='./', BL_time1 = 24, normalized_ratio = 1, win = 1, time_window = 6, norm_Tx = 1, norm_YFP = 1):
    
#     file_name_list = []
    time_list = []
    TxRed_list = []
    YFP_list = []
    ratio_list = []
    elongation_list =[]
    
    Result = Result_all[Result_all['id_number'] == id_number].copy()
    Result['Time'] = (Result['Slice']-BL_time1)/time_window
    Result[ratio] = Result[fluorescence1]/(Result[fluorescence2]*normalized_ratio)
    Result['Time2'] = Result['Time'].shift(1)
    Result['Area2'] = Result['Area'].shift(1)
    Result[rate] = np.log(Result['Area2']/Result['Area'])/(Result['Time2']-Result['Time'])
    Result.loc[Result[rate] <= np.log(0.75)/(Result['Time2']-Result['Time']), rate]=np.nan
    Result = Result[(Result['Time']<-0.5)|(Result['Time']>=0)]

    for t in np.arange(-8.0,73.0,1/time_window):
#210129 revise
#             file_name_list += [i] 
#         time_list += [t]
#         TxRed_list += [Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), fluorescence1].mean()]
#         YFP_list += [Result[(Result['Time']>=t-win)&(Result['Time']<t+win)][fluorescence2].mean()]
#         ratio_list += [Result[(Result['Time']>=t-win)&(Result['Time']<t+win)][ratio].mean()]
#         elongation_list += [Result[(Result['Time']>=t-win)&(Result['Time']<t+win)][rate].mean()]
        
        file_name_number = len(Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), 'file'].unique())
        time_list += [t]*file_name_number
        TxRed_list += Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), ['file',fluorescence1]].groupby('file').mean()[fluorescence1].tolist()
        YFP_list += Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), ['file',fluorescence2]].groupby('file').mean()[fluorescence2].tolist()
        ratio_list += Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), ['file',ratio]].groupby('file').mean()[ratio].tolist()
        elongation_list += Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), ['file',rate]].groupby('file').mean()[rate].tolist()

    Result_all = pd.DataFrame(data = {'Time':time_list,'fluorescence_ratio':ratio_list,'elongation_rate':elongation_list,'RplS-mCherry':TxRed_list,'RpsB-YFP':YFP_list})
    return Result_all.dropna(subset = [rate]).reset_index(drop=True)

def Median_95percent_area(Result_all, index_list = ['fluorescence_ratio', 'elongation_rate', 'RplS-mCherry', 'RpsB-YFP'], bootstrapping_number = 1000):
    #local lists
    time_list = []
    median_list = []
    bottom_list = []
    top_list =[]

    Result = Result_all.dropna()
    index_number_list = [Result.columns.get_loc(x) for x in index_list]
    for i in np.unique(Result['Time'].values):
        rng = np.random.default_rng()
        Result_time_i = Result[Result['Time'] == i]
        number = len(Result_time_i)
        Result_2d = Result_time_i.iloc[rng.choice(number, size = number*bootstrapping_number),index_number_list].to_numpy(copy = True)
        Result_3d = Result_2d.reshape(bootstrapping_number, number, len(index_number_list))
        Result_med_sort = np.sort(np.median(Result_3d, axis = 1), axis = 0)

        time_list += [i]
        bottom_list += [list(Result_med_sort[24])]
        top_list += [list(Result_med_sort[974])]
        median_list += [list(np.median(Result_time_i.iloc[:,index_number_list].to_numpy(copy = True), axis = 0))]


    time_array = np.array(time_list)
    bottom_array = np.array(bottom_list)
    median_array = np.array(median_list)
    top_array = np.array(top_list)

    df = pd.DataFrame({'Time':time_array})    
    for j in np.arange(len(index_list)):
        index = index_list[j]
        index_top = index + '_top'
        index_median = index + '_median'
        index_bottom = index + '_bottom'
        df_j = pd.DataFrame({index_top:top_array[:,j], index_median:median_array[:,j], index_bottom:bottom_array[:,j]})
        df = pd.concat([df, df_j], axis = 1)
        
    return df

def Figure_generator(Data_list, dr = '.', id_list = [], index_list = [], ylabel_list = [], color_list = [], label_list = [], max_list = [], min_list = [], w = 1/12, separate = True, xlim_min = -5, xlim_max = 75):
    plt.rcParams['text.usetex'] = True 
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}', r'\sansmath'] 
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Helvetica'

    plt.rcParams["font.size"] = 16
    
    if not os.path.exists(dr):
        os.mkdir(dr)
    for i in np.arange(len(index_list)):
        index = index_list[i]
        ylabel = ylabel_list[i]
        max_i = max_list[i]
        min_i = min_list[i]
        index_top = index + '_top'
        index_median = index + '_median'
        index_bottom = index + '_bottom'
        if separate == True:
            fig, ax = plt.subplots(figsize = (12,9))
            for data_index, id_index in enumerate(id_list):
                if len(Data_list) == 1:
                    Data = Data_list[0]
                else:
                    Data = Data_list[data_index]
                color = color_list[data_index]
                label = label_list[data_index]
                ax.plot(Data['Time'],Data[index_median], color = color, label = label)
                ax.fill_between(Data['Time'],Data[index_top],Data[index_bottom],color = color,alpha = 0.25,zorder = 7)
            ax.legend()
            ax.set_ylabel(ylabel)
            ax.set_xlabel('Time (hr)')
            ax.vlines(0,-5,max_i,zorder = 10)
            ax.hlines(0,xlim_min,xlim_max,zorder = 11)
            ax.set_ylim(-0.3,max_i)
            ax.set_xlim(xlim_min,xlim_max)
            ax.set_ylabel(ylabel)
            if w == 1/12:
                fig.savefig(dr + '/'+index+'.pdf', transparent=True, bbox_inches = 'tight')
            else:
                rolling_number ='_rolling_'+str(2*w).zfill(4)+'hr' 
                fig.savefig(dr + '/'+index+rolling_number+'.pdf', transparent=True, bbox_inches = 'tight')
            plt.clf()
            plt.close()
            
def calculate_max_mim_list(Data_list, index_list):
    def roundup(number, digit):
        if number >= 0:
            return (int(number/(10.0**(digit-1))+1))*(10.0**(digit-1))
        else:
            return (int(number/(10.0**(digit-1))))*(10.0**(digit-1))

    def rounddown(number, digit):
        if number >= 0:
            return (int(number/(10.0**(digit-1))))*(10.0**(digit-1))
        else:
            return (int(number/(10.0**(digit-1))-1))*(10.0**(digit-1))

    index_top_list = [x + '_top' for x in index_list]
    index_bottom_list = [x + '_bottom' for x in index_list]
    maximum_list = pd.concat([Data_list[i][index_top_list].max() for i in np.arange(len(Data_list))],axis = 1).T.max().tolist()
    minimum_list = pd.concat([Data_list[i][index_bottom_list].min() for i in np.arange(len(Data_list))],axis = 1).T.max().tolist()
    round_max_list = [int(math.log10(abs(x))) for x in maximum_list]
    round_min_list = [int(math.log10(abs(x))) for x in minimum_list]
    round_list = [max(round_max_list[i], round_min_list[i]) for i in np.arange(len(round_max_list))]
    return list(map(roundup, maximum_list, round_list)), list(map(rounddown, minimum_list, round_list))

def All_Result_time_elongation_rate(Result_all, id_number, directory='./', BL_time1 = 24.0, win = 1.0, time_window = 6.0):
    
    time_list = []
    TxRed_list = []
    ratio_list = []
    elongation_list =[]
    
    Result = Result_all[Result_all['id_number'] == id_number].copy()
    Result['Time'] = (Result['Slice']-BL_time1)/time_window
    Result['Time2'] = Result['Time'].shift(1)
    Result['Area2'] = Result['Area'].shift(1)
    Result[rate] = np.log(Result['Area2']/Result['Area'])/(Result['Time2']-Result['Time'])
    Result.loc[Result[rate] <= np.log(0.75)/(Result['Time2']-Result['Time']), rate]=np.nan
    Result = Result[(Result['Time']<-0.5)|(Result['Time']>=0)]

    for t in np.arange(-8.0, 73.0, 1.0/time_window):
        file_name_number = len(Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), 'file'].unique())
        time_list += [t]*file_name_number
        TxRed_list += Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), ['file',fluorescence1]].groupby('file').mean()[fluorescence1].tolist()
        elongation_list += Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), ['file',rate]].groupby('file').mean()[rate].tolist()

    Result_all = pd.DataFrame(data = {'Time':time_list,'elongation_rate':elongation_list,'mCherry-CAT':TxRed_list})
    return Result_all.dropna(subset = [rate]).reset_index(drop=True)  

def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None, maxasterix=None):
    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.01
        # ** is p < 0.001
        # *** is p < 0.0001
        # etc.
        text = ''
        p = .01

        while data < p:
            text += '$\ast$'
            p /= 10.
            if maxasterix and len(text) == maxasterix:
                break
        if len(text) == 0:
            text = 'n.s.'
    lx, ly = center[num1], height[num1]
    rx, ry = center[num2], height[num2]
    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]
    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)
    y = max(ly, ry) + dh
    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh)
    plt.plot(barx, bary, c='black')
    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs
    plt.text(*mid, text, **kwargs)
    
def All_Result_time_elongation_fluorescence_ratio_norm(Result_all, id_number, directory='./', BL_time1 = 24, normalized_ratio = 1, win = 1, time_window = 6, norm_Tx = 1, norm_YFP = 1):
    
#     file_name_list = []
    time_list = []
    TxRed_list = []
    YFP_list = []
    ratio_list = []
    elongation_list =[]
    
    Result = Result_all[Result_all['id_number'] == id_number].copy()
    Result['Time'] = (Result['Slice']-BL_time1)/time_window
    Result[ratio] = Result[fluorescence1]/(Result[fluorescence2]*normalized_ratio)
    Result['TxRed'] /= norm_Tx
    Result['YFP'] /= norm_YFP
    Result['Time2'] = Result['Time'].shift(1)
    Result['Area2'] = Result['Area'].shift(1)
    Result[rate] = np.log(Result['Area2']/Result['Area'])/(Result['Time2']-Result['Time'])
    Result.loc[Result[rate] <= np.log(0.75)/(Result['Time2']-Result['Time']), rate]=np.nan
    Result = Result[(Result['Time']<-0.5)|(Result['Time']>=0)]

    for t in np.arange(-8.0,73.0,1/time_window):
#210129 revise
#             file_name_list += [i] 
#         time_list += [t]
#         TxRed_list += [Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), fluorescence1].mean()]
#         YFP_list += [Result[(Result['Time']>=t-win)&(Result['Time']<t+win)][fluorescence2].mean()]
#         ratio_list += [Result[(Result['Time']>=t-win)&(Result['Time']<t+win)][ratio].mean()]
#         elongation_list += [Result[(Result['Time']>=t-win)&(Result['Time']<t+win)][rate].mean()]
        
        file_name_number = len(Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), 'file'].unique())
        time_list += [t]*file_name_number
        TxRed_list += Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), ['file',fluorescence1]].groupby('file').mean()[fluorescence1].tolist()
        YFP_list += Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), ['file',fluorescence2]].groupby('file').mean()[fluorescence2].tolist()
        ratio_list += Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), ['file',ratio]].groupby('file').mean()[ratio].tolist()
        elongation_list += Result.loc[(Result['Time']>=t-win)&(Result['Time']<t+win), ['file',rate]].groupby('file').mean()[rate].tolist()

    Result_all = pd.DataFrame(data = {'Time':time_list,'fluorescence_ratio':ratio_list,'elongation_rate':elongation_list,'RplS-mCherry':TxRed_list,'RpsB-YFP':YFP_list})
    return Result_all.dropna(subset = [rate]).reset_index(drop=True)

def Modify_generation_data(data, BL_time):
    def div_func(data):
        return pd.DataFrame({'Division_number': data.loc[:,'Division_number'] - data.loc[data.Slice == BL_time,'Division_number'].iloc[-1]})
    
    test_Result = data.copy()
    test_Result.loc[:,'Division']= test_Result.loc[:, 'Area']*0.75 > test_Result.loc[:,'Area'].shift(1)
    test_Result = test_Result[test_Result['Slice'].shift(1)!=1]
    test_Result.loc[:,'Division_number']=test_Result['Division'].sum()-test_Result['Division'].cumsum()
    test_Result.loc[:,'Division_number'] = test_Result.groupby('file').apply(div_func)
    return test_Result

def calculation_elongation_rate(data):
    initial_area = data.Area.iloc[-1]
    final_area = data.Area.iloc[0]
    time = (data.Slice.iloc[0]-data.Slice.iloc[-1])/6
    if time == 0:
        return np.nan
    else:
        return np.log(final_area/initial_area)/time
