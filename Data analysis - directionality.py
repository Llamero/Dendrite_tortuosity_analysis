#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Ben
#
# Created:     14/10/2020
# Copyright:   (c) Ben 2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import pandas
import numpy
from tkinter import Tk
from tkinter.filedialog import askopenfilename, asksaveasfilename #https://stackoverflow.com/questions/3579568/choosing-a-file-in-python-with-simple-dialog
import os
import re

import scipy
import xlrd #reminder to install updated library as it is needed for pandas Excel import
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    #Import ID and score files
    def importExcel(title_string):
        Tk().withdraw()
        path = askopenfilename(title=title_string, filetypes=[("Excel files", "*.xls *.xlsx *.xlsm *.xlsb *.odf *.ods *.odt *.csv"), ("CSV spreadsheet", "*.csv")]) #https://stackoverflow.com/questions/44403566/add-multiple-extensions-in-one-filetypes-mac-tkinter-filedialog-askopenfilenam
        extension = os.path.splitext(path)[1] #https://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
        if(extension in [".xls", ".xlsx", ".xlsm", ".xlsb", ".odf", ".ods", ".odt"]):
            df = pandas.read_excel(path)
        elif(extension == ".csv"):
            df = pandas.read_csv(path)
        else:
            print("Error: " + str(extension) + " is not a valid file extension.  How did you get here?")
            sys.exit()

        return (df, path)

    def concatenateDataframes(local_df, global_df, score_df):
        #Get sample name from file path
        score_df["sample"] = score_df["path"].str.extract(r"(?P<sample>[\\/][^\\/]*$)") #Use regex to get file name - https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.str.extract.html

        #remove file extensions from score sample name column
        score_df["sample"] = score_df["sample"].str.replace(r"\.[^.]+$", "") #https://stackoverflow.com/questions/22235518/regex-for-any-file-extension

        #remove file extensions from score sample name column
        score_df["sample"] = score_df["sample"].str[1:] #https://stackoverflow.com/questions/42349572/remove-first-x-number-of-characters-from-each-row-in-a-column-of-a-python-datafr

        #Normalize dF column names to lower case - https://cmdlinetips.com/2020/07/cleaning_up_pandas-column-names/
        local_df = local_df.rename(columns=str.lower)
        global_df = global_df.rename(columns=str.lower)
        score_df = score_df.rename(columns=str.lower)

        #sort all dataframes by "sample"
        local_df = local_df.sort_values(["sample"], ascending=(True))
        global_df = global_df.sort_values(["sample"], ascending=(True))
        score_df = score_df.sort_values(["sample"], ascending=(True))

        #add prefix to dataframe columns
        local_df = local_df.add_prefix("local_")
        global_df = global_df.add_prefix("global_")
        score_df = score_df.add_prefix("score_")

        #Horiztonally concatenate dataframes
        concat_df = pandas.concat([local_df, global_df, score_df], axis=1)

        #Copy scores and IDs to ID_df and verify that the scores match

#        concat_df.at["1", "local_sample"] = 0 #Test code to force ID columns to not match
#        concat_df.at["2", "global_sample"] = 0 #Test code to force ID columns to not match
#        concat_df.at["3", "score_sample"] = 0 #Test code to force ID columns to not match
        if(not (concat_df["local_sample"].equals(concat_df["global_sample"]) and concat_df["score_sample"].equals(concat_df["global_sample"]))):
            print("Error: The sample IDs don't match.")
            sys.exit()

        return concat_df
    def createPublicationPlot(df, main_title):

        #Get tortuosity as aspect ratio
        df.local_weighted = [1/x for x in df.local_weighted]
        df.global_weighted = [1/x for x in df.global_weighted]
        genotype_labels = df.score_genotype.unique();  # Get list of all unique genotypes - https://chrisalbon.com/python/data_wrangling/pandas_list_unique_values_in_column/
        genotype_df_list = []

        df = df[df.score_genotype.str.contains("C57")]
        print(df[df.score_location.str.contains("center")].local_weighted.median())
        print(df[df.score_location.str.contains("peri")].local_weighted.median())
        print(df[df.score_location.str.contains("center")].global_weighted.median())
        print(df[df.score_location.str.contains("peri")].global_weighted.median())
        print(scipy.stats.kstest(df[df.score_location.str.contains("center")].local_weighted, df[df.score_location.str.contains("peri")].local_weighted))
        print(scipy.stats.kstest(df[df.score_location.str.contains("center")].global_weighted, df[df.score_location.str.contains("peri")].global_weighted))
       
        #Sort data by genotype and location
        df = df.sort_values(["score_genotype", "score_location"], ascending=(True, True)) #Sort df by genotype and then location - https://www.kite.com/python/answers/how-to-sort-a-pandas-dataframe-by-multiple-columns-in-python
        genotype_labels = df.score_genotype.unique(); #Get list of all unique genotypes - https://chrisalbon.com/python/data_wrangling/pandas_list_unique_values_in_column/
        location_labels = df.score_location.unique();
        center_df = df[df.score_location.eq("center")] #Filter dataframe for only center scores - https://cmdlinetips.com/2018/02/how-to-subset-pandas-dataframe-based-on-values-of-a-column/
        peri_df = df[df.score_location.eq("peri")]

        #local directionality plots
        local_dir_fig, local_dir_axes = plt.subplots(1, 2, figsize=(16,8), sharey=True) #Plot graphs as subplots - https://dev.to/thalesbruno/subplotting-with-matplotlib-and-seaborn-5ei8
        fixed_ylim = (0, 1)
        plt.setp(local_dir_axes, ylim=fixed_ylim) #Keep the y-axis for all subplots fixed to the max score range - https://stackoverflow.com/questions/31006971/setting-the-same-axis-limits-for-all-subplots-in-matplotlib
        local_dir_fig.suptitle(main_title + "\nLocal Directionality Analysis")
        local_dir_loc_bplot = sns.boxplot(ax=local_dir_axes[0], y="local_weighted", x="score_location", data=df, notch=False, color="white")
        plt.setp(local_dir_loc_bplot.artists, edgecolor = 'k', facecolor='w') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        plt.setp(local_dir_loc_bplot.lines, color='k') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        local_dir_loc_bplot.set_title("Location") #Add title - https://stackoverflow.com/questions/42406233/how-to-add-title-to-seaborn-boxplot
        local_dir_loc_bplot = sns.swarmplot(ax=local_dir_axes[0], y="local_weighted", x="score_location", data=df, hue="score_genotype", s=3)

        #Global directionality plots
        global_dir_loc_bplot = sns.boxplot(ax=local_dir_axes[1], y="global_weighted", x="score_location", data=df, notch=False, color="white")
        plt.setp(global_dir_loc_bplot.artists, edgecolor = 'k', facecolor='w') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        plt.setp(global_dir_loc_bplot.lines, color='k') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        global_dir_loc_bplot.set_title("Location") #Add title - https://stackoverflow.com/questions/42406233/how-to-add-title-to-seaborn-boxplot
        global_dir_loc_bplot = sns.swarmplot(ax=local_dir_axes[1], y="global_weighted", x="score_location", data=df, hue="score_genotype", s=3)

        plt.show()



    def createScatterBoxPlot(df, main_title): #https://towardsdatascience.com/scattered-boxplots-graphing-experimental-results-with-matplotlib-seaborn-and-pandas-81f9fa8a1801

        #Sort data by genotype and location
        df = df.sort_values(["score_genotype", "score_location"], ascending=(True, True)) #Sort df by genotype and then location - https://www.kite.com/python/answers/how-to-sort-a-pandas-dataframe-by-multiple-columns-in-python
        genotype_labels = df.score_genotype.unique(); #Get list of all unique genotypes - https://chrisalbon.com/python/data_wrangling/pandas_list_unique_values_in_column/
        location_labels = df.score_location.unique();
        center_df = df[df.score_location.eq("center")] #Filter dataframe for only center scores - https://cmdlinetips.com/2018/02/how-to-subset-pandas-dataframe-based-on-values-of-a-column/
        peri_df = df[df.score_location.eq("peri")]

        #Control plots
        control_fig, control_axes = plt.subplots(2, 2, figsize=(15,15), sharey=True) #Plot graphs as subplots - https://dev.to/thalesbruno/subplotting-with-matplotlib-and-seaborn-5ei8
        control_fig.suptitle(main_title + "\nControl plots")

        local_channel_comp_bplot = sns.boxplot(ax=control_axes[0,0], y="local_weighted", x="score_channels", data=df, notch=True, color="white")
        plt.setp(local_channel_comp_bplot.artists, edgecolor = 'k', facecolor='w') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        plt.setp(local_channel_comp_bplot.lines, color='k') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        local_channel_comp_bplot.set_title("Local channel comparison") #Add title - https://stackoverflow.com/questions/42406233/how-to-add-title-to-seaborn-boxplot
        loc_bplot = sns.swarmplot(ax=control_axes[0,0], y="local_weighted", x="score_channels", data=df, color="black", s=3)

        global_channel_comp_bplot = sns.boxplot(ax=control_axes[0,1], y="global_weighted", x="score_channels", data=df, notch=True, color="white")
        plt.setp(global_channel_comp_bplot.artists, edgecolor = 'k', facecolor='w') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        plt.setp(global_channel_comp_bplot.lines, color='k') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        global_channel_comp_bplot.set_title("global channel comparison") #Add title - https://stackoverflow.com/questions/42406233/how-to-add-title-to-seaborn-boxplot
        loc_bplot = sns.swarmplot(ax=control_axes[0,1], y="global_weighted", x="score_channels", data=df, color="black", s=3)

        local_dir_v_dens = sns.regplot(ax=control_axes[1,0], y="local_weighted", x="local_mask_density", data=df, ci=95)
        local_dir_v_dens.set_title("Local directionality vs mask density") #Add title - https://stackoverflow.com/questions/42406233/how-to-add-title-to-seaborn-boxplot

        global_dir_v_dens = sns.regplot(ax=control_axes[1,1], y="global_weighted", x="global_mask_density", data=df, ci=95)
        global_dir_v_dens.set_title("Global directionality vs mask density") #Add title - https://stackoverflow.com/questions/42406233/how-to-add-title-to-seaborn-boxplot

        #Mask density "score" plots
        mask_density_fig, mask_density_axes = plt.subplots(2, 2, figsize=(15,15), sharey=True) #Plot graphs as subplots - https://dev.to/thalesbruno/subplotting-with-matplotlib-and-seaborn-5ei8
        fixed_ylim = (df["local_mask_density"].min(), df["local_mask_density"].max())
        plt.setp(mask_density_axes, ylim=fixed_ylim) #Keep the y-axis for all subplots fixed to the max score range - https://stackoverflow.com/questions/31006971/setting-the-same-axis-limits-for-all-subplots-in-matplotlib
        mask_density_fig.suptitle(main_title + "\nMask Density Analysis")
        mask_density_loc_bplot = sns.boxplot(ax=mask_density_axes[0,0], y="local_mask_density", x="score_location", data=df, notch=True, color="white")
        plt.setp(loc_bplot.artists, edgecolor = 'k', facecolor='w') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        plt.setp(loc_bplot.lines, color='k') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        loc_bplot.set_title("Location") #Add title - https://stackoverflow.com/questions/42406233/how-to-add-title-to-seaborn-boxplot
        loc_bplot = sns.swarmplot(ax=mask_density_axes[0,0], y="local_mask_density", x="score_location", data=df, color="black", s=3)
        geno_bplot = sns.boxplot(ax=mask_density_axes[0,1], y="local_mask_density", x="score_genotype", data=df, notch=True, color="white")
        plt.setp(geno_bplot.artists, edgecolor = 'k', facecolor='w')
        plt.setp(geno_bplot.lines, color='k')
        geno_bplot.set_title("Genotpye")
        geno_bplot = sns.swarmplot(ax=mask_density_axes[0,1], y="local_mask_density", x="score_genotype", data=df, color="black", s=3)
        center_bplot = sns.boxplot(ax=mask_density_axes[1,0], y="local_mask_density", x="score_genotype", data=center_df, notch=True, color="white")
        plt.setp(center_bplot.artists, edgecolor = 'k', facecolor='w')
        plt.setp(center_bplot.lines, color='k')
        center_bplot.set_title("Center - Genotpye")
        center_bplot = sns.swarmplot(ax=mask_density_axes[1,0], y="local_mask_density", x="score_genotype", data=center_df, color="black", s=3)
        peri_bplot = sns.boxplot(ax=mask_density_axes[1,1], y="local_mask_density", x="score_genotype", data=peri_df, notch=True, color="white")
        plt.setp(peri_bplot.artists, edgecolor = 'k', facecolor='w')
        plt.setp(peri_bplot.lines, color='k')
        peri_bplot.set_title("Peri - Genotpye")
        peri_bplot = sns.swarmplot(ax=mask_density_axes[1,1], y="local_mask_density", x="score_genotype", data=peri_df, color="black", s=3)

        #Correlation between scores and density
        scores_vs_density_fig, scores_vs_density_axes = plt.subplots(1, 2, figsize=(20,10), sharey=True, sharex=True)
        scores_vs_density_fig.suptitle(main_title + "\nScore/Mask Density Correlation")
        sns.swarmplot(ax=scores_vs_density_axes[0], y="local_mask_density", x="score_score", data=df, hue="score_genotype")
        sns.swarmplot(ax=scores_vs_density_axes[1], y="local_mask_density", x="score_score", data=df, hue="score_location")

        #Correlation local and global directionality
        local_vs_global_fig, local_vs_global_axes = plt.subplots(1, 2, figsize=(20,10), sharey=True)
        local_vs_global_fig.suptitle(main_title + "\nLocal/Global Directionality Correlation")
        fixed_ylim = (df[["local_weighted", "global_weighted"]].min().min(), df[["local_weighted", "global_weighted"]].max().max())
        plt.setp(local_vs_global_axes, ylim=fixed_ylim, xlim=fixed_ylim)
        sns.scatterplot(ax=local_vs_global_axes[0], y="local_weighted", x="global_weighted", data=df, hue="score_genotype")
        sns.scatterplot(ax=local_vs_global_axes[1], y="local_weighted", x="global_weighted", data=df, hue="score_location")

        #local directionality plots
        local_dir_fig, local_dir_axes = plt.subplots(2, 2, figsize=(15,15), sharey=True) #Plot graphs as subplots - https://dev.to/thalesbruno/subplotting-with-matplotlib-and-seaborn-5ei8
        fixed_ylim = (df["local_weighted"].min(), df["local_weighted"].max())
        plt.setp(local_dir_axes, ylim=fixed_ylim) #Keep the y-axis for all subplots fixed to the max score range - https://stackoverflow.com/questions/31006971/setting-the-same-axis-limits-for-all-subplots-in-matplotlib
        local_dir_fig.suptitle(main_title + "\nLocal Directionlity Analysis")
        local_dir_loc_bplot = sns.boxplot(ax=local_dir_axes[0,0], y="local_weighted", x="score_location", data=df, notch=True, color="white")
        plt.setp(loc_bplot.artists, edgecolor = 'k', facecolor='w') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        plt.setp(loc_bplot.lines, color='k') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        loc_bplot.set_title("Location") #Add title - https://stackoverflow.com/questions/42406233/how-to-add-title-to-seaborn-boxplot
        loc_bplot = sns.swarmplot(ax=local_dir_axes[0,0], y="local_weighted", x="score_location", data=df, color="black", s=3)
        geno_bplot = sns.boxplot(ax=local_dir_axes[0,1], y="local_weighted", x="score_genotype", data=df, notch=True, color="white")
        plt.setp(geno_bplot.artists, edgecolor = 'k', facecolor='w')
        plt.setp(geno_bplot.lines, color='k')
        geno_bplot.set_title("Genotpye")
        geno_bplot = sns.swarmplot(ax=local_dir_axes[0,1], y="local_weighted", x="score_genotype", data=df, color="black", s=3)
        center_bplot = sns.boxplot(ax=local_dir_axes[1,0], y="local_weighted", x="score_genotype", data=center_df, notch=True, color="white")
        plt.setp(center_bplot.artists, edgecolor = 'k', facecolor='w')
        plt.setp(center_bplot.lines, color='k')
        center_bplot.set_title("Center - Genotpye")
        center_bplot = sns.swarmplot(ax=local_dir_axes[1,0], y="local_weighted", x="score_genotype", data=center_df, color="black", s=3)
        peri_bplot = sns.boxplot(ax=local_dir_axes[1,1], y="local_weighted", x="score_genotype", data=peri_df, notch=True, color="white")
        plt.setp(peri_bplot.artists, edgecolor = 'k', facecolor='w')
        plt.setp(peri_bplot.lines, color='k')
        peri_bplot.set_title("Peri - Genotpye")
        peri_bplot = sns.swarmplot(ax=local_dir_axes[1,1], y="local_weighted", x="score_genotype", data=peri_df, color="black", s=3)

        #Global directionality plots
        global_dir_fig, global_dir_axes = plt.subplots(2, 2, figsize=(15,15), sharey=True) #Plot graphs as subplots - https://dev.to/thalesbruno/subplotting-with-matplotlib-and-seaborn-5ei8
        fixed_ylim = (df["global_weighted"].min(), df["global_weighted"].max())
        plt.setp(global_dir_axes, ylim=fixed_ylim) #Keep the y-axis for all subplots fixed to the max score range - https://stackoverflow.com/questions/31006971/setting-the-same-axis-limits-for-all-subplots-in-matplotlib
        global_dir_fig.suptitle(main_title + "\nGlobal Directionlity Analysis")
        global_dir_loc_bplot = sns.boxplot(ax=global_dir_axes[0,0], y="global_weighted", x="score_location", data=df, notch=True, color="white")
        plt.setp(loc_bplot.artists, edgecolor = 'k', facecolor='w') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        plt.setp(loc_bplot.lines, color='k') #Draw boxplot lines as black - https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
        loc_bplot.set_title("Location") #Add title - https://stackoverflow.com/questions/42406233/how-to-add-title-to-seaborn-boxplot
        loc_bplot = sns.swarmplot(ax=global_dir_axes[0,0], y="global_weighted", x="score_location", data=df, color="black", s=3)
        geno_bplot = sns.boxplot(ax=global_dir_axes[0,1], y="global_weighted", x="score_genotype", data=df, notch=True, color="white")
        plt.setp(geno_bplot.artists, edgecolor = 'k', facecolor='w')
        plt.setp(geno_bplot.lines, color='k')
        geno_bplot.set_title("Genotpye")
        geno_bplot = sns.swarmplot(ax=global_dir_axes[0,1], y="global_weighted", x="score_genotype", data=df, color="black", s=3)
        center_bplot = sns.boxplot(ax=global_dir_axes[1,0], y="global_weighted", x="score_genotype", data=center_df, notch=True, color="white")
        plt.setp(center_bplot.artists, edgecolor = 'k', facecolor='w')
        plt.setp(center_bplot.lines, color='k')
        center_bplot.set_title("Center - Genotpye")
        center_bplot = sns.swarmplot(ax=global_dir_axes[1,0], y="global_weighted", x="score_genotype", data=center_df, color="black", s=3)
        peri_bplot = sns.boxplot(ax=global_dir_axes[1,1], y="global_weighted", x="score_genotype", data=peri_df, notch=True, color="white")
        plt.setp(peri_bplot.artists, edgecolor = 'k', facecolor='w')
        plt.setp(peri_bplot.lines, color='k')
        peri_bplot.set_title("Peri - Genotpye")
        peri_bplot = sns.swarmplot(ax=global_dir_axes[1,1], y="global_weighted", x="score_genotype", data=peri_df, color="black", s=3)





        plt.show()
#        plt.draw()
#        plt.pause(0.1)

    global_df, dummy = importExcel("Select global directionality file")
    local_df, dummy = importExcel("Select local directionality file")
    score_df, input_path = importExcel("Select score file")
    concat_df = concatenateDataframes(local_df, global_df, score_df)

    score_file_name = re.findall("[\\/][^\\/]*$", input_path)[0][1:] #Get file name
    score_file_name = re.sub(r"\.[A-Za-z0-9_]+$", "", score_file_name) #remove extension - https://stackoverflow.com/questions/11475885/python-replace-regex/11475905
    combined_df = concatenateDataframes(local_df, global_df, score_df) #Combine dataframes

#    plt.ion()
    plt.show()
#    createScatterBoxPlot(combined_df, score_file_name)
    createPublicationPlot(combined_df, score_file_name)

    #Save parsed csv file
#    output_path = asksaveasfilename(title="Save spreadsheet", filetypes=[("CSV spreadsheet", "*.csv")], defaultextension=".csv", initialfile="test" + " - sorted scores") #http://effbot.org/tkinterbook/tkinter-file-dialogs.htm
#    concat_df.to_csv(output_path, index=False, header=True) #https://datatofish.com/export-dataframe-to-csv/

##    #Save plot image
#    output_path = asksaveasfilename(title="Save plots", filetypes=[("PNG image", "*.png")], defaultextension=".png", initialfile=score_file_name + " - plots") #http://effbot.org/tkinterbook/tkinter-file-dialogs.htm
#    plt.savefig(output_path) #https://datatofish.com/export-dataframe-to-csv/
###    input("Press [enter] to continue.")

if __name__ == '__main__':
    main()
