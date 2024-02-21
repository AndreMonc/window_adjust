# !/usr/bin/env python
# encoding: utf-8

"""
Goal: To convert variable-length windows of mean recombination rate into 10kb 
windows.

Copyright 2024 Andre E. Moncrieff. All rights reserved.

"""

import argparse
from distutils.command import clean
import pandas
from collections import Counter
import numpy

pandas.options.mode.chained_assignment = None  # default='warn'


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataFile", required=True,
                        help="Enter the file name + extension",
                        type=str)
    parser.add_argument("--windowFile", required=True,
                        help="Enter the file name (+.bed)",
                        type=str)
    args = parser.parse_args()
    return args


def read_in_data(txt_file_dataframe):
    raw_dataframe = pandas.read_csv(txt_file_dataframe,
                                    sep='\t',
                                    encoding = "ISO-8859-1",
                                    dtype=str)
    return raw_dataframe



def read_in_windows(windowFile, column_names):
    raw_dataframe = pandas.read_csv(windowFile,
                                    sep='\t',
                                    encoding = "ISO-8859-1",
                                    dtype=str,
                                    names=column_names)
    return raw_dataframe


def window_adj(data_df, window_df):
    # window bed file columns to lists
    chrom_window = window_df['chrom'].tolist()  
    start_window = window_df['start'].tolist()
    end_window = window_df['end'].tolist()
    # ReLERNN output file columns to lists
    chrom_data = data_df['chrom'].tolist()
    start_data = data_df['start'].tolist()
    end_data = data_df['end'].tolist()
    nSites_data = data_df['nSites'].tolist()
    recombRate_data = data_df['recombRate'].tolist()
    # list to store final recombination values for 10kb windows in window file
    final_recomb_values = []
    for line in range(len(chrom_window)):
        # iterate over window file
        temp_list = [] # may store a recombination value or 'NA' for the 10kb 
        # window under several circumstances
        overlap_list = [] # stores recombination values for the 10kb window 
        # when this window overlaps two recombination windows   
        x = range(int(start_window[line]), int(end_window[line]))  
        for i in range(len(chrom_data)):
            #iterate over the ReLERNN recombination windows
            y = range(int(start_data[i]), int(end_data[i]))  
            window_overlap = range(max(x[0], y[0]), min(x[-1], y[-1])+1)
            if chrom_window[line] == chrom_data[i] and len(window_overlap) > 0:
                # within this loop, we know that the windows overlap
                if int(nSites_data[i]) < 200:
                    temp_list.append('NA')
                    break
                elif len(window_overlap) == 10000:
                    # if pre-defined window entirely within recomb rate window
                    temp_list.append(float(recombRate_data[i]))
                    break
                elif len(window_overlap) > 0 and len(window_overlap) < 10000:
                    # if pre-defined window not entireley within recomb rate
                    # window. This could only be due to pre-defined window
                    # overlapping two recombination windows. (Recomb windows 
                    # with 200+ nSites are all over 10kb in my datasets).
                    rec_rate_wind = len(window_overlap) * float(
                        recombRate_data[i])
                    overlap_list.append(float(rec_rate_wind))
                    if len(overlap_list) < 2:
                        continue
                    elif len(overlap_list) == 2:   
                        window_length = int(end_window[line]) - int(start_window[line])
                        # in almost all cases window length will be 10kb, exception is the 
                        # end of a chromosome
                        overlap_list = [float(i) for i in overlap_list]
                        sum_of_recomb_values = sum(overlap_list)
                        average_window_recomb_rate = sum_of_recomb_values/window_length
                        temp_list.append(float(average_window_recomb_rate))
                        break
            else:
                pass
        if len(temp_list) > 0:
            final_recomb_values.append(temp_list[0])
        if len(temp_list) == 0:
            final_recomb_values.append('NA')
    return final_recomb_values


def main():
    args = parser()    
    data_df = read_in_data(args.dataFile)    
    column_names = ['chrom', 'start', 'end']
    window_df = read_in_windows(args.windowFile, column_names) 
    final_recomb_values = window_adj(data_df, window_df)
    window_df['recombination_rate'] = final_recomb_values  
    window_df.to_csv('recomb_rate_by_10kb_windows.csv', sep=',', index=False)
    

if __name__ == '__main__':
    main()
