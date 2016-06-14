#-*- coding: utf-8 -*-
# MS2_Hunter Ver beta 0.2
#
# This script reads the .mzML MS data, and extracts the MS/MS information of SLGGDSIMGIQLVSR peptide
# and its modified analogue according to the tag ions
#
# Tag ions include y9, y8, y7, y5, and y4
#
# The result writes in a .csv file, which has the same name to .mzML file
# The recoding information include spectrumID, precursor's retention time, precursor's charge, precursor's m/z and
# the 5 highest peaks of MS/MS spectrum
#
# Based on hasPeak.py from example script of pymzML
# Author: Bo Pang (SIOC)
# Bless!
# Update V0.21- Based on V 0.2, adding the FileList Function to read the filename in current directory, return a list

import sys
import pymzml
import csv
import glob

def FileList():
    allfiles = glob.glob(r'*.mzML')
    return allfiles

def open_ms_file(ms_file_main):
    try:
        ms_file_name = ms_file_main
        return pymzml.run.Reader(ms_file_name, MS1_Precision = 10e-6, MSn_Precision = 20e-6)
    except IOError as ioerr:
        print('File error: '+ioerr)
        return(None)

def write_csv_file(ms_file_main, data_list):
    with open(ms_file_main+'-v21.csv','a',newline='') as csv_file:
        write_data = [data_list]
        csvout=csv.writer(csv_file)
        print(write_data)
        csvout.writerows(write_data)

def main(ms_file):
    run = open_ms_file(ms_file)
    print('Now processing %s ...' %(ms_file+'.mzML'))
    for spectrum in run:
        if isinstance(spectrum['id'],str):continue
        #if spectrum['id'] < 5000: continue
        print(spectrum['id'])
        if spectrum['ms level'] == 2:
            test = []
            test.append(spectrum.hasPeak(1016.5921))
            test.append(spectrum.hasPeak(903.5080))
            test.append(spectrum.hasPeak(772.4676))
            test.append(spectrum.hasPeak(715.4461))
            test.append(spectrum.hasPeak(602.3620))
            test.append(spectrum.hasPeak(474.3035))
            test.append(spectrum.hasPeak(361.2194))
            test.append(spectrum.hasPeak(262.1510))
            test.append(spectrum.hasPeak(175.1190))
            

            tests = 0
            for test_Peak in test:
                if test_Peak:
                    tests += 1
           

            hightP=sorted(spectrum.highestPeaks(5),key = lambda x:x[1],reverse = True)
            highestPeak_mz = []
            write=[]
            if tests >= 6:
                for mz, i in hightP:
                    highestPeak_mz.append(mz)
                write.append(spectrum['id'])
                write.append(spectrum['scan start time'])
                write.append(spectrum["precursors"][0]['charge'])
                write.append(spectrum["precursors"][0]['mz'])
                write.extend(highestPeak_mz)
                print(spectrum['id'])
                write_csv_file(ms_file, write)

mzfiles = FileList()
print("We will process this file:", mzfiles)
for mzfile in mzfiles:
    try:
        main(mzfile)
    except:
        continue