#-*- coding: utf-8 -*-
# AlFinder strategy 1
# Version 0.11
# This script reads the *.mzML (all mzML file in current folder) MS data, and extracts the MS/MS information of
# SLGGDSIMGIQLVSR peptide and its modified analogue according to the tag ions
#
# Tag ions include y9, y8, y7, y5, and y4
#
# The result writes in a *_S1.csv file, * is the name of .mzML file
# The recoding information include spectrumID, precursor's retention time, precursor's charge, precursor's m/z,
# the 5 highest peaks of MS/MS spectrum , the theoretical ppant ejection m/z (based on the m/z of precursors),
# whether the theoretical ppant ejection is found in ms/ms spectrum and the relative intension of ppant ejection
# ion according to base peak
#
# Also try to write a new mzML file include the selective ms/ms spectrums, the name is *_S1.mzML (ongoing)
#
# Based on hasPeak.py from example script of pymzML
# Author: Bo Pang (SIOC)
# Bless!



#import sys
import pymzml
import csv
#import os
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
    with open(ms_file_main+'_S1.csv','a',newline='') as csv_file:
        write_data = [data_list]
        csvout=csv.writer(csv_file)
        print(write_data)
        csvout.writerows(write_data)

def theo_ppant_ejection(precursor_mz, precursor_charge):
    return precursor_mz * precursor_charge - 1.00783 * precursor_charge - 1610.7619

#def writemzML():
    #pass

def main(ms_file):
    run = open_ms_file(ms_file)
    #run_write = pymzml.run.Writer(filename="1.mzML", run=run, overwrite=True)
    print('Now processing %s' %(ms_file))
    for spectrum in run:
        if isinstance(spectrum['id'],str):continue
        #if spectrum['id'] < 5000: continue
        print(spectrum['id'])
        if spectrum['ms level'] == 2:
            test1 = spectrum.hasPeak(1016.5921)
            test2 = spectrum.hasPeak(903.5080)
            test3 = spectrum.hasPeak(772.4676)
            test4 = spectrum.hasPeak(602.3620)
            test5 = spectrum.hasPeak(474.3035)
            

            tests = [False, False, False, False, False]
            if test1 != []:
                tests[0] = True
            if test2 != []:
                tests[1] = True
            if test3 != []:
                tests[2] = True
            if test4 != []:
                tests[3] = True
            if test5 != []:
                tests[4] = True


            all = True
            for i in range(len(tests)):
                if not tests[i]:
                    all = False

            hightP=sorted(spectrum.highestPeaks(5),key = lambda x:x[1],reverse = True)
            highestPeak_mz = []
            writecsv=[]
            if all:
                for mz, i in hightP:
                    highestPeak_mz.append(mz)
                writecsv.append(spectrum['id'])
                writecsv.append(spectrum['scan start time'])
                writecsv.append(spectrum["precursors"][0]['charge'])
                writecsv.append(spectrum["precursors"][0]['mz'])
                writecsv.extend(highestPeak_mz)
                theo_ppant = round(theo_ppant_ejection(spectrum["precursors"][0]['mz'],spectrum["precursors"][0]['charge']),5)
                writecsv.append(theo_ppant)
                if spectrum.hasPeak(theo_ppant) !=[]:
                    writecsv.append('Eject')
                    #print(theo_ppant)
                    #print(spectrum.hasPeak(theo_ppant)[0][1])
                    #print(hightP[0][1])
                    writecsv.append(spectrum.hasPeak(theo_ppant)[0][1]/hightP[0][1]*100)
                else:
                    writecsv.append('No Eject')
                    writecsv.append('0')
                write_csv_file(ms_file, writecsv)
                #run_write.addSpec(spectrum)
                #run_write.save()

mzfiles = FileList()
print("We will process this file:", mzfiles,", based on precursor selection strategy 1")
for mzfile in mzfiles:
    try:
        main(mzfile)
    except:
        continue


