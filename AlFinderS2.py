#-*- coding: utf-8 -*-
# AlFinder strategy 2
# Version 0.12
# This script reads the *.mzML (all mzML file in current folder) MS data, and extracts the MS/MS information of
# SLGGDSIMGIQLVSR peptide and its modified analogue according to the tag ions
#
# Tag ions include y1-y9, at least 6 hits make the precursor selected
#
# The result writes in a *_S2.csv file, * is the name of .mzML file
# The recoding information include spectrumID, precursor's retention time, precursor's charge, precursor's m/z,
# the 5 highest peaks of MS/MS spectrum , the theoretical ppant ejection m/z (based on the m/z of precursors),
# whether the theoretical ppant ejection is found in ms/ms spectrum and the relative intension of ppant ejection
# ion according to base peak
#
# Also try to write a new mzML file include the selective ms/ms spectrums, the name is *_S2.mzML (ongoing)
#
# Based on hasPeak.py from example script of pymzML
# Author: Bo Pang (SIOC)
# Bless!


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
    with open(ms_file_main+'_S2.csv','a',newline='') as csv_file:
        write_data = [data_list]
        csvout=csv.writer(csv_file)
        print(write_data)
        csvout.writerows(write_data)

def theo_ppant_ejection(precursor_mz, precursor_charge):
    return precursor_mz * precursor_charge - 1.00783 * precursor_charge - 1610.7619

def writemzML():
    pass

def main(ms_file):
    run = open_ms_file(ms_file)
    print('Now processing %s ...' %(ms_file))
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
                if test_Peak !=[]:
                    tests += 1
           

            hightP=sorted(spectrum.highestPeaks(5),key = lambda x:x[1],reverse = True)
            highestPeak_mz = []
            writecsv=[]
            if tests >= 6:
                for mz, i in hightP:
                    highestPeak_mz.append(mz)
                writecsv.append(spectrum['id'])
                writecsv.append(spectrum['scan start time'])
                writecsv.append(spectrum["precursors"][0]['charge'])
                writecsv.append(spectrum["precursors"][0]['mz'])
                writecsv.extend(highestPeak_mz)
                theo_ppant = round(theo_ppant_ejection(spectrum["precursors"][0]['mz'], spectrum["precursors"][0]['charge']), 5)
                writecsv.append(theo_ppant)
                if spectrum.hasPeak(theo_ppant) != []:
                    writecsv.append('Eject')
                    print(theo_ppant)
                    print(spectrum.hasPeak(theo_ppant)[0][1])
                    print(hightP[0][1])
                    writecsv.append(spectrum.hasPeak(theo_ppant)[0][1] / hightP[0][1] * 100)
                else:
                    writecsv.append('No Eject')
                    writecsv.append('0')
                write_csv_file(ms_file, writecsv)

mzfiles = FileList()
print("We will process this file:", mzfiles, ", based on precursor selection strategy 2")
for mzfile in mzfiles:
    try:
        main(mzfile)
    except:
        continue