#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 16:22:47 2020

@author: alexranieri
"""

# Imports
import argparse, os

# Functions

def read_input(inputfile):
    """ 
        Read the input file passed by args.gff_file
        Return: file object
    """
    file = open(inputfile, 'r')
    return file

def write_headlines_tmp(file):
    """
        Write the lines starting with # to a temp_file
    """
    temp_file = open('temp.gff', 'a')
    for line in file:
        if (line[0] == '#'):
            temp_file.write(line)
        else:
            break
        
def write_for_sort(file):
    """
        Write the rest of the lines to a temp_file to be sorted
    """
    for_sort = open('toSort.gff', 'w')
    for_sort.writelines(file.readlines())
    
def sort_gff():
    """ 
        Sort the temp GFF file from write_for_sort function using bash sort
        Input: toSort.gff
        Output: temp sorted GFF file
    """
    os.system('sort -k1,1 -k4n,4 -k5n,5 toSort.gff > sorted.gff')
    os.system('rm toSort.gff')

def get_lineFeatures(line):
    """
        Get features of the GFF line. 
        Columns are described in: http://gmod.org/wiki/GFF3 
        Return an array
    """
    feature = line.split("\t") # array with features
    return feature

    feature_type = features[2]
    if (feature_type in desired_features):
            # store data for eventual polycistron 
            currentPC, id_array = storePolycistronData(features)
    else: # write line in temp_out and continue to read next line
            temp_out.write(line)
            line = sorted_input.readline()
def storePolycistronData(features):
    """
        Get features array of line. 
        Columns are described in: http://gmod.org/wiki/GFF3 
        Return:
            1 - array containing polycistron data
            2 - an array containing IDs in that polycistron
    """
    polyFeatures = features[0:-1] 
    idContents = features[-1].split(';') # get feature ID
    idContents = [idContents[0].replace('ID=','')]
    return polyFeatures, idContents
    
def processEventualPolycistron(currentPC, polType, id_array, file):
    """
        Print current polycistron to file
    """
    polType = 'Polycistron-' + currentPC[2]
    currentPC[2] = polType
    if (polType == 'Polycistron-CDS'):
        global CDS
        CDS += 1
        description = 'ID=' + polType + '_' + str(CDS) + ';Note='
    elif(polType == 'Polycistron-ncRNA'):
        global ncRNA
        ncRNA += 1
        description = 'ID=' + polType + '_' + str(ncRNA) + ';Note='
    elif(polType == 'Polycistron-rRNA'):
        global rRNA
        rRNA += 1
        description = 'ID=' + polType + '_' + str(rRNA) + ';Note='
    elif(polType == 'Polycistron-snoRNA'):
        global snoRNA
        snoRNA += 1
        description = 'ID=' + polType + '_' + str(snoRNA) + ';Note='
    else:
        global tRNA
        tRNA += 1
        description = 'ID=' + polType + '_' + str(tRNA) + ';Note='
    currentPC[1] = 'annotatePolycistron'
    id_array = ','.join(id_array)
    description = description + id_array
    currentPC.append(description)
    currentPC = '\t'.join(currentPC)
    file.write(currentPC + '\n')
        
# main program starts here
if __name__ == '__main__':
    
    # Defining arguments, usage and help message 
    message = 'This script will annotate the polycistrons along with dSSR and cSSR\
        (divergent and convergent Strand Switch Region) in trypanosomatid genomes.\
            Type of annotated Polycistrons: CDS, ncRNA, rRNA, snoRNA and tRNA'
    parser = argparse.ArgumentParser(prog = 'annotatePolycistron_SSR.py',
                                     description= message)
    # Set arguments
    parser.add_argument('-g', '--gff', action = 'store', dest = 'gff_file', 
                        required = True, help = 'GFF input file name', type=str)
    help_message = 'GFF output file name. \
        Default is: {input}_polycistronAnnotated.gff'
    parser.add_argument('-o', '--output', action = 'store', dest = 'gff_out', 
                        required = False, help = help_message, type=str,
                        default = None)
    # Create variables
    args = parser.parse_args() 
    if args.gff_out is None:
        args.gff_out = args.gff_file[:-4] + '_polycistronAnnotated.gff'
    
    # gff file content
    print("Reading input file...")
    gff = read_input(args.gff_file)
    
    # separate the lines stating with '#' to a temp file
    write_headlines_tmp(gff) 
    
    # write the rest of the file (which contain the features) for sort
    write_for_sort(gff)
    gff.close() # close input file - no longer needed
    
    # sort the file created with the write_for_sort function
    print("Sorting input file...")
    sort_gff()  
    
    # start read the sorted temp file
    sorted_input = open('sorted.gff', 'r')
    temp_out = open('temp.gff', 'a') # open temp_output file
    line = sorted_input.readline() # read initial line 
    desired_features = ['CDS','ncRNA','rRNA','snoRNA','tRNA']
    # counters
    CDS = 0
    ncRNA = 0
    rRNA = 0
    snoRNA = 0
    tRNA = 0
    print("Processing...")
    # get initial polycistron
    features = get_lineFeatures(line) # recover features of a line
    # check if line contains desired features to process
    feature_type = features[2]        
    if (feature_type in desired_features):
        # store data for eventual polycistron 
        processingPC, processingId_array = storePolycistronData(features)
        processingPolType = processingPC[2]
        processingPolStrand = processingPC[6]
        temp_out.write(line)
        line = sorted_input.readline()                     
    else: # write line in temp_out and continue to read next line
        temp_out.write(line)
        line = sorted_input.readline()
        
    while line: # loop trough all lines in sorted input file
        features = get_lineFeatures(line) # recover features of a line
        # check if line contains desired features to process
        feature_type = features[2]        
        if (feature_type in desired_features):
            # check if lines contains same chromosome ID, feature type or strand
            # If one of those is different, thre is an eventual polycistron boundary.
            chrom = features[0]
            strand = features[6]
            temp_out.write(line)
            if (chrom != processingPC[0] or feature_type != processingPolType or strand != processingPolStrand):
                processEventualPolycistron(processingPC, processingPolType, processingId_array, temp_out) # write to file 
                # store data for new eventual polycistron and read next line
                processingPC, processingId_array = storePolycistronData(features)
                processingPolType = processingPC[2]
                processingPolStrand = processingPC[6]
                line = sorted_input.readline()
            else: # we need to update the processing polycistron and read next line
                processingPC[4] = features[4] # update end
                # append feature id to polycistron Note
                idToAppend = features[-1].split(';')
                idToAppend = idToAppend[0].replace('ID=','')
                processingId_array.append(idToAppend)
                temp_out.write(line)
                line = sorted_input.readline()
        else: # write line in temp_out and continue to read next line
            temp_out.write(line)
            line = sorted_input.readline()
            
    processEventualPolycistron(processingPC, processingPolType, processingId_array, temp_out)
    print("Done!")
    # Create final output file and cleaning
    temp_out.close()
    sorted_input.close()
    os.system('mv temp.gff ' + args.gff_out)
    os.system('rm sorted.gff')
    print("Number of annotated features:")
    print("CDS: " + str(CDS))
    print("tRNA: " + str(tRNA))
    print("ncRNA: " + str(ncRNA))
    print("rRNA: " + str(rRNA))
    print("snoRNA: " + str(snoRNA))