#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 16:22:47 2020

@author: alexranieri
"""

# Imports
import argparse, os
import pandas as pd
from shutil import copyfile

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
def getChromLengrh():
    """
        Get chromosomes length from temp_file
        Return a dictionary: ID: [start, end]
        Use field ##sequence-region; details in 
        https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    """ 
    Dict = {}
    input = open('temp.gff', 'r')
    for line in input:
        line = line.rstrip().split(' ')
        if (line[0] == '##sequence-region'):
            Dict[line[1]] = [int(line[2]), int(line[3])]
    return Dict   
    
def write_for_sort(file):
    """
        Write the rest of the lines to a temp_file to be sorted
    """
    for_sort = open('toSort.gff', 'w')
    for_sort.writelines(file.readlines())
    
def sort_gff(file):
    """ 
        Sort the temp GFF file using pandas
        Input: GFF file
        Output: temp sorted GFF file
    """
    df = pd.read_csv(file, sep='\t', header=None) # open file as DataFrame
    sorted_df = df.sort_values([0,3,4]) # sort DataFrame using ChromID, ChromStart and ChromEnd
    sorted_df.to_csv('sorted.gff', sep ='\t', header = False, index = False) # write to file
      
def get_lineFeatures(line):
    """
        Get features of the GFF line. 
        Columns are described in: http://gmod.org/wiki/GFF3 
        Return an array
    """
    feature = line.split("\t") # array with features
    return feature

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
    polType = 'polycistron-' + currentPC[2]
    currentPC[2] = polType
    count_id = str(len(id_array))
    if (polType == 'polycistron-CDS'):
        global CDS
        CDS += 1
        description = 'ID=' + polType + '_' + str(CDS) + ';contentCount=' + count_id + ';content='
    elif(polType == 'polycistron-ncRNA'):
        global ncRNA
        ncRNA += 1
        description = 'ID=' + polType + '_' + str(ncRNA) + ';contentCount=' + count_id + ';content='
    elif(polType == 'polycistron-rRNA'):
        global rRNA
        rRNA += 1
        description = 'ID=' + polType + '_' + str(rRNA) + ';contentCount=' + count_id + ';content='
    elif(polType == 'polycistron-snoRNA'):
        global snoRNA
        snoRNA += 1
        description = 'ID=' + polType + '_' + str(snoRNA) + ';contentCount=' + count_id + ';content='
    else:
        global tRNA
        tRNA += 1
        description = 'ID=' + polType + '_' + str(tRNA) + ';contentCount=' + count_id + ';content='
    currentPC[1] = 'annotatePolycistron'
    id_array = ','.join(id_array)
    description = description + id_array
    currentPC.append(description)
    currentPC = '\t'.join(currentPC)
    file.write(currentPC + '\n')

def processSSR(poly1, poly2, file):
    """
        Process SSR and then print to file
    """
    # check if it is dSSR, cSSR or Intergenic Polycistron Region 
    if (poly1[6] == '+' and poly2[6] == '-'): # then it is cSSR
        global cSSR
        cSSR += 1
        id = 'ID=cSSR_' + str(cSSR)
        type = 'cSSR'
    elif (poly1[6] == '-' and poly2[6] == '+'): # then it is dSSR
        global dSSR
        dSSR += 1
        id = 'ID=dSSR_' + str(dSSR)
        type = 'dSSR'
    else: # the it is Intergenic Polycistron Region
        global intergenicPoly
        intergenicPoly += 1
        id = 'ID=Intergenic_Polycistron_Region_' + str(intergenicPoly)
        type = 'Intergenic_Polycistron_Region'
    # SSR features
    chrom = poly1[0]
    prog = 'annotatePolycistron'
    start = str(int(poly1[4]) + 1)
    end = str(int(poly2[3]) - 1)
    poly1Id = poly1[8].split(';')
    poly1Id = poly1Id[0][3:]
    poly2Id = poly2[8].split(';')
    poly2Id = poly2Id[0][3:]
    description = id + ';adjacentPol=' + poly1Id + ',' + poly2Id
    SSR = [chrom,prog,type,start,end,'.\t.\t.',description]
    SSR = '\t'.join(SSR)
    file.write(SSR + '\n')
    
    # insert TSS if feature is dSSR
    if (type == 'dSSR'):
        # TSS features
        global TSS
        TSS += 1
        midpoint = (int(end) + int(start))/2
        size = int(end) - int(start)
        tss_start = str(int(midpoint - size/4))  
        tss_end = str(int(midpoint + size/4))
        tss_id = 'ID=TSS_' + str(TSS) + ';adjacentPol=' + poly1Id + ',' + poly2Id
        TSS_desc = [chrom,prog,'TSS',tss_start,tss_end,'.\t.\t.',tss_id]
        TSS_desc = '\t'.join(TSS_desc)
        file.write(TSS_desc + '\n')
        
    
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
    toSort = 'toSort.gff'
    sort_gff(toSort)  
    os.remove('toSort.gff')
    
    # start read the sorted temp file
    sorted_input = open('sorted.gff', 'r')
    temp_out = open('temp_2.gff', 'a') # open temp_output file
    line = sorted_input.readline() # read initial line 
    desired_features = ['CDS','ncRNA','rRNA','snoRNA','tRNA']
    # counters
    CDS = 0
    ncRNA = 0
    rRNA = 0
    snoRNA = 0
    tRNA = 0
    print("Processing eventual polycistrons...")
    # get initial polycistron
    features = get_lineFeatures(line) # recover features of a line
    # check if line contains desired features to process
    feature_type = features[2]   
    while (feature_type not in desired_features):
        features = get_lineFeatures(line) 
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
            if (chrom != processingPC[0] or feature_type != processingPolType or strand != processingPolStrand):
                temp_out.write(line)
                processEventualPolycistron(processingPC, processingPolType, processingId_array, temp_out) # write to file 
                # store data for new eventual polycistron and read next line
                processingPC, processingId_array = storePolycistronData(features)
                processingPolType = processingPC[2]
                processingPolStrand = processingPC[6]
                line = sorted_input.readline()
            else: # we need to update the processing polycistron and read next line
                temp_out.write(line)
                processingPC[4] = features[4] # update end
                # append feature id to polycistron Note
                idToAppend = features[-1].split(';')
                idToAppend = idToAppend[0].replace('ID=','')
                processingId_array.append(idToAppend)
                line = sorted_input.readline()
        else: # write line in temp_out and continue to read next line
            temp_out.write(line)
            line = sorted_input.readline()
            
    processEventualPolycistron(processingPC, processingPolType, processingId_array, temp_out)
    print("Done!")
    temp_out.close()
    sorted_input.close()   
    
    # Start to annotate SSR, TSS and TTS
    print("Processing SSR regions...")
    # Counters
    dSSR = 0
    cSSR = 0
    intergenicPoly = 0
    TSS = 0 
    # Get chroms length
    chromLength = getChromLengrh()
    # Copy temp_2.gff to temp_3.gff to read and reopen temp_2.gff
    copyfile('temp_2.gff', 'temp_3.gff')
    temp_input = open('temp_3.gff', 'r')
    temp_out = open('temp_2.gff', 'a') 
    line = temp_input.readline() # read initial line 
    # Read file line by line processing only lines containing 'annotatePolycistron'
    # Get first Polycistron data
    for line in temp_input:
        if ('annotatePolycistron' in line):
            # get poly info
            poly1 = get_lineFeatures(line)
            poly1Chrom = poly1[0]
            line = temp_input.readline()
            break
    
    # Read rest of the file, comparing polycistron data 
    while line:
        if ('annotatePolycistron' in line): # process
            # get polycistron 2 data to compare    
            poly2 = get_lineFeatures(line)
            poly2Chrom = poly2[0]
            # check for same chrom
            if (poly1Chrom == poly2Chrom):
                processSSR(poly1, poly2, temp_out)
                poly1 = poly2
                poly1Chrom = poly1[0]
                line = temp_input.readline()
            else: # reached chrom border
            # pass features to poly1 in order to process next polycistron
                poly1 = poly2
                poly1Chrom = poly1[0]
                line = temp_input.readline()
        else:
           line = temp_input.readline() 
    
    temp_out.close()
    temp_input.close()
    print("Done!")
    
    # Create final output file and cleaning
    finalSort = 'temp_2.gff'
    sort_gff(finalSort) # sorting output
    file = open('sorted.gff', 'r')
    output = open('temp.gff', 'a')
    output.writelines(file.readlines()) # merging sorted output with headlines
    file.close()
    output.close()
    os.rename('temp.gff', args.gff_out)
    os.remove('sorted.gff')
    os.remove('temp_2.gff')
    os.remove('temp_3.gff')
    # Log
    print("Number of annotated features:")
    print("Polycistron CDS: " + str(CDS))
    print("Polycistron tRNA: " + str(tRNA))
    print("Polycistron ncRNA: " + str(ncRNA))
    print("Polycistron rRNA: " + str(rRNA))
    print("Polycistron snoRNA: " + str(snoRNA))
    print("dSSR: " + str(dSSR))
    print("cSSR: " + str(cSSR))
    print("TSS: " + str(TSS))
    print("Intergenic Polycistron Region: " + str(intergenicPoly))