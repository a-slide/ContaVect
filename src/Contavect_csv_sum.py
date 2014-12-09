#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import argv as argv
import os, csv, glob

""" Draft function to summarize the number of reads mapped in distribution files from several runs
of Contavect in a single comprehensive csv file"""

def main ():
    """Find csv files with a common patern fetch the 2 first row of the first alphabetical
    file and the 3rd row of all files. Merge everything in a common csv file"""

    ref_all = []

    # Find ref matching the patern and sorting the list alphabetically
    ref = list(glob.iglob("*"+argv[1]))
    ref.sort()

    print (ref)

    # Fetch first column (ref names)
    with open(ref[0], newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        ref_all.append([row[0] for row in reader])

    # Fetch second column (ref len)
    with open(ref[0], newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        ref_all.append([row[1] for row in reader])

    # Fetch Third column of all reference(ref len)
    for r in ref:
        with open(r, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            # Parse the nread row
            myrow = [row[2] for row in reader]
            # Replace nread by the name of the ref minus the common suffix
            myrow[0] = r.replace(argv[1], '')
            ref_all.append(myrow)

    # Transpose the table
    t_ref_all = [[x[i] for x in ref_all] for i in range(len(ref_all[0]))]

    # Finally write a new table
    with open(argv[2], 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in t_ref_all:
            writer.writerow(i)

    print("done")

    exit (0)


def usage():
    """Simple usage function"""

    print ("Usage: ", argv[0], "<Patern of csv files to match> <output name of the csv file>")
    print ("\tExample : ", argv[0], " _Reference_distribution.csv  ALL_Reference_distribution.csv")

if __name__ == '__main__':
    if len(argv) < 2:        # if not enought arg call usage function
        usage()
    else:
        main()              # else call the main function
