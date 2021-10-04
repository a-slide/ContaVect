#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import argv as argv
import os, csv, glob

""" Draft function to summarize the coverage over a ref in multiple bed files created by Contavect
in a single comprehensive csv file"""

def main ():
    """Find bed files with a common pattern fetch the 2nd row of the first alphabetical
    file and the 4rth row of all files. Merge everything in a common csv file
    Create a additional """

    # Find ref matching the patern and sorting the list alphabetically
    ref = list(glob.iglob("*"+argv[1]))
    ref.sort()

    print (ref)

    # Fetch second column (ref position)
    with open(ref[0], newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')

        # Get rid of the bed header
        title = next(reader)

        # Fetch the position row and store it in 2 list for raw count and norm count
        pos_row = [row[1] for row in reader]
        all_count = [["", "Position"] + pos_row]
        all_norm =  [["", "Position"] + pos_row]

    # Fetch fourth column of all reference(coverage)
    for r in ref:
        with open(r, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')

            # Get rid of the bed header
            title = next(reader)

            # Extract ref name
            ref_name = r.replace(argv[1], '')

            # Parse the count row
            count_row = [int(row[3]) for row in reader]
            all_count.append([ref_name, "Count"] + count_row)

            # Create a normalized count row
            row_sum = sum (count_row)
            norm_row = [i/row_sum for i in count_row]
            all_norm.append([ref_name, "Norm Count"] + norm_row)


    # Fuse count and norm count tables and transpose the table
    all_data = all_count+all_norm
    all_data = [[x[i] for x in all_data] for i in range(len(all_data[0]))]

    # Finally write a new table
    with open(argv[2], 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for i in all_data:
            writer.writerow(i)

    print("done")

    exit (0)


def usage():
    """Simple usage function"""

    print(("Usage: ", argv[0], "<Pattern of bed files to match> <output name of the csv file>"))
    print(("\tExample : ", argv[0], " .AAV.csv  ALL_AAV.csv"))
    exit(1)

if __name__ == '__main__':
    if len(argv) <= 2:        # if not enough arg call usage function
        usage()
    else:
        main()              # else call the main function
