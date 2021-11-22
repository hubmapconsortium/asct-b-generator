#!/usr/bin/python3

import sys
import re

DEBUG = False

# This should be an argument
#HEAD = "ovary"
HEAD = "head"

# AS_LEVELS should be arguments or set programmatically. The CT_LEVELS
# should also be 1, based on my understanding of the ASCT+B tables.
AS_LEVELS = 3
CT_LEVELS = 1

# compute these based on the max number of elements across cells
BGene_levels = 0
BProtein_levels = 0

input_dict = {}

def print_details(entities, min_level, level):
    if entities:
        out_string = ""
        while level < min_level:
            out_string += "\t\t\t"
            level += 1

        for entity in entities:
            out_string += "\t" + entity + "\t" + input_dict[entity]['label'] + "\t" + input_dict[entity]['id']
            level += 1

        out_file.write(out_string)
    return level

def output_struct(parents):
    # we want to go from the top level structure to the bottom level
    parents.reverse()

    level = 0
    for key in parents:
        if DEBUG:
            print(key, end="\t")

        s_type = input_dict[key]['s_type']
        if s_type == "AS":
            # process anatomical structures
            if level == 0:
                out_string = key + "\t" + input_dict[key]['label'] + "\t" + input_dict[key]['id']
            else:
                out_string = "\t" + key + "\t" + input_dict[key]['label'] + "\t" + input_dict[key]['id']
            level += 1

        elif s_type == "CT":
            # process cell types
            out_string = ""
            while level < AS_LEVELS:
                out_string += "\t\t\t"
                level += 1
            out_string += "\t" + key + "\t" + input_dict[key]['label'] + "\t" + input_dict[key]['id']
            level += 1

        else:
            # children should only be of type AS or CT
            print("ERROR: erroneous child value", input_dict[key])
            continue

        out_file.write(out_string)
        level = print_details(input_dict[key]['genes'], (AS_LEVELS + CT_LEVELS), level)
        level = print_details(input_dict[key]['proteins'], (AS_LEVELS + CT_LEVELS + BGene_levels), level)

    if DEBUG:
        print()
    out_file.write("\n")

    # restore list order
    parents.reverse()


# this is very inefficient but optimized for easy inputs and assumed processing speed is irrelevant.
def get_parents(leaf, parents):
    # loop through all stuctures and print any that contain this leaf as a child
    for key in input_dict.keys():
        children = input_dict[key]["children"]
        if not children:
            continue
        if leaf in children:
            parents.append(key)
            get_parents(key, parents)

    # if at the top level, then output the list
    if leaf == HEAD:
        output_struct(parents)

    # remove latest structure as we back up the recursion.
    parents.pop()


# execute script
if __name__ == "__main__":
    # open input file. argv[0] is the program name
    in_filename = sys.argv[1]
    in_file = open(in_filename, "r")

    # this will overwrite any existing file
    out_filename = sys.argv[2]
    out_file = open(out_filename, "w")

    contents = in_file.readlines()
    for line in contents:
        # the first line might be a header line. That's intrisically resolved as the header "shouldn't" match to any children
        name, label, reference, s_type, children_string, genes_string, proteins_string = re.split(r'\t', line.rstrip('\n'))

        # skip header line
        if name.lower() == "name":
            continue

        # remove the quotes and extra spaces from the TSV file
        children_string = children_string.replace('"', '').replace(', ', ',')
        genes_string = genes_string.replace('"', '').replace(', ', ',')
        proteins_string = proteins_string.replace('"', '').replace(', ', ',')

        # convert from a string into a list
        children = []
        if children_string:
            children = children_string.split(',')
            
        genes = []
        if genes_string:
            genes = genes_string.split(',')
            if BGene_levels < len(genes):
                BGene_levels = len(genes)
                
        proteins = []
        if proteins_string:
            proteins = proteins_string.split(',')
            if BProtein_levels < len(proteins):
                BProtein_levels = len(proteins)

        # add anatomical structure to our dictionary
        input_dict.update({name:{"label":label, "id":reference, "s_type":s_type, "children":children, "genes":genes, "proteins":proteins}})

    # all structures loaded so now we step though each structure and print out all leaves
    for key in input_dict.keys():
        children = input_dict[key]["children"]
        s_type = input_dict[key]['s_type']
        # build hierarchical structure with anatomical structures and cells.
        if not children and (s_type in ("AS", "CT")):
            # no children so end point and need to print upstream structures
            parents = [key]
            get_parents(key, parents)

    in_file.close()
    out_file.close()

