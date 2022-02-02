#!/usr/bin/python3

import sys
import re

DEBUG = False

# issues
# 1. The user needs to know how many levels for the anatomical structure or at least over estimate the number
# 2. The program doesn't insert a header line

# usage
# ./process.py ovary 10 Ovaries-v2.txt Ovaries-v2-ASCTB.xls

# Ovary header
# AS/1	AS/1/LABEL	AS/1/ID	AS/2	AS/2/LABEL	AS/2/ID	AS/3	AS/3/LABEL	AS/3/ID	AS/4	AS/4/LABEL	AS/4/ID	AS/5	AS/5/LABEL	AS/5/ID	AS/6	AS/6/LABEL	AS/6/ID	AS/7	AS/7/LABEL	AS/7/ID	AS/8	AS/8/LABEL	AS/8/ID	AS/9	AS/9/LABEL	AS/9/ID	AS/10	AS/10/LABEL	AS/10/ID	CT/1	CT/1/LABEL	CT/1/ID	BGene/1	BGene/1/LABEL	BGene/1/ID	BProtein/1	BProtein/1/LABEL	BProtein/1/ID	BProtein/2	BProtein/2/LABEL	BProtein/2/ID	BProtein/3	BProtein/3/LABEL	BProtein/3/ID	BProtein/4	BProtein/4/LABEL	BProtein/4/ID	BProtein/5	BProtein/5/LABEL	BProtein/5/ID	BProtein/6	BProtein/6/LABEL	BProtein/6/ID	BProtein/7	BProtein/7/LABEL	BProtein/7/ID	BProteoform/1	BProteoform/1/LABEL	BProteoform/1/ID	BLipid/1	BLipid/1/LABEL	BLipid/1/ID	BMetabolites/1	BMetabolites/1/LABEL	BMetabolites/1/ID	FTU/1	FTU/1/LABEL	FTU/1/ID	REF/1	REF/1/DOI	REF/1/NOTES	REF/2	REF/2/DOI	REF/2/NOTES	REF/3	REF/3/DOI	REF/3/NOTES	REF/4	REF/4/DOI	REF/4/NOTES	REF/5	REF/5/DOI	REF/5/NOTES	REF/6	REF/6/DOI	REF/6/NOTES

# Fallopian tube header
# AS/1	AS/1/LABEL	AS/1/ID	AS/2	AS/2/LABEL	AS/2/ID	AS/3	AS/3/LABEL	AS/3/ID	AS/4	AS/4/LABEL	AS/4/ID	AS/5	AS/5/LABEL	AS/5/ID	CT/1	CT/1/LABEL	CT/1/ID	BGene/1	BGene/1/LABEL	BGene/1/ID	BGene/2	BGene/2/LABEL	BGene/2/ID	BProtein/1	BProtein/1/LABEL	BProtein/1/ID	BProtein/2	BProtein/2/LABEL	BProtein/2/ID	BProtein/3	BProtein/3/LABEL	BProtein/3/ID	BProteoform/1	BProteoform/1/LABEL	BProteoform/1/ID	BLipid/1	BLipid/1/LABEL	BLipid/1/ID	BMetabolites/1	BMetabolites/1/LABEL	BMetabolites/1/ID	FTU/1	FTU/1/LABEL	FTU/1/ID	REF/1	REF/1/DOI	REF/1/NOTES	REF/2	REF/2/DOI	REF/2/NOTES	REF/3	REF/3/DOI	REF/3/NOTES	REF/4	REF/4/DOI	REF/4/NOTES	REF/5	REF/5/DOI	REF/5/NOTES	REF/6	REF/6/DOI	REF/6/NOTES

# Uterus header
# AS/1	AS/1/LABEL	AS/1/ID	AS/2	AS/2/LABEL	AS/2/ID	AS/3	AS/3/LABEL	AS/3/ID	AS/4	AS/4/LABEL	AS/4/ID	AS/5	AS/5/LABEL	AS/5/ID	AS/6	AS/6/LABEL	AS/6/ID	AS/7	AS/7/LABEL	AS/7/ID	CT/1	CT/1/LABEL	CT/1/ID	BGene/1	BGene/1/LABEL	BGene/1/ID	BGene/2	BGene/2/LABEL	BGene/2/ID	BGene/3	BGene/3/LABEL	BGene/3/ID	BGene/4	BGene/4/LABEL	BGene/4/ID	BGene/5	BGene/5/LABEL	BGene/5/ID	BGene/6	BGene/6/LABEL	BGene/6/ID	BProtein/1	BProtein/1/LABEL	BProtein/1/ID	BProtein/2	BProtein/2/LABEL	BProtein/2/ID	BProtein/3	BProtein/3/LABEL	BProtein/3/ID	BProtein/4	BProtein/4/LABEL	BProtein/4/ID	BProtein/5	BProtein/5/LABEL	BProtein/5/ID	BProtein/6	BProtein/6/LABEL	BProtein/6/ID	REF/1	REF/1/DOI	REF/1/NOTES	REF/2	REF/2/DOI	REF/2/NOTES

# this will contain the top level anatomical structure (e.g., "ovary")
head = ""

# AS_levels needs to be provided by the user as I don't think we can
# compute this without generating the full tree twice, which we might
# eventually need to do.
AS_levels = 0

# The CT_LEVELS hould also be 1, based on my understanding of the ASCT+B tables.
CT_LEVELS = 1

# compute these based on the max number of elements across cells
BGene_levels = 0
BProtein_levels = 0
BProteoform_levels = 0
BLipid_levels = 0
BMetabolites_levels = 0
# we assume FTU has the same three columns as the other descriptors.
FTU_levels = 0
# we don't actually need the REF levels since REF is the last set,
# however, it's a hack to allow us to handle REF like the other
# descriptors.
REF_levels = 0

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
            while level < AS_levels:
                out_string += "\t\t\t"
                level += 1
            out_string += "\t" + key + "\t" + input_dict[key]['label'] + "\t" + input_dict[key]['id']
            level += 1

        else:
            # children should only be of type AS or CT
            print("ERROR: erroneous child value", input_dict[key])
            continue

        out_file.write(out_string)
        level = print_details(input_dict[key]['genes'],
                              (AS_levels + CT_LEVELS), level)
        level = print_details(input_dict[key]['proteins'],
                              (AS_levels + CT_LEVELS + BGene_levels), level)
        level = print_details(input_dict[key]['proteoforms'],
                              (AS_levels + CT_LEVELS + BGene_levels + BProtein_levels), level)
        level = print_details(input_dict[key]['lipids'],
                              (AS_levels + CT_LEVELS + BGene_levels + BProtein_levels + BProteoform_levels), level)
        level = print_details(input_dict[key]['metabolites'],
                              (AS_levels + CT_LEVELS + BGene_levels + BProtein_levels + BProteoform_levels + BLipid_levels), level)
        level = print_details(input_dict[key]['ftu'],
                              (AS_levels + CT_LEVELS + BGene_levels + BProtein_levels + BProteoform_levels + BLipid_levels + BMetabolites_levels), level)
        level = print_details(input_dict[key]['refs'],
                              (AS_levels + CT_LEVELS + BGene_levels + BProtein_levels + BProteoform_levels + BLipid_levels + BMetabolites_levels + FTU_levels), level)

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
    if leaf.lower() == head.lower():
        output_struct(parents)

    # remove latest structure as we back up the recursion.
    parents.pop()


def split_string(array_string, max_level):
    tmp = []
    if array_string:
        tmp = array_string.split(',')
        if max_level < len(tmp):
            max_level = len(tmp)

    return tmp, max_level


# execute script
if __name__ == "__main__":
    # open input file. argv[0] is the program name

    # this should be the top level anatomical structure (e.g., "ovary")
    head = sys.argv[1]

    # number AS levels
    AS_levels = int(sys.argv[2])
    
    in_filename = sys.argv[3]
    in_file = open(in_filename, "r")

    # this will overwrite any existing file
    out_filename = sys.argv[4]
    out_file = open(out_filename, "w")

    # stop if input file is incorrect
    found_input_error = False
    
    contents = in_file.readlines()
    headerLine = True
    for line in contents:
        # the first line might be a header line. That's intrisically resolved as the header "shouldn't" match to any children
        name, label, reference, s_type, children_string, genes_string, proteins_string, proteoforms_string, lipids_string, metabolites_string, ftu_string, refs_string = re.split(r'\t', line.rstrip('\n'))

        # first line is the header, so skip it
        if headerLine:
            headerLine = False
            continue

        # remove the quotes and extra spaces from the TSV file
        name = name.rstrip()
        children_string = children_string.replace('"', '').replace(', ', ',')
        genes_string = genes_string.replace('"', '').replace(', ', ',')
        proteins_string = proteins_string.replace('"', '').replace(', ', ',')
        proteoforms_string = proteoforms_string.replace('"', '').replace(', ', ',')
        lipids_string = lipids_string.replace('"', '').replace(', ', ',')
        metabolites_string = metabolites_string.replace('"', '').replace(', ', ',')
        ftu_string = ftu_string.replace('"', '').replace(', ', ',')
        refs_string = refs_string.replace('"', '').replace(', ', ',')

        # convert from a string into a list
        children = []
        if children_string:
            children = children_string.split(',')

        # process descriptor arrays
        genes, BGene_levels = split_string(genes_string, BGene_levels)
        proteins, BProtein_levels = split_string(proteins_string, BProtein_levels)
        proteoforms, BProteoform_levels = split_string(proteoforms_string, BProteoform_levels)
        lipids, BLipid_levels = split_string(lipids_string, BLipid_levels)
        metabolites, BMetabolites_levels = split_string(metabolites_string, BMetabolites_levels)
        ftu, FTU_levels = split_string(ftu_string, FTU_levels)
        refs, REF_levels = split_string(refs_string, REF_levels)

        if children and (genes or proteins or proteoforms or lipids or metabolites or ftu or refs):
            print("ERROR: genes, proteins, proteoforms, lipids, metabolites, ftu, and refs can only be applied to structures without children or to cell types")
            print("Structure: ", name)
            print("Children: ", children)
            print("Genes: ", genes)
            print("Proteins: ", proteins)
            print("Proteoforms: ", proteoforms)
            print("Lipids: ", lipids)
            print("Metabolites: ", metabolites)
            print("FTU: ", ftu)
            print("References: ", refs)
            found_input_error = True
                    
        # add anatomical structure to our dictionary
        input_dict.update({name:{"label":label, "id":reference, "s_type":s_type, "children":children, "genes":genes, "proteins":proteins, "proteoforms":proteoforms, "lipids":lipids, "metabolites":metabolites, "ftu":ftu, "refs":refs }})

    if not found_input_error:
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

    if DEBUG:
        print(AS_levels, CT_LEVELS, BGene_levels, BProtein_levels, BProteoform_levels, BLipid_levels, BMetabolites_levels, FTU_levels)
