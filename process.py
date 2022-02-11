#!/usr/bin/python3

import sys
import re
from anytree import Node, RenderTree, AsciiStyle, LevelOrderGroupIter, Walker

'''
https://pypi.org/project/anytree/

Diffs from version 1:
  1. don't need to know the name of organ
  2. don't need to know the number of AS levels
  3. added "cells" column to the input file
  4. allow biomarkers and references to be allied to any anatomical structure

Assumptions:
  1. Inner nodes can contain biomarkers or references, not just leaves and cells.
  2. All anatomical structures must be uniquely named, for example,
     there can not be two structures called "ovary" but there can be
     "left ovary" and "right ovary".
  3. Cell type is only one level.

Error testing:
  1. make sure only one root structure
  2. make sure every structure has only one parent
  3. make sure no two structures or features have the same name
'''

################################################
# Globals

DEBUG = False

# all anatomical structures
nodes = {}

# the top level node
tree_root = ""

# all cells, biomarkers and references
features = {}

# compute these based on the max number of biomarkers across nodes and
# cells
max_genes = 0
max_proteins = 0
max_proteoforms = 0
max_lipids = 0
max_metabolites = 0
max_ftu = 0
# we don't actually need to track the max level for references as it's
# last data block, however, this is a hack to allow us to handle
# references like the biomarkers.
max_references = 0

# the file handles
in_file = ""
out_file = ""


################################################
# Classes

# used to store cells as well as biomarkers and references
class Feature:
    def __init__(self, name, feature_type, label="", id="", **kwargs):
        self.name = name
        self.feature_type = feature_type
        self.label = label
        self.id = id
        self.genes = kwargs.get('genes', None)
        self.proteins = kwargs.get('proteins', None)
        self.proteoforms = kwargs.get('proteoforms', None)
        self.lipids = kwargs.get('lipids', None)
        self.metabolites = kwargs.get('metabolites', None)
        self.ftu = kwargs.get('ftu', None)
        self.references = kwargs.get('references', None)

    def __str__(self):
        out = "[" + self.name
        out += " type:" + self.feature_type
        if self.genes:
            out += " genes:" + str(self.genes)
        out += "]"
        return out


################################################
# Functions

# test if a node contains any cells or features.
def no_features(node):
    # test for all biomarkers and references
    if node.cells or node.genes or node.proteins or node.proteoforms or node.lipids or node.metabolites or node.ftu or node.references:
        return False
    return True


# generates a triplet of column headers as exemplified below:
# AS/1	AS/1/LABEL	AS/1/ID
def get_header_block(content, depth):
    return content + "/" + str(depth) + "\t" + content + "/" + str(depth) + "/LABEL\t" + content + "/" + str(depth) + "/ID"


def get_biomarker_header(biomarker, max_depth):
    output = ""
    for depth in range(1, max_depth+1):
        output += get_header_block(biomarker, depth) + "\t"
    return output

def print_ASCTB_header(max_AS_depth):
    global max_genes, out_file

    header = ""

    # add anatomical structures
    for depth in range(1, max_AS_depth+1):
        header += get_header_block("AS", depth) + "\t"

    # assume cell type is only one depth
    header += get_header_block("CT", 1) + "\t"

    # add biomarkers
    header += get_biomarker_header("BGene", max_genes)
    header += get_biomarker_header("BProteins", max_proteins)
    header += get_biomarker_header("BProteoforms", max_proteoforms)
    header += get_biomarker_header("BLipids", max_lipids)
    header += get_biomarker_header("BMetabolites", max_metabolites)
    header += get_biomarker_header("FTU", max_ftu)

    # we handle references separately because their header details
    # differ from the biomarkers.
    for depth in range(1, max_references+1):
        header += "REF/" + str(depth) + "\t" + "REF/" + str(depth) + "/DOI\t" + "REF/" + str(depth) + "/NOTES"
        # prevent a tab at the end of the header line
        if depth < max_references:
            header += "\t"
    
    header += "\n"
    out_file.write(header)


# Generates a triplet for a data element
def get_data_block(element):
    return element.name + "\t" + element.label + "\t" + element.id


# Outputs the biomarker data and/or tabs as appropriate
def add_biomarkers(elements, max_depth):
    output = ""
    count = 0
    if elements:
        for element in elements:
            output += get_data_block(features[element]) + "\t"
            count += 1
    if count < max_depth:
        output += "\t\t\t" * (max_depth-count)
    return output


# The argument (element) can be either an anatomical structure ("node")
# or a cell as they contain the same features.
def add_features(record, element):
    global features, out_file

    # need to track how many entities we output per feature type, to
    # appropriately pad the output.
    count = 0

    # add biomarkers
    record += add_biomarkers(element.genes, max_genes)
    record += add_biomarkers(element.proteins, max_proteins)
    record += add_biomarkers(element.proteoforms, max_proteoforms)
    record += add_biomarkers(element.lipids, max_lipids)
    record += add_biomarkers(element.metabolites, max_metabolites)
    record += add_biomarkers(element.ftu, max_ftu)
    record += add_biomarkers(element.references, max_references)

    record += "\n"
    out_file.write(record)


# Print the ASCT+B table, stepping through the tree, one level at a
# time and applying each of the relevant cells and features to the
# node.
def print_ASCTB_table():
    global nodes, tree_root, features

    # get a list of lists where the sub-lists are the nodes per level
    # of the tree
    levels = [[node.name for node in children] for children in LevelOrderGroupIter(tree_root)]

    # total number of anatomical structure levels
    max_AS_depth = len(levels)

    print_ASCTB_header(max_AS_depth)

    # we need to track which AS level we're at, so we can pad as
    # needed when we add cells.
    AS_depth = 0
    w = Walker()
    for level in levels:
        AS_depth += 1
        for node_name in level:
            node = nodes[node_name]
            
            # check for features to output and if no features, then
            # skip to the next node.
            if no_features(node) and not node.is_leaf:
                continue

            # walk the tree up from the node to the root. The first
            # tuple is for upward walking and the second is just the
            # tree root. So we only care about the third tuple which
            # is the downward walking.
            walked = w.walk(tree_root, node)[2]
            # add root, since it was in the middle tuple that we
            # ignore.
            walked = (tree_root,) + walked

            # add all ancestral anatomical structures to the data
            # record being output.
            anatomical_structure = ""
            for ancestor in walked:
                anatomical_structure += get_data_block(ancestor) + "\t"

            # pad record, to account for structures that aren't as
            # deep as the maximum possible depth (max_AS_depth).
            anatomical_structure += "\t\t\t" * (max_AS_depth-AS_depth)

            # add cell-independent features here
            if node.genes or node.proteins or node.proteoforms or node.lipids or node.metabolites or node.ftu or node.references:
                # skip the cell block
                record = anatomical_structure + "\t\t\t"
                # add features
                add_features(record, node)

            # add cells and cell-dependent features
            cells = node.cells
            if cells:
                for cell in cells:
                    cell_feature = features[cell]
                    record = anatomical_structure + get_data_block(cell_feature) + "\t"

                    # add cell-specific features
                    add_features(record, cell_feature)

                    
# Prepare to quit
def close_files():
    global in_file, out_file

    in_file.close()
    out_file.close()


# Quit with an error
def exit_with_error(error):
    close_files()
    sys.exit(error)


# Convert the string to an array and also track the longest array, so
# we know how many levels for the feature type.
def split_string(array_string, max_depth):
    tmp = []
    if array_string:
        tmp = array_string.split(',')
        if max_depth < len(tmp):
            max_depth = len(tmp)
    return tmp, max_depth


################################################
# Main

# Cxecute script
if __name__ == "__main__":
    # open input file. argv[0] is the program name

    in_filename = sys.argv[1]
    in_file = open(in_filename, "r")

    # this will overwrite any existing file
    out_filename = sys.argv[2]
    out_file = open(out_filename, "w")

    contents = in_file.readlines()
    headerLine = True
    for line in contents:
        # XXX: need to catch if this throws an exception due to incorrect number of values
        name, label, id, feature_type, children_string, cells_string, genes_string, proteins_string, proteoforms_string, lipids_string, metabolites_string, ftu_string, references_string = re.split(r'\t', line.rstrip('\n'))

        # the first line is a header line, which we skip here
        if headerLine:
            headerLine = False
            continue

        # remove the quotes and extra spaces from the TSV file
        name = name.rstrip()
        children_string = children_string.replace('"', '').replace(', ', ',')
        cells_string = cells_string.replace('"', '').replace(', ', ',')
        genes_string = genes_string.replace('"', '').replace(', ', ',')
        proteins_string = proteins_string.replace('"', '').replace(', ', ',')
        proteoforms_string = proteoforms_string.replace('"', '').replace(', ', ',')
        lipids_string = lipids_string.replace('"', '').replace(', ', ',')
        metabolites_string = metabolites_string.replace('"', '').replace(', ', ',')
        ftu_string = ftu_string.replace('"', '').replace(', ', ',')
        references_string = references_string.replace('"', '').replace(', ', ',')

        # convert strings into lists and, if relevant, get the number
        # of levels per feature type
        children = []
        if len(children_string) > 0:
            children = children_string.split(',')
        cells = []
        if len(cells_string) > 0:
            cells = cells_string.split(',')
        genes, max_genes = split_string(genes_string, max_genes)
        proteins, max_proteins = split_string(proteins_string, max_proteins)
        proteoforms, max_proteoforms = split_string(proteoforms_string, max_proteoforms)
        lipids, max_lipids = split_string(lipids_string, max_lipids)
        metabolites, max_metabolites = split_string(metabolites_string, max_metabolites)
        ftu, max_ftu = split_string(ftu_string, max_ftu)
        references, max_references = split_string(references_string, max_references)

        '''
        # TEST: make sure biomarkers and references are only applied to leaves or cell types
        if (feature_type == "AS") and children and (genes, lipids...):
            error = "ERROR: biomarkers and references can only be applied to anatomical structures without children or to cell types."
            error += "\n\tStructure: " + name
            error += "\n\tChildren: " + str(children)
            error += "\n\tCells: " + str(cells)
            error += "\n\tGenes: " + str(genes)
            exit_with_error(error)
        '''
        
        # TEST: check that biomarkers and references are only applied to AS and CT.
        if (feature_type not in ("AS", "CT")) and (genes or proteins or proteoforms or lipids or metabolites or ftu or references):
            error = "ERROR: biomarkers and references can only be applied to anatomical structures or cell types."
            error += "\n\tEntity: " + name
            error += "\n\tLabel: " + label
            error += "\n\tID: " + id
            error += "\n\tChildren: " + str(children)
            error += "\n\tCells: " + str(cells)
            error += "\n\tGenes: " + str(genes)
            error += "\n\tProteins: " + str(proteins)
            error += "\n\tProteoforms: " + str(proteoforms)
            error += "\n\tLipids: " + str(lipids)
            error += "\n\tMetabolites: " + str(metabolites)
            error += "\n\tFTU: " + str(ftu)
            error += "\n\tReferences: " + str(references)
            exit_with_error(error)

        # add record to the node or feature dictionaries
        if feature_type == "AS":
            # "children" is a reserved word in Node() and here we want
            # to store the string of children, so we use "kids"
            node = Node(name, label=label, id=id, kids=children, cells=cells, genes=genes, proteins=proteins, proteoforms=proteoforms, lipids=lipids, metabolites=metabolites, ftu=ftu, references=references)
            if name in nodes:
                # TEST: make sure all anatomical structures are uniquely named
                error = "ERROR: two anatomical structures with the same name."
                error += "\n\t" + str(node)
                error += "\n\t" + str(nodes[name])
                exit_with_error(error)
            else:
                nodes[name] = node
        else:
            # Feature class is used to store both cells and biomarkers.
            feature = Feature(name, feature_type, label=label, id=id, genes=genes, proteins=proteins, proteoforms=proteoforms, lipids=lipids, metabolites=metabolites, ftu=ftu, references=references)
            if name in features:
                # TEST: make sure all features are uniquely named
                error = "ERROR: two features with the same name."
                error += "\n\t" + str(feature)
                error += "\n\t" + str(features[name])
                exit_with_error(error)
            else:
                features[name] = feature

    # everything loaded, so now build the tree
    for key in nodes.keys():
        node = nodes[key]
        children = node.kids
        if children:
            for child in children:
                if child not in nodes:
                    # TEST: make sure the child exists
                    error = "ERROR: anatomical structure hasn't been defined."
                    error += "\n\tStructure: " + child
                    exit_with_error(error)
                child_node = nodes[child]
                if child_node.parent:
                    # TEST: make sure no nodes have multiple parents
                    error = "ERROR: anatomical structure has two parents."
                    error += "\n\tStructure: " + child
                    error += "\n\tParent: " + child_node.parent.name
                    error += "\n\tParent: " + key
                    exit_with_error(error)
                child_node.parent = node

    # TEST: make sure only one root node
    for key in nodes.keys():
        if nodes[key].is_root:
            if tree_root:
                error = "ERROR: multiple anatomical structures lacking parents."
                error += "\n\tStructure: " + tree_root.name
                error += "\n\tStructure: " + key
                exit_with_error(error)
            tree_root = nodes[key]

    # print the ASCT+B table
    print_ASCTB_table()

    if DEBUG:
        # Height returns the number of edges, but we want the number of node levels.
        print("Anatomical structure levels: ", tree_root.height+1)
        # print each of the levels.
        # print([[node.name for node in children] for children in LevelOrderGroupIter(tree_root)])

        # print tree (simplified and with all details)
        print("\n", RenderTree(tree_root, style=AsciiStyle()).by_attr())

        print("\nNodes:")
        for node in nodes:
            print(nodes[node])

        print("\nFeatures:")
        for feature in features:
            print(features[feature])

    close_files()
