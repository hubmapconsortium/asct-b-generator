#!/usr/bin/python3

import sys
import argparse
import csv
from typing import List, Tuple, Dict, TextIO
from anytree import Node, SymlinkNode, RenderTree, AsciiStyle, LevelOrderGroupIter, Walker
from anytree.exporter import DotExporter

"""
TODO:
- need better demo and more clear instructions. Use color coding in demo and documentation.
- use a real world demo, perhaps defining a family tree.
- label = Uberon RDFS label
- "input table" is the input file for the Generator
- note in docs that revisions need to happen to the input table, not the ASCT+B table
- expect 10 row header and duplicate into output
"""

'''
Requires anytree:
https://pypi.org/project/anytree/

Assumptions:

  1. Inner nodes can contain biomarkers or references, not just leaves
  and cell types.

  2. All anatomical structures must be uniquely named, for example,
     there can not be two structures called "ovary" but there can be
     "left ovary" and "right ovary".

  3. Cell type is only one level; that is, cell types can not have
  child cell types.

  4. Commas can not be used in names for anatomical structures, cells,
  or features.
'''
################################################
# Classes

class Feature:
    """
    Used to store cells as well as biomarkers and references
    """
    def __init__(self, name: str, feature_type: str, label="", id="", note="",\
                abbr="", **kwargs)->None:
        self.name = name
        self.feature_type = feature_type
        self.label = label
        self.id = id
        self.note = note
        self.abbr = abbr
        self.genes = kwargs.get('genes', None)
        self.proteins = kwargs.get('proteins', None)
        self.proteoforms = kwargs.get('proteoforms', None)
        self.lipids = kwargs.get('lipids', None)
        self.metabolites = kwargs.get('metabolites', None)
        self.ftu = kwargs.get('ftu', None)
        self.references = kwargs.get('references', None)

    def __str__(self)->str:
        out = "[" + self.name
        out += " type:" + self.feature_type
        if self.genes:
            out += " genes:" + str(self.genes)
        out += "]"
        return out

################################################
# Functions

def no_features(node: Node)->bool:
    """
    Check if a node contains any cells or features.
    """
    # check for all biomarkers and references
    if node.cells or node.genes or node.proteins or node.proteoforms \
        or node.lipids or node.metabolites or node.ftu or node.references:
        return False
    return True

def get_header_block(content: str, depth: int)->List[str]:
    """
    Generates the five column headers as exemplified below:
    AS/1	AS/1/LABEL	AS/1/ID	AS/1/NOTE	AS/1/ABBR
    """
    headers = [content + "/" + str(depth), content + "/" + str(depth) + "/LABEL", \
    content + "/" + str(depth) + "/ID", content + "/" + str(depth) + "/NOTE", \
    content + "/" + str(depth) + "/ABBR"]
    return headers

def get_biomarker_header(header: List[str])->List[str]:
    """
    Takes the anatomical structures header list as input
    and appends biomarker headers to the same input list
    """
    biomarker_headers = []
    biomarkers_list = ["BGene", "BProtein", "BProteoform", "BLipid", "BMetabolites", "FTU"]
    depth_list = [max_genes, max_proteins, max_proteoforms, max_lipids, max_metabolites, max_ftu]
    for i, _ in enumerate(biomarkers_list):
        for depth in range(1, depth_list[i]+1):
            headers = get_header_block(biomarkers_list[i], depth)
            biomarker_headers.append(headers)
    for biomarkers in biomarker_headers:
        for biomarker in biomarkers:
            header.append(biomarker)
    return header

def print_column_header(max_as_depth: int)->List[str]:
    """
    Returns the column headers according to the depth of anatomical structure tree.
    """
    header = []
    col_name_list = []
    # add anatomical structures
    for depth in range(1, max_as_depth+1):
        col_name_list = get_header_block("AS", depth)
        for as_header in col_name_list:
            header.append(as_header)
    col_name_list = []
    # assume cell type is only one depth
    col_name_list = get_header_block("CT", 1)
    for cell_headers in col_name_list:
        header.append(cell_headers)
    col_name_list = []
    # add the headers for biomarkers
    get_biomarker_header(header)
    # we handle references separately because their header details
    # differ from the biomarkers.
    for depth in range(1, max_references+1):
        col_name_list.append("REF/" + str(depth))
        col_name_list.append("REF/" + str(depth) + "/DOI")
        col_name_list.append("REF/" + str(depth) + "/NOTES")
    for reference in col_name_list:
        header.append(reference)
    return header

def get_data_block(feature: Feature)->List[str]:
    """
    Returns a set of values for a non-reference feature.
    """
    data_block_list = [feature.name, feature.label, feature.id, feature.note, feature.abbr]
    return data_block_list

def get_reference_block(feature: Feature)->List[str]:
    """
    Returns a set of values for a reference.
    References only have 3 columns where as biomarkers contain 5 columns.
    """
    reference_block_list = [feature.name, feature.label, feature.id]
    return reference_block_list

def add_biomarkers(
        elements: List[str], max_depth: int, is_reference: bool, in_file: TextIO,\
        out_file: TextIO, dot_file: TextIO, missing_feature_ok: bool,\
        features: Dict[str, 'Feature'])->List[str]:
    """
    Returns the biomarker data and/or empty cells as appropriate
    """
    biomarkers = []
    count = 0
    if elements:
        for element in elements:
            if element not in features:
                # TEST: feature isn't defined
                if missing_feature_ok:
                    # add a place holder cell type
                    features[element] = Feature(element, "", label="", id="", note="", abbr="")
                else:
                    error = "ERROR: feature hasn't been defined. Please add a row to the input that\
                            defines this feature."
                    error += "\n\tFeature: " + element
                    exit_with_error(error, in_file, out_file, dot_file)
            if is_reference:
                reference_list = get_reference_block(features[element])
                if reference_list:
                    for reference in reference_list:
                        biomarkers.append(reference)
            else:
                data_block_list = get_data_block(features[element])
                if data_block_list:
                    for data_block in data_block_list:
                        biomarkers.append(data_block)
            count += 1
    if count < max_depth:
        for _ in range((max_depth-count)):
            for _ in range(5):
                biomarkers.append("")
    return biomarkers

def add_features(
        record: list, element: 'Feature', in_file: TextIO, out_file: TextIO, dot_file: TextIO,
        missing_feature_ok: bool, features: Dict[str, 'Feature'])->List[str]:
    """
    The argument (element) can be either an anatomical structure ("node")
    or a cell as they contain the same features.
    """
    # need to track how many entities we output per feature type, to
    # appropriately pad the output.
    # add biomarkers
    all_biomarkers = [
        add_biomarkers(element.genes, max_genes, False, in_file, out_file, dot_file,\
            missing_feature_ok, features),
        add_biomarkers(element.proteins, max_proteins, False, in_file, out_file,\
            dot_file, missing_feature_ok, features),
        add_biomarkers(element.proteoforms, max_proteoforms, False, in_file, out_file,\
            dot_file, missing_feature_ok, features),
        add_biomarkers(element.lipids, max_lipids, False, in_file, out_file, dot_file,\
            missing_feature_ok, features),
        add_biomarkers(element.metabolites, max_metabolites, False, in_file, out_file,\
            dot_file, missing_feature_ok, features),
        add_biomarkers(element.ftu, max_ftu, False, in_file, out_file, missing_feature_ok,\
            dot_file, features),
        add_biomarkers(element.references, max_references, True, in_file, out_file,\
            dot_file, missing_feature_ok, features)]
    for bio_markers in all_biomarkers:
        for bio_marker in bio_markers:
            record.append(bio_marker)
    return record

def print_asctb_table(
        in_file: TextIO, out_file: TextIO, dot_file: TextIO,\
        missing_feature_ok: bool, nodes: Dict[str, Node],\
        features: Dict[str, 'Feature'], root: Node)->Tuple[List[str], List[str]]:
    """
    Return the ASCT+B table data, stepping through the tree, one level at a
    time and applying each of the relevant cells and features to the
    node.
    """
    features_data = []
    # get a list of lists where the sub-lists are the nodes per level
    # of the tree. The list returned will appropriately include
    # Symlink'ed Nodes but it won't actually reference the Symlink but
    # rather the original Node. So we need to manually process the
    # Symlinks, when we've seen a Node more than once.
    levels = [[node.name for node in children] for children in LevelOrderGroupIter(root)]
    # total number of anatomical structure levels
    max_as_depth = len(levels)
    column_headers = print_column_header(max_as_depth)
    #print column headers in main

    walker = Walker()
    for level in levels:
        for node_name in level:
            node = nodes[node_name]
            # if we've already processed this node once, then we need
            # to use one of the Symlinks. LevelOrderGroupIter() only
            # returns the original Node in place of each of the
            # Symlinks.
            if node.processed > 0:
                node = node.symlinks[node.processed-1]
            node.processed += 1
            # if inner node and doesn't have any features, then skip
            # to the next node.
            if no_features(node) and not node.is_leaf:
                continue

            # walk the tree up from the node to the root. The first
            # tuple is for upward walking and the second is just the
            # tree root. So we only care about the third tuple which
            # is the downward walking.
            walked = walker.walk(root, node)[2]
            # add root, since it was in the middle tuple that we
            # ignore.
            walked = (root,) + walked

            # compute the AS level (ie depth of tree), so we can pad
            # as needed when we add cells.
            as_depth = len(walked)

            # add all ancestral anatomical structures to the data
            # record being output.
            anatomical_structure = []
            for ancestor in walked:
                data_block_list = get_data_block(ancestor)
                for data_block in data_block_list:
                    anatomical_structure.append(data_block)
                anatomical_structure_data = anatomical_structure[:]
            # pad record, to account for structures that aren't as
            # deep as the maximum possible depth (max_as_depth).
            for _ in range((max_as_depth - as_depth)):
                for _ in range(5):
                    anatomical_structure.append("")
                anatomical_structure_data = anatomical_structure[:]
            # track if we've printed the node in some form
            node_output = False

            # add cell-independent features here
            if node.genes or node.proteins or node.proteoforms or node.lipids or node.metabolites\
            or node.ftu or node.references:
                # skip the cell block
                for _ in range(5):
                    anatomical_structure_data.append("")
                # add features
                features_data.append(add_features(anatomical_structure_data, node, in_file, out_file,\
                                            dot_file, missing_feature_ok, features))
                #print featuresdata in main
                node_output = True

            # add cells and cell-dependent features
            cells = node.cells
            if cells:
                for cell in cells:
                    if cell not in features:
                        # TEST: cell isn't defined
                        if missing_feature_ok:
                            # add a place holder cell type
                            features[cell] = Feature(cell, "CT", label="", id="", note="", abbr="")
                        else:
                            error = "ERROR: cell hasn't been defined. Please add a row to the input\
                                    that defines this cell."
                            error += "\n\tCell: " + cell
                            error += "\n\tAnatomical Structure: " + node.name
                            exit_with_error(error, in_file, out_file, dot_file)

                    cell_feature = features[cell]
                    cell_feature_data_block_list = get_data_block(cell_feature)
                    for cell_data in cell_feature_data_block_list:
                        anatomical_structure_data.append(cell_data)
                    #add cell-specific features
                    features_data.append(add_features(anatomical_structure_data, cell_feature,\
                                                    in_file, out_file, dot_file,\
                                                    missing_feature_ok, features))
                    node_output = True
                    anatomical_structure_data = anatomical_structure[:]

            # if node is leaf and doesn't have any assigned cells or
            # features then we output it here.
            if not node_output:
                for _ in range((1 + max_genes + max_proteins + max_proteoforms + max_lipids \
                                + max_metabolites + max_ftu + max_references)):
                    anatomical_structure.append("")
                features_data.append(anatomical_structure)
    return column_headers, features_data

def close_files(in_file: TextIO, out_file: TextIO, dot_file: TextIO)->None:
    """
    Prepare to quit
    """
    in_file.close()
    out_file.close()
    if dot_file:
        dot_file.close()

def exit_with_error(error: str, in_file: TextIO, out_file: TextIO, dot_file: TextIO)->None:
    """
    Quit with an error
    """
    close_files(in_file, out_file, dot_file)
    sys.exit(error)

def process_input(input_string: str, max_depth: int)->Tuple[List[str], int]:
    """
    Clean up the input, convert it to an array and compute the longest
    array, per feature type.
    """

    # remove the quotes and extra spaces from the input string
    input_string = input_string.replace('"', '').replace(', ', ',').strip()

    # convert the string to an array and also track the longest array, so
    # we know how many levels for the feature type.
    tmp = []
    if input_string:
        tmp = input_string.split(',')
        if max_depth < len(tmp):
            max_depth = len(tmp)

    # return the array and the depth
    return tmp, max_depth

def build_tree(
        nodes_to_process: Dict[str, Node], dup_nodes: dict, in_file: TextIO, out_file: TextIO,
        dot_file: TextIO, only_one_parent: bool, nodes: Dict[str, Node])->Dict[Node, Node]:
    """
    Builds a tree
    """
    for node in nodes_to_process:
        children = nodes[node].kids
        if children:
            for child in children:
                if child not in nodes:
                    error = "ERROR: anatomical structure hasn't been defined. Please add a row to \
                    the input that defines this structure."
                    error += "\n\tStructure: " + child
                    error += "\n\tParent: " + node
                    exit_with_error(error, in_file, out_file, dot_file)

                child_node = nodes[child]
                if child_node.parent:
                    # TEST: a node has more than one parent
                    if only_one_parent:
                        error = "ERROR: anatomical structure has two parents."
                        error += "\n\tStructure: " + child
                        error += "\n\tParent: " + child_node.parent.name
                        error += "\n\tParent: " + node
                        exit_with_error(error, in_file, out_file, dot_file)

                    # child already has a parent, so need to create a
                    # duplicate (Symlink) child, to allow for multiple
                    # parents.
                    new_child = SymlinkNode(child_node, parent=nodes[node])
                    dup_nodes[new_child] = new_child
                    # We need to track duplicates because
                    # LevelOrderGroupIter() doesn't differentiate
                    # between a Node and a Symlink. Hence when we
                    # export the tree Symlinks get lost.
                    if child_node.symlinks:
                        child_node.symlinks.append(new_child)
                    else:
                        child_node.symlinks = [new_child]
                else:
                    child_node.parent = nodes[node]
    return dup_nodes

def process_arguments()->Tuple[TextIO, TextIO, bool, bool, bool]:
    """
    Handle command line arguments.
    """
    # the file handles
    dot_file = None

    parser = argparse.ArgumentParser(description="Generate ASCT+B table.")
    parser.add_argument("-m", "--missing", help="Ignore missing cell types, biomarkers and\
                        references. For example, if a cell type is marked as containing a biomarker that wasn't defined,\
                        this flag would prevent the program from exiting with an error and instead the ASCT+B table would\
                        be generated. When the flag isn't used, all features must be defined.",\
                        action="store_true")
    parser.add_argument("-u", "--unique", help="Make sure all anatomical structures have one and,\
                        only one parent.", action="store_true")
    parser.add_argument("-d", "--dot", help="Output tree as a DOT file for plotting with Graphviz."\
                        , action="store_true")
    parser.add_argument("-v", "--verbose", help="Print the tree to the terminal."\
                        , action="store_true")
    parser.add_argument("input", type=str, help="Input file")
    parser.add_argument("output", type=str, help="Output file (CSV)")

    args = parser.parse_args()

    # open input file and check the delimiter in the file.
    #return the reader object to access it in the main function
    in_file = open(args.input, "r")
    dialect = csv.Sniffer().sniff(in_file.read(1024))
    in_file.seek(0)
    if dialect.delimiter == ',':
        file_reader = csv.reader(in_file, delimiter=',')
    elif dialect.delimiter == '\t':
        file_reader = csv.reader(in_file, delimiter='\t')

    # this will overwrite any existing file
    out_file = open(args.output, "w", newline='')

    if args.dot:
        print("indotfile")
        dot_file = open(args.output + ".dot", "w")

    print_tree = args.verbose
    only_one_parent = args.unique
    missing_feature_ok = args.missing
    return in_file, out_file, dot_file, print_tree, only_one_parent, missing_feature_ok, file_reader

################################################
# Main

# Cxecute script
if __name__ == "__main__":

    DEBUG = False

    # number of columns in the input file
    INPUT_COLUMNS: int = 15

    # print the tree to the command line
    PRINT_TREE: bool = False

    # requiring each anatomical structure to have one and only one
    # parent. When true, if the same sub-structure exists in multiple
    # parent-structures, then each child structure would need to be
    # uniquely named.
    SINGLE_PARENT: bool = False

    # if true, we will automatically generate missing cells, biomarkers
    # and references. The generated ones will contain no secondary
    # details. Anatomical structures must be defined.
    IS_MISSING_FEATURE: bool = False

    # the top level node
    tree_root: Node = ""

    # all cells, biomarkers and references
    all_features: Dict[str, 'Feature'] = {}

    # compute these based on the max number of biomarkers across nodes and
    # cells
    max_genes: int = 0
    max_proteins: int = 0
    max_proteoforms: int = 0
    max_lipids: int = 0
    max_metabolites: int = 0
    max_ftu: int = 0
    # we don't actually need to track the max level for references as it's
    # last data block, however, this is a hack to allow us to handle
    # references like the biomarkers.
    max_references: int = 0

    INPUT_FILE, OUTPUT_FILE, DOT_FILE, PRINT_TREE, SINGLE_PARENT,\
    IS_MISSING_FEATURE, CONTENTS = process_arguments()
    CSVWRITER = csv.writer(OUTPUT_FILE)
    LINECOUNT = 0

    # number of lines expected as the header. This row count includes the
    # descriptive text at the top of the input file and also the row of
    # column-specific headers
    HEADER_LENGTH = 11

    # all anatomical structures
    anatomical_structure_nodes: Dict[str, Node] = {}

    for line in CONTENTS:
        if LINECOUNT < HEADER_LENGTH:
            if LINECOUNT < HEADER_LENGTH-1:
                CSVWRITER.writerow(line)
            LINECOUNT += 1
            continue

        line_as_list = line
        # make sure the line contains the appropriate number of fields
        if len(line_as_list) != INPUT_COLUMNS:
            error_message = "ERROR: incorrect number of fields in line. The tab-delimited line\
                    should contain the following " + str(INPUT_COLUMNS) + " fields: "
            error_message += "\n\tname, label, ID, node, abbreviation, feature type, children,\
                    cells, genes, proteins, proteoforms, lipids, metabolites, FTU, references"
            error_message += "\n\tNumber of fields found in line: " + str(len(line_as_list))
            error_message += "\n\tLine: " + str(line)
            exit_with_error(error_message, INPUT_FILE, OUTPUT_FILE, DOT_FILE)
        names, labels, id_, notes, abbrs, feature_types, children_string, cells_string,\
        genes_string, proteins_string, proteoforms_string, lipids_string, metabolites_string,\
        ftu_string, references_string = line_as_list
        # clean up white spaces
        names = names.strip()

        # TEST: make sure names don't contain commas
        if names.find(',') > -1:
            error_message = "ERROR: names can not include commas (',')."
            error_message += "\n\tName: " + names
            exit_with_error(error_message, INPUT_FILE, OUTPUT_FILE, DOT_FILE)

        # convert strings into lists and get the number of levels per
        # feature type
        children_string = children_string.replace('"', '').replace(', ', ',').strip()
        children_list = []
        if children_string:
            children_list = children_string.split(',')
        cells_string = cells_string.replace('"', '').replace(', ', ',').strip()
        cells_list = []
        if cells_string:
            cells_list = cells_string.split(',')
        genes, max_genes = process_input(genes_string, max_genes)
        proteins, max_proteins = process_input(proteins_string, max_proteins)
        proteoforms, max_proteoforms = process_input(proteoforms_string, max_proteoforms)
        lipids, max_lipids = process_input(lipids_string, max_lipids)
        metabolites, max_metabolites = process_input(metabolites_string, max_metabolites)
        ftu, max_ftu = process_input(ftu_string, max_ftu)
        references, max_references = process_input(references_string, max_references)

        '''
        # TEST: make sure biomarkers and references are only applied to leaves or cell types
        if (feature_type == "AS") and children and (genes, lipids...):
            error_message = "ERROR: biomarkers and references can only be applied to anatomical\
                    structures without children or to cell types."
            error_message += "\n\tStructure: " + name
            error_message += "\n\tChildren: " + str(children)
            error_message += "\n\tCells: " + str(cells)
            error_message += "\n\tGenes: " + str(genes)
            exit_with_error(error_message,in_file,out_file)
        '''

        # TEST: check that biomarkers and references are only applied to AS and CT.
        if (feature_types not in ("AS", "CT")) and (genes or proteins or proteoforms or lipids or\
            metabolites or ftu or references):
            error_message = "ERROR: biomarkers and references can only be applied to anatomical\
            structures or cell types."
            error_message += "\n\tEntity: " + names
            error_message += "\n\tLabel: " + labels
            error_message += "\n\tID: " + id_
            error_message += "\n\tnote: " + notes
            error_message += "\n\tABBR: " + abbrs
            error_message += "\n\tChildren: " + str(children_list)
            error_message += "\n\tCells: " + str(cells_list)
            error_message += "\n\tGenes: " + str(genes)
            error_message += "\n\tProteins: " + str(proteins)
            error_message += "\n\tProteoforms: " + str(proteoforms)
            error_message += "\n\tLipids: " + str(lipids)
            error_message += "\n\tMetabolites: " + str(metabolites)
            error_message += "\n\tFTU: " + str(ftu)
            error_message += "\n\tReferences: " + str(references)
            exit_with_error(error_message, INPUT_FILE, OUTPUT_FILE, DOT_FILE)

        # add record to the node or feature dictionaries
        if feature_types == "AS":
            # "children" is a reserved word in Node() and here we want
            # to store the string of children, so we use "kids"
            as_node = Node(names, label=labels, id=id_, note=notes, abbr=abbrs, kids=children_list,\
                        cells=cells_list, genes=genes, proteins=proteins, proteoforms=proteoforms,\
                        lipids=lipids, metabolites=metabolites, ftu=ftu, references=references,\
                        symlinks=None, processed=0)

            if names in anatomical_structure_nodes:
                # TEST: anatomical structures' name already used
                error_message = "ERROR: two anatomical structures can not have the same name.\
                        However, by default, a single anatomical structure can be the child of\
                        multiple parents."
                error_message += "\n\t" + str(as_node)
                error_message += "\n\t" + str(anatomical_structure_nodes[names])
                exit_with_error(error_message, INPUT_FILE, OUTPUT_FILE, DOT_FILE)
            else:
                anatomical_structure_nodes[names] = as_node
        else:
            # Feature class is used to store both cells and biomarkers.
            feature_class = Feature(names, feature_types, label=labels, id=id_, note=notes,\
                            abbr=abbrs, genes=genes, proteins=proteins, proteoforms=proteoforms,\
                            lipids=lipids, metabolites=metabolites, ftu=ftu, references=references)
            if names in all_features:
                # TEST: feature name already used
                error_message = "ERROR: all features must have a unique name."
                error_message += "\n\t" + str(feature_class)
                error_message += "\n\t" + str(all_features[names])
                exit_with_error(error_message, INPUT_FILE, OUTPUT_FILE, DOT_FILE)
            else:
                all_features[names] = feature_class

    # everything loaded, so now build the tree, allowing for duplicate nodes
    DUPLICATE_NODES = build_tree(anatomical_structure_nodes, {}, INPUT_FILE, OUTPUT_FILE, DOT_FILE,\
                                SINGLE_PARENT, anatomical_structure_nodes)
    # add duplicate nodes to our global node dictionary
    anatomical_structure_nodes.update(DUPLICATE_NODES)

    # keep running build_tree() until all duplicate nodes have been
    # processed as each run of build_tree() might generate more
    # duplicate nodes needing to be processed.
    while DUPLICATE_NODES:
        DUPLICATE_NODES = build_tree(DUPLICATE_NODES, {}, INPUT_FILE, OUTPUT_FILE, DOT_FILE,\
                                    SINGLE_PARENT, anatomical_structure_nodes)
        anatomical_structure_nodes.update(DUPLICATE_NODES)

    # set tree_root to the root of the tree
    for as_node in anatomical_structure_nodes:
        if anatomical_structure_nodes[as_node].is_root:
            if tree_root:
                # TEST: more than one node without a parent (ie., "root").
                error_message = "ERROR: anatomical structure(s) lacking parents."
                error_message += "\n\tStructure: " + tree_root.name
                error_message += "\n\tStructure: " + as_node
                exit_with_error(error_message, INPUT_FILE, OUTPUT_FILE, DOT_FILE)
            tree_root = anatomical_structure_nodes[as_node]

    # print the ASCT+B table
    COLUMN_HEADER, FEATURE_DATA = print_asctb_table(INPUT_FILE, OUTPUT_FILE, DOT_FILE,\
                                                    IS_MISSING_FEATURE, anatomical_structure_nodes,\
                                                    all_features, tree_root)
    CSVWRITER.writerow(COLUMN_HEADER)
    for feature_row in FEATURE_DATA:
        CSVWRITER.writerow(feature_row)

    if DEBUG:
        # Height returns the number of edges, but we want the number of node levels.
        print("Anatomical structure levels: ", tree_root.height+1)
        # print each of the levels.
        print([[node.name for node in children] for children in LevelOrderGroupIter(tree_root)])

        # print tree with all details
        print("\n", RenderTree(tree_root))

        print("\nNodes:")
        for as_node in anatomical_structure_nodes:
            if anatomical_structure_nodes[as_node].parent:
                print(as_node, anatomical_structure_nodes[as_node].parent.name)

        print("\nFeatures:")
        for feature_ in all_features:
            print(all_features[feature_])

    if PRINT_TREE:
        print("\n", RenderTree(tree_root, style=AsciiStyle()).by_attr())

    if DOT_FILE:
        for line in DotExporter(tree_root):
            print("End line is:", line)
            #dot_file.write(line)

    close_files(INPUT_FILE, OUTPUT_FILE, DOT_FILE)
