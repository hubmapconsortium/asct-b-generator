# Changelog

Changelog for the ASCT+B Generator

## 1.1.1 - 2022-08-08

* Fixed alignment issues for the data with cell-independent features.
* Fixed alignment issues where the node is leaf and doesn't have any assigned cells or features.
* Fixed computing the depth of the anatomical structure.

## 1.1.0 - 2022-06-25

* Added a feature to generate CSV files from TSV/CSV files.
* Added typings and docstrings.
* Modified the code according to PEP-8 style guide.

## 1.0.3 - 2022-05-25

### Added in 1.0.3

* Fixed a bug causing each reference to incorrectly span 5 columns. References now correctly span 3 columns.

## 1.0.2 - 2022-02-20

### Added in 1.0.2

* Added "note" and "ABBR" (abbreviation) fields for each feature. These columns are now required in the input file.
* The first 10 lines of the input file are now assumed to contain a descriptive header. This is duplicated to the output file. The 11th line in the input file is assumed to be a per-column header that is ignored.

## 1.0.1 - 2022-02-15: 

### Added in 1.0.1

* Test to be sure entities (anatomical structures, cell types, biomarkers and references) don't have commas in their names. 
* Added a command line argument users can set to force anatomical structures to be unique (ie., each has one and only one parent). When set, if the same sub-structure exists in multiple parent-structures, then each child structure would need to be uniquely named. By default, a child structure can have multiple parent structures.
* Added a command line argument users can set to cause the program to automatically create missing features, when features are used. For example, if a biomarker is assigned to a cell type, by default, the biomarker must be independently defined but now users can optionally disable this requirement. Anatomical structures must always be defined.

## 1.0 - 2022-02-14

### Added in 1.0

* Removed the (erroneous) assumption that an anatomical structure can only have a single parent, added more validation of the inputs, added debugging output options, and better handle command line arguments. Also various bug fixes.
* These are the significant differences from version v0.1 and v1.0.

    1. The command line arguments have been greatly simplified.
    2. The number of AS levels is automatically computed.
    3. The input file has changed with this release. Cells are listed in a separate column from the children column.
    4. Biomarkers and references can now be added to any anatomical structure.
    5. The anytree Python module is required.
    6. A header is autogenerate in the output file.
    7. A DOT file can be generated to display the tree in Graphviz.
    8. Lots of tests to validate the input file

## 1.0-beta - 2022-02-11

### Added in 1.0-beta

* This is a complete rewrite of the program. This version has none of the limitations from the alpha version, it includes more data validation, and requires less user intervention. 

## 0.1-alpha - 2022-02-02

### Added in 0.1-alpha

* This is a proof of concept and not meant for production use.