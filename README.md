# ASCT+B Generator

This program converts a simple TSV file into a HuBMAP ASCT+B table.

The included file "demo-input.txt" was generated by Excel using the "demo-input.xlsx" file (Save As "Tab delimited Text"). The generated output is a TSV file. The example file "demo-output.xls" was generated by the program.

## Version

May 25, 2022

- Fixed a bug causing each reference to incorrectly span 5 columns. References now correctly span 3 columns.

Feb 20, 2022

- Added "note" and "ABBR" (abbreviation) fields for each feature. These columns are now required in the input file.
- The first 10 lines of the input file are now assumed to contain a descriptive header. This is duplicated to the output file. The 11th line in the input file is assumed to be a per-column header that is ignored.

Feb 15, 2022: 

- Test to be sure entities (anatomical structures, cell types, biomarkers and references) don't have commas in their names. 

- Added a command line argument users can set to force anatomical structures to be unique (ie., each has one and only one parent). When set, if the same sub-structure exists in multiple parent-structures, then each child structure would need to be uniquely named. By default, a child structure can have multiple parent structures.
- Added a command line argument users can set to cause the program to automatically create missing features, when features are used. For example, if a biomarker is assigned to a cell type, by default, the biomarker must be independently defined but now users can optionally disable this requirement. Anatomical structures must always be defined.

**v1.0** - Feb 14, 2022: Removed the (erroneous) assumption that an anatomical structure can only have a single parent, added more validation of the inputs, added debugging output options, and better handle command line arguments. Also various bug fixes.

These are the significant differences from version v0.1 and v1.0.

    1. The command line arguments have been greatly simplified.
    2. The number of AS levels is automatically computed.
    3. The input file has changed with this release. Cells are listed in a separate column from the children column.
    4. Biomarkers and references can now be added to any anatomical structure.
    5. The anytree Python module is required.
    6. A header is autogenerate in the output file.
    7. A DOT file can be generated to display the tree in Graphviz.
    8. Lots of tests to validate the input file

**v1.0-beta** - Feb 11, 2022: This is a complete rewrite of the program. This version has none of the limitations from the alpha version, it includes more data validation, and requires less user intervention. 

**v0.1-alpha** - Feb 2, 2022: This is a proof of concept and not meant for production use.

## Assumptions

The following assumptions are built into the program.
  1. The ASCT+B table format allows anatomical structures that are not "leaves" to contain biomarkers or references.
  2. All anatomical structures must be uniquely named, for example, there can not be two structures called "ovary" but there can be "left ovary" and "right ovary".
  3. Cell type is only one level.
  3. Commas can not be used in names for anatomical structures, cells, or features.
  3. ***It is assumed that the "author preferred name" is unique across anatomical structures and ontology IDs.***

## Data validation

The program performs the following data validation checks.
  1. Check that there is only one root to the anatomical structure.
  1. Enforce the parent requirements for anatomical structure. By default an anatomical structure can have multiple parents. For example, the primary ovarian follicle and the primordial ovarian follicle both have a granulosa cell layer. A command line argument can change this behavior such that anatomical structures can have only one parent.
  2. Check that anatomical structures, cells, biomarkers and references are appropriately defined. By default the program requires all features be explicitely defined, although a command line argument can disable this requirement.
  3. Check that anatomical structures, cell types, biomarkers and references all have unique names.
  3. Check that names do not contain commas.
  3. Check that biomarkers and references are only applied to anatomical structures and cell types.

## Requirements

This program has only been tested on a Mac OS using Python 3. Although it should work on a Linux system.

The program requires the anytree Python package.
```
https://pypi.org/project/anytree/

```
The anytree package can be installed as follows.
```
python3 -m pip install anytree --user
```

## Usage

```
usage: process.py [-h] [-m] [-u] [-d] [-v] input output

Generate ASCT+B table.

positional arguments:
  input          Input file
  output         Output file (TSV)

optional arguments:
  -h, --help     show this help message and exit
  -m, --missing  Ignore missing cell types, biomarkers and references. For example, if a cell type is marked as containing a biomarker that wasn't defined, this flag would prevent the program from exiting with an error and instead the ASCT+B table would be generated. When the flag isn't used, all features must be defined.
  -u, --unique   Make sure all anatomical structures have one and only one parent.
  -d, --dot      Output tree as a DOT file for plotting with Graphviz.
  -v, --verbose  Print the tree to the terminal.
```

To process the demo input file and generate a TSV file that can be opened by Excel

```
process.py <input TSV file> <output TSV file>
```


```
process.py demo-input.txt demo-output.xls
```

## Input file (TSV)

The tab delimited file must contain a header line and the following twelve columns:

NAME (REF DOI)	LABEL (REF DETAILS)	ID (REF NOTES)	NOTE	ABBR	TYPE	CHILDREN	CELLS	GENES	PROTEINS	PROTEOFORMS	LIPIDS	METABOLITES	FTUs	REFERENCES

The **Type** value needs to be "AS" for anatomical structures and "CT" for cell types. It doesn't matter what type values are used for the other items, so long as it's not either AS or CT.

**Children** is a comma separated list of child anatomical structure (AS) objects. These children need to be either anatomical structures (AS). The **Cells**, **Genes**, **Proteins**, **Proteoforms**, etc fields should be comma separated lists of the appropriate objects (e.g., **Cells**, should be a comma separated list of relevant cells). In all cases the objects **Name** or **Ref DOI** should be used. 

The first line in the input file is assumed to contain a header and is ignored.

The following example is incomplete and just included to exemplify the field values and usage:

```
NAME (REF DOI)	LABEL (REF DETAILS)	ID (REF NOTES)	NOTE	ABBR	TYPE	CHILDREN	CELLS	GENES	PROTEINS	PROTEOFORMS	LIPIDS	METABOLITES	FTU	REFERENCES (NAME/DOI)
ovary		UBERON:0000992	AS	central ovary, lateral ovary, medial ovary, mesovarium, ovarian ligament	hilum of ovary
central ovary			AS	central inferior ovary, central superior ovary	
lateral ovary			AS	lateral inferior ovary, lateral superior ovary	
medial ovary			AS	medial inferior ovary, medial superior ovary	
mesovarium		UBERON:0001342	AS		
ovarian ligament		UBERON:0008847	AS		
hilum of ovary			AS	ovarian artery, ovarian vein, pampiniform plexus, rete ovarii	hilar cell
corona radiata		CL:0000713	CT									doi:10.1093/oxfordjournals.humrep.a136365
hilar cell		CL:0002095	CT				alkaline phosphatase, acid phosphatase, non-specific esterase, inhibin, calretinin, melan-A, cholesterol esters					McKay et al 1961, Boss et al 1965, Mills et al 2020, Jungbluth et al 1998, Pelkey et al 1998
mural granulosa cell			CT									doi:10.1093/oxfordjournals.humrep.a136365
primary oocyte		CL:0000654	CT									doi:10.1093/oxfordjournals.humrep.a136365
secondary oocyte		CL:0000655	CT									doi:10.1093/oxfordjournals.humrep.a136365
columnar ovarian surface epithelial columnar cell			CT				calretinin, mesothelin					Mills et al 2020, Reeves et al 1971, Hummitzsch et al 2013, Blaustein et al 1979, McKay et al 1961
flattened cuboidal ovarian surface epithelial cell			CT				oviduct-specific glycoprotein-1, E-cadherin					Mills et al 2020, Reeves et al 1971, Hummitzsch et al 2013, Blaustein et al 1979, McKay et al 1961
oviduct-specific glycoprotein-1			Protein									
mesothelin			Protein									
E-cadherin			Protein									
doi:10.1093/oxfordjournals.humrep.a136365	PMID: 3558758		Reference									
McKay et al 1961	McKay, D., Pinkerton, J., Hertig, A. & Danziger, S. (1961). The Adult Human Ovary: A Histochemical Study. Obstetrics & Gynecology, 18(1), 13-39. 		Reference									
```

## Known problems and limitations

1. The program should validate the biomarkers using the TYPE field designation.

1. Export ASCT+B table as CSV file.

1. Need to allow for case-independence. At present if a cell type is defined with upper cases and applied to a structure in lower case then the program will consider these different entities and throw an error.

1. Need better example and docs.

1. The program should allow for non-unique "author preferred name" field values.

1. Test if a parent has a child which is actually the parent.

1. Need to strip blank space for the left/right of each text field.

1. If there are dulicate references in a comma separated list of references (column 15) this causes an error claiming a feature hasn't been defined when actually the error is about duplications in the comma separated list of features. The program needs to test for duplications in feature lists.

1. The program requires UTF-8 encoding.

1. There is a bug in Excel where by when generating TSV or CSV files, it may incorrectly include a lot of empty COLUMNS. For example, if the input file only has 15 columns of data, Excel may generate a TSV or CSV file that includes 30 columns that correctly includes the 15 columns of data and another 15 empty columns. This causes the program to error. The program needs to include a workaround for this issue. 

   ```
   ERROR: incorrect number of fields in line. The tab-delimited line should contain the following 15 fields: 
   	name, label, ID, node, abbreviation, feature type, children, cells, genes, proteins, proteoforms, lipids, metabolites, FTU, references
   	Number of fields found in line: 
   ```

1. There is a bug in Excel where by when generating TSV or CSV files, it may incorrectly include a lot of empty ROWS. For example, if the input file only has 100 rows of data, Excel may generate a TSV or CSV file that includes 500 rows that correctly includes the 100 rows of data and another 400 empty rows. This causes the program to error. The program needs to include a workaround for this issue. 

   ```
   ERROR: all features must have a unique name.
   	[ type:]
   	[ type:]
   ```

