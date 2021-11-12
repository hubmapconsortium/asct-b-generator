# ASCT+B Generator

This program converts a simple TSV file into a HuBMAP ASCT+B table.

# Input file (TSV)

The tab delimited file should contain six columns: name, label, id, type, and children

**Name**, **label**, and **ID** are the HuBMAP values. 

**Type** is one of the following:
AS - anatomical structure
CT - cell type
BG - gene

**Children** is a common separated list of child objects (AS, CT, and/or BG)

A header line is optional and will be ignored.

**Example**:

``NAME	LABEL	ID	TYPE	CHILDREN`
`ovary		UBERON:0000992	AS	"central ovary, lateral ovary, medial ovary, mesovarium, ovarian ligament, hilum of ovary"`
`central ovary			AS	"central inferior ovary, central superior ovary"`
`lateral ovary			AS	"lateral inferior ovary, lateral superior ovary"`
`medial ovary			AS	"medial inferior ovary, medial superior ovary"`
`columnar ovarian surface epithelial columnar cell			CT`	
`flattened cuboidal ovarian surface epithelial cell			CT`	
corpus luteum endothelial cell			BG`	
