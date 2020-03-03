# Installation/Running instructions.

instantiate virtual environment, activate, then install requirements:
```
virtualenv .venv
source .venv/bin/activate
pip3 install -r requirements.txt
```

This repo contains 2 scripts, the first adds gene regions to a gff file (five_prime_utr, three_prime_utr, intron, exon) and the seconds classifies regions given in csv with the gff file that has the added regions. For example, a csv file containing transposons will be mapped to the gff file, with the output being the gene it crosses, the regions within the gene it crosses, and the coverage in all of these regions:

```
14		-			Intergenic			717580/717581	LTR/Copia	19878	21495	1
15		-			Intergenic			717582/717582	LTR/Copia	21496	22349	1
16	LOC104590429	+	23349	57568	Promoter	22349	23349	717583	LTR/Copia	23238	23886	0.111
17	LOC104590429	+	23349	57568	Gene			717583	LTR/Copia	23238	23886	0.0157
18	LOC104590429	+	23349	57568	exon0	23349	23786	717583	LTR/Copia	23238	23886	1
19	LOC104590429	+	23349	57568	intron0	23787	24656	717583	LTR/Copia	23238	23886	0.1139
20		-			Intergenic			717584	LTR/Copia	24106	24400	1
21		-			Intergenic			717585	DNAnona/MULE	24819	24984	1
22	LOC104590429	+	23349	57568	Gene			717586	SINE/unknown	25963	26088	0.0037
23	LOC104590429	+	23349	57568	intron1	25758	28403	717586	SINE/unknown	25963	26088	0.0473
24	LOC104590429	+	23349	57568	Gene			717587	DNAnona/MULE	27356	27474	0.0034
25	LOC104590429	+	23349	57568	intron1	25758	28403	717587	DNAnona/MULE	27356	27474	0.0446
26		-			Intergenic			717588	DNAnona/MULE	27748	27871	1
27		-			Intergenic			717589	LINE/L1	29339	29500	1
28		-			Intergenic			717590/717590/717590/717590	LTR/Copia	29520	34102	1

```

## Running add_gene_regions.py

python3 add_gene_regions.py --i file.gff

There are also some other options - you can check them by running:

`python3 add_gene_region.py -h` 

## Running classify_regions.py

`python3 classify_region.py --csv csv_file.csv --gff gff_file.gff`

Example:

python3 classify_region.py --csv transposons.csv --gff chomosome1.gff

There are also some other options - you can check them by running:

`python3 classify_region.py -h`

