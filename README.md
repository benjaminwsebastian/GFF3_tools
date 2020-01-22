Installation/Running instructions.

install reqs with:

pip3 -r install requirements.txt

There are 2 scripts, the first adds utr and introns to the gff and the seconds classifies the regions in the csv file.

## 1st. ##

python3 --i file.gff

If you don't want the script to overwrite the original file you can specify an output name with the "--o" tag

## 2nd ##

python3 classify_region.py --csv csv_file.csv --gff gff_file.gff

Example:

python3 classify_region.py --csv WGS.fasta.out_GK000068.1.elem_sorted.csv --gff Nelumbo_nucifera_sacred_lotus.gff

I also added a promotor region class with the promotor defined as 1000bp upstream. If you want to change this you can tune it with the "--prom" tag (--prom 100 for 100 bp upstream, --prom 120 for 120 bp upstream, and so on). You can check the syntax for this by running "python3 classify_region.py -h" or "python3 classify_region.py --help."

There are also some other options - you can check them by running:

python3 classify_region.py -h

