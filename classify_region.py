from collections import namedtuple
import gzip
import re
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pandas as pd

# region constants and variables
GFF_FIELDS = ["seq_id", "source", "seq_type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", GFF_FIELDS)
cur_gene = None
CHROM_LEGEND = {
    "GK000068.1" : 1,
    "GK000069.1" : 2,
    "GK000070.1" : 3,
    "GK000071.1" : 4,
    "GK000072.1" : 5,
    "GK000073.1" : 6,
    "GK000074.1" : 7,
    "GK000075.1" : 8
}
# For the intergenic classification case:
intergenic_tup = namedtuple("intergenic", ["start", "end", "name"])
default = {"start": None, "end": None, "name": None} 
intergenic = intergenic_tup(**default)

# endregion


#############
## CLASSES ##
#############

class GFFGene(object):
    """ Create data structure for storing genes. """

    def __init__(self, gff_record, transcripts=None, regions=None, misc_features=None):
        self.gff_record = gff_record
        self.start = gff_record.start
        self.end = gff_record.end
        self.length = self.end - (self.start - 1)
        self.name = parse_gff_attributes(self.gff_record.attributes)["Name"]
        self.id = parse_gff_attributes(self.gff_record.attributes)["ID"]
        self.transcripts = transcripts if isinstance(transcripts, list) else []
        self.regions = regions if isinstance(regions, GFFRecord) else []
        self.processed = False


class GFFTranscript(object):
    """ Create data structure for storing transcripts. """

    def __init__(self, gff_record, exons=None):
        self.gff_record = gff_record
        self.id = parse_gff_attributes(self.gff_record.attributes)["ID"]
        self.exons = exons or []
        # self.transcript = transcript


#################
## SUBROUTINES ##
#################

## Begin reading GFF functions ##

def parse_gff_attributes(attr_str):
    """
    Parse GFF3 attribute column and return a dict with all attributes.
    :param attr_str:
    :type attr_str:
    :return: dictionary of attributes
    :rtype: dict
    """
    retval = dict()     # return value
    if attr_str == ".":
        return {}
    for attr in attr_str.split(";"):
        if not attr: continue
        key, value = attr.split("=")
        # if we want decoding of the URL-encoded string we can use
        # retval[urllib.unquote(key)] = urllib.unquote(value)
        # but we do not need it
        retval[key] = value
    return retval

def create_gff_record(seq_id, source, seq_type, start, end, score, strand, phase, attributes):
    """ Create GFF record from the provided variables. """
    start = None if start == "." else int(start)
    end = None if end == "." else int(end)
    if end is not None and start is not None:
        if end < start:
            tmp = end
            end = start
            start = end
            
    if strand == "-" or strand == "C":
        if start is not None and end is not None:
            tmp = end
            end = -1*start
            start = -1*tmp
    
    gff_dict = {
        # if we want decoding of the URL-encoded strings we can use
        # "seq_id": None if seq_id == "." else urllib.unquote(seq_id),
        # but we do not need it
        "seq_id": None if seq_id == "." else seq_id,
        "source": None if source == "." else source,
        "seq_type": None if seq_type == "." else seq_type,
        "start": start,
        "end": end,
        "score": None if score == "." else float(score),
        "strand": None if strand == "." else strand,
        "phase": None if phase == "." else phase,
        "attributes": attributes
    }
    gff_rec = GFFRecord(**gff_dict)
    return gff_rec

def reset_gene():
    """ Reset global variables that are used to keep track of the current gene data. """
    global cur_gene
    cur_gene = None

def create_gene_gff_record_from_line(line):
    """ Create and return GFF record from the line (string) only if line contains records that we need. """
    fields = line.strip().split("\t")
    if len(fields) != len(GFF_FIELDS):
        return None
    if fields[2] not in ("gene", "mRNA", "CDS", "exon", "intron", "five_prime_utr", "three_prime_utr", "C_gene_segment", "V_gene_segment", "D_gene_segment",
                         "J_gene_segment"):
        return None
    gff_record = create_gff_record(*fields)
    return gff_record


def analyze_gff_record(gff_record):
    """ Analyze GFF record. If this is a new gene then analyze and process previous gene. """
    global cur_gene
    # check type of the record and act accordingly
    if gff_record.seq_type == "gene":
        # this is gene: process old gene or create new empty (first case) gene
        if cur_gene is not None:
            gene = cur_gene
            cur_gene.processed = True
                
            reset_gene()
            # create new gene based on gff record
            cur_gene = GFFGene(gff_record=gff_record)
            
            return gene
        else:
            reset_gene()
            cur_gene = GFFGene(gff_record=gff_record)
            
            return None
        
    elif gff_record.seq_type in ("exon", "intron", "CDS", "five_prime_utr", "three_prime_utr"):
        # this is a category we want: add to the current gene's regions
        cur_gene.regions.append(gff_record)
    elif gff_record.seq_type in ("mRNA", "C_gene_segment", "V_gene_segment", "D_gene_segment", "J_gene_segment"):
        # this is a category we don't want: pass
        pass
    else:
        raise Exception("ERROR: Something went wrong while processing GFF record:" + str(gff_record))


def read_gff(gff_in, params):
    """
    Reads in genes from gff file
    :param gff_in: input file handle
    :type gff_in: file
    """
    genes = {
        "+" : {},
        "-" : {}
    }
    
    starts = {"+" : [] , "-" : []}
    gff_record = None
    # parse lines
    for line in gff_in:
        if line.startswith("#"):
            pass # pass the preamble
        elif line.strip() != '':
            # this is data line, process it
            gff_record = create_gene_gff_record_from_line(line)     # treat only if record type is what we want
            if isinstance(gff_record, GFFRecord):
                gene = analyze_gff_record(gff_record)
                # add to appropriate list and filter on chromosome
                if gene is not None and gene.gff_record.seq_id == params["chrom"]:
                    if gene.gff_record.strand == "+":
                        genes["+"][gene.start] = gene
                        starts["+"].append(gene.start)
                    else:
                        genes["-"][gene.start] = gene
                        starts["-"].append(gene.start)
    starts["+"].sort()
    starts["-"].sort()
    
    return genes, starts

## End reading GFF functions

## Begin classifying regions functions ##

def get_coverage(region, TE):
    len_te = abs(TE["End"] - TE["Beg"])
    post, pre = 0, 0
    if TE["End"] > region.end:
        post = abs(TE["End"] - region.end)
    if TE["Beg"] < region.start:
        pre = abs(region.start - TE["Beg"])
    coverage = abs(len_te - post - pre)/abs(region.end - region.start)
    return coverage

def add_region(classified_regions, gene, TE, region_type, region_coverage, **kwargs):
    region_start, region_end, gene_start, gene_end = None, None, None, None
    try:
        region_start = min(abs(kwargs["region_start"]), abs(kwargs["region_end"]))
        region_end = max(abs(kwargs["region_start"]), abs(kwargs["region_end"]))
    except:
        pass
    
    try:
        gene_start = min(abs(gene.start), abs(gene.end))
        gene_end = max(abs(gene.start), abs(gene.end))
    except:
        pass
    
    TE_beg = min(abs(TE["Beg"]), abs(TE["End"]))
    TE_end = max(abs(TE["Beg"]), abs(TE["End"]))
    
    classified_regions = classified_regions.append({
        "Gene ID": gene.name,
        "Sense": TE["Sense"],
        "Gene start": gene_start if gene.start is not None else None,
        "Gene end": gene_end if gene.end is not None else None,
        "Region type": region_type,
        "Region start": region_start if region_start is not None else None,
        "Region end": region_end if region_end is not None else None,
        "TE ID": TE["ID"],
        "TE family": TE["Family"],
        "TE start": TE_beg,
        "TE stop": TE_end,
        "Region coverage": round(abs(region_coverage), 4)
    }, ignore_index = True)
               
    return classified_regions
               
def add_promoter(classified_regions, right_start, genes, TE, params):
    # Calculate percent of promoter region covered by the TE
    # The promoter region coverage is the length of the TE within the given promoter region length (default 1000) dividied by the given promoter region length
    start_of_promoter = genes[TE["Sense"]][right_start].start - params["promoter"]
    post, pre = 0, 0
    if TE["End"] > genes[TE["Sense"]][right_start].start:
        post =  abs(TE["End"] - genes[TE["Sense"]][right_start].start)
    if TE["Beg"] < start_of_promoter:
        pre = abs(start_of_promoter - TE["Beg"])
    len_in_prom = abs(TE["End"] - TE["Beg"]) - post - pre
    coverage = len_in_prom / params["promoter"]
    
    # Is the TE in the area before the start of the promoter region?
    pre_promoter = False
    if TE["Beg"] < start_of_promoter:
        pre_promoter = True
        
    # Add the intergenic region before the start of the promoter region that the TE covers
    if pre_promoter:
        classified_regions = add_region(classified_regions, intergenic, TE, "Intergenic", 1, region_start = TE["Beg"], region_end = abs(start_of_promoter-1))
    classified_regions = add_region(classified_regions, genes[TE["Sense"]][right_start], TE, "Promoter", coverage, region_start = start_of_promoter, region_end = genes[TE["Sense"]][right_start].start)
    
    return classified_regions
    
def where_in_genes(left_start, right_start, starts, genes, TE, classified_regions, params):
    # Define starting point - does the current point begin in a gene or does TE end after the beginning of the next gene (exclusive "or" used)
    covered_genes = []
    already_classified = False
    if left_start is not None and left_start + genes[TE["Sense"]][left_start].length > TE["Beg"]:
        # LEFT OVERLAP start
        already_classified = True
        covered_genes.append(left_start)
        
    if right_start is not None and TE["End"] > genes[TE["Sense"]][right_start].start:
        # RIGHT OVERLAP start
        
        if not already_classified:
            # It must be the case that the gene starts in the promoter region - or crosses through it - before entering a downstream (with respect to sense) gene
            classified_regions = add_promoter(classified_regions, right_start, genes, TE, params)
        
        already_classified = True
        covered_genes.append(right_start)
        
    if not already_classified:
        # This case, the gene must be intergenic or just in a promoter region
        return classified_regions, already_classified
    
    all_genes_found = False
    while not all_genes_found:
        # Get the location of the next nearest upstream gene
        try:
            right_start = min([x for x in starts[TE["Sense"]] if x > right_start])
        except:
            pass
        
        if right_start is not None and TE["End"] > genes[TE["Sense"]][right_start].start:
            covered_genes.append(right_start)
        else:
            all_genes_found = True
    
    for gene_num, gene in enumerate(covered_genes):
        gene = genes[TE["Sense"]][gene]
        regions = {}
        region_starts = []
            
        for region in gene.regions:
            if region.seq_type.upper() not in ("CDS",):
                regions[region.start] = region
                region_starts.append(region.start)
                
        region_starts.sort()
        covered_regions = []
        try:
            covered_regions = [x for x in region_starts if x > TE["Beg"] and x < TE["End"]]
        except:
            pass
        
        # First add the gene and the TE, with the percent of gene covered by the TE
        coverage = get_coverage(gene, TE)
        classified_regions = add_region(classified_regions, gene, TE, "Gene", coverage)
        
        # Next add each component
        
        # Get the starting region
        try:
            left = max([x for x in region_starts if x < TE["Beg"]])
        except:
            left = min(region_starts) # The case where the TE begins outside of the gene
        try:
            i = region_starts.index(min([x for x in region_starts if x > TE["End"]])) - 1
            right = region_starts[i]
        except:
            right = None
         
        # Add the starting region
        coverage = get_coverage(regions[left], TE)
        classified_regions = add_region(classified_regions, gene, TE, regions[left].seq_type, coverage, region_start = regions[left].start, region_end = regions[left].end)
        
        # Add the middle regions
        try:
            covered_regions.remove(left)
        except:
            pass
        try:
            covered_regions.remove(right)
        except:
            pass
        
        for region in covered_regions:
            coverage = get_coverage(regions[region], TE)
            classified_regions = add_region(classified_regions, gene, TE, regions[region].seq_type, coverage, region_start = regions[region].start, region_end = regions[region].end)

            # Add the ending region
        # If there are no covered regions, then a region is entirely covering the TE or the TE is between 2 regions
        if left != right and right is not None:
            coverage = get_coverage(regions[right], TE)
            classified_regions = add_region(classified_regions, gene, TE, regions[right].seq_type, coverage, region_start = regions[right].start, region_end = regions[right].end)

        # Ending of TE is outside of the gene
        elif right is None:
            # Not in promoter region of next gene
            if gene_num == len(covered_genes) - 1:
                classified_regions = add_region(classified_regions, intergenic, TE, "Intergenic", 1)
                
            # In promoter region of next gene
            elif abs(genes[TE["Sense"]][covered_genes[gene_num + 1]].start - TE["End"]) < params["promoter"]:
                classified_regions = add_promoter(classified_regions, genes[TE["Sense"]][covered_genes[gene_num + 1]].start, genes, TE, params)        
        
    return classified_regions, already_classified
    
def classify_regions(genes, starts, TEs, params):
    """ Classifies coordinates into regions """
    
    # Columns for output file
    columns = [
        "Gene ID",
        "Sense",
        "Gene start",
        "Gene end",
        "Region type",
        "Region start",
        "Region end",
        "TE ID",
        "TE family",
        "TE start",
        "TE stop",
        "Region coverage"
    ]
    classified_regions = pd.DataFrame(columns = columns)
    for index, TE in TEs.iterrows():
        if TE["Sense"] == "C":
            TE["Sense"] = "-"
        if TE["Sense"] == "-":
            tmp = -1*TE["Beg"]
            TE["Beg"] = -1*TE["End"]
            TE["End"] = tmp
        
        right_start = None
        left_start = None
        
        # Get the location of the nearest upstream and downstream genes to the start of the TE
        try:
            right_start = min([x for x in starts[TE["Sense"]] if x > TE["Beg"]])
        except:
            pass
        try:
            left_start = max([x for x in starts[TE["Sense"]] if x < TE["Beg"]])
        except:
            pass

        classified_regions, already_classified = where_in_genes(left_start, right_start, starts, genes, TE, classified_regions, params) 
        
        if not already_classified and right_start is not None and genes[TE["Sense"]][right_start].start - TE["End"] < params["promoter"]:
            # PROMOTER REGION - The TE is just in the promoter region
            already_classified = True
            
            # Need to classify intergenic region before promoter
            classified_regions = add_promoter(classified_regions, right_start, genes, TE, params)
            
        if not already_classified:
            # INTERGENIC
            
            classified_regions = add_region(classified_regions, intergenic, TE, "Intergenic", 1)
    return classified_regions
        
# End classifying regions section

# Begin post processing section

def mark_regions(classified_regions):
    for region in ("exon", "intron", "five_prime_utr", "three_prime_utr"):
        
        current_region = classified_regions[classified_regions["Region type"] == region].copy()
        grouped_regions = current_region.groupby(["Region start"])
        
        for i, group in enumerate(grouped_regions):
            for index in group[1].index:
                classified_regions.iloc[index]["Region type"] = region + str(i)
            
    return classified_regions

# End post processing section

###########
## INPUT ##
###########
          
def get_args():
    parser = ArgumentParser(description="""
    Application to add gene regions to csv files from based
    on gene locations in the csv file
    The input file should be a GFF (version 3) or gziped GFF (version 3).
    and csv file with locations as given in fiverr chat

    Usage:
    python3 classify_region.py --i input.gff[.gz, .gff3]
    """, formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument("--csv", type = str, required = True, help = "csv file path (just put it's name if it is in the same directory)")
    parser.add_argument("--gff", type = str, required = True, help = "gff file (or path)")
    parser.add_argument("--o", type = str, required = False, default = None, help = "output name")
    parser.add_argument("--chrom", type = int, required = False, default = None, help = "chromosome number")
    parser.add_argument("--delim", type = str, required = False, default = "tab", help = "delimiter of the csv file. Options: tab, comma")
    parser.add_argument("--p", type = int, required = False, default = 1000, help = "length of promotor region")

    args = parser.parse_args()
    
    return args
          
##########
## MAIN ##
##########

def main(gff_in, csv_file_name, params):
    
    print("Reading in gff file...")
    genes, starts = read_gff(gff_in, params)
    print("Reading in TE file...")
    TEs = pd.read_csv(csv_file_name, sep=params["sep"])
    
    # remove rows that don't begin with "###"
    TEs = TEs[TEs['Score'].str.contains("###")]
    
    # rename columns so they work with python
    TEs.columns = ['Score', '%_Div', '%_Del', '%_Ins', 'Query', 'Beg', 'End',
       'Length', 'Sense', 'Element', 'Family', 'Pos_Repeat_Beg',
       'Pos_Repeat_End', 'Pos_Repeat_Left', 'ID', 'Num_Assembled',
       '%_of_Ref']
    
    print("Classifying regions...")
    classified_regions = classify_regions(genes, starts, TEs, params)
    classified_regions = classified_regions.drop_duplicates()
    classified_regions = classified_regions.reset_index(drop = True)
    
    classified_regions = mark_regions(classified_regions)
    
    print("Saving to file...", params["out_file"])
    classified_regions.to_csv(params["out_file"], sep='\t')
    
if __name__ == "__main__":
    args = get_args()
    
    gff_file_name = args.gff   # input file name
    file_open = gzip.open if gff_file_name.endswith(".gz") else open
    gff_input = file_open(gff_file_name, 'rt', encoding = 'utf-8')   # input file handle
    
    csv_file_name = args.csv
    
    # infor chrom from file name if not specified
    chrom_str = csv_file_name.split('_')[1].split('.')[0] + ".1"
    chrom = CHROM_LEGEND[chrom_str]
    if args.chrom is not None:
        chrom = args.chrom
        
    f_out = "chrom_" +  chrom_str + "_labeled_gene_regions.csv"
    if args.o is not None:
        f_out = args.o
        
    sep = '\t'
    if args.delim == "comma":
        sep = ',' 
    
    params = {
        "chrom": str(chrom),
        "sep": sep,
        "promoter": args.p,
        "out_file": f_out
    }
    
    main(gff_input, csv_file_name, params)
    gff_input.close()
