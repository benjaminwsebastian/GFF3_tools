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
    gff_dict = {
        # if we want decoding of the URL-encoded strings we can use
        # "seq_id": None if seq_id == "." else urllib.unquote(seq_id),
        # but we do not need it
        "seq_id": None if seq_id == "." else seq_id,
        "source": None if source == "." else source,
        "seq_type": None if seq_type == "." else seq_type,
        "start": None if start == "." else int(start),
        "end": None if end == "." else int(end),
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
        "+" : [],
        "-" : []
    }
    
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
                        genes["+"].append(gene)
                    else:
                        genes["-"].append(gene)
                    
    return genes

## End reading GFF functions

## Begin classifying regions functions ##

def check_start(genes, sample):
    closest_gene = (None, 10000000)
        
    for gene in genes[sample["Sense"]]:
        diff = sample["Beg"] - gene.start
            
        if diff < closest_gene[1] and diff >= 0: # Check if it could be in a gene
            closest_gene = (gene, diff)
            
    return closest_gene

def check_end(genes, sample):
    closest_gene = (None, 10000000)
       
    for gene in genes[sample["Sense"]]:
        diff = sample["End"] - gene.start
        
        if diff < closest_gene[1] and diff >= 0: # check if it could be in a gene
            closest_gene = (gene, diff)
            
    return closest_gene

def where_in_gene(gene, sample):
    closest_region_beg = (None, 10000000)
    for region in gene.regions:
        diff = sample["Beg"] - region.start
            
        if diff <= closest_region_beg[1] and diff >= 0:
            closest_region_beg = (region, diff)
    
    closest_region_end = (None, 10000000)
    for region in gene.regions:
        diff = region.end - sample["End"]
            
        if diff <= closest_region_end[1] and diff >= 0:
            closest_region_end = (region, diff)

    if closest_region_beg[0] is not None and closest_region_end[0] is not None:
        if closest_region_beg[0].start == closest_region_end[0].start:
            return closest_region_beg[0].seq_type + ',' + gene.name
        else:
            return closest_region_beg[0].seq_type + ',' + closest_region_end[0].seq_type + ',' + gene.name
    elif closest_region_beg[0] is not None:
        return closest_region_beg[0].seq_type + ',' + gene.name
    else:
        return closest_region_end[0].seq_type + ',' + gene.name
        
def classify_regions(genes, locations, params):
    """ Classifies coordinates into regions """
    classified_regions = []
    for index, sample in locations.iterrows():
        if sample["Sense"] == "C":
            sample["Sense"] = "-"
        
        closest_gene_start = check_start(genes, sample)
        closest_gene_end = check_end(genes, sample)
        
        if closest_gene_start[0] is not None:
            if sample["Beg"] < closest_gene_start[0].end: # sample beg is contained in a gene
                if sample["End"] <= closest_gene_start[0].end: # sample is entirely contained within a gene
                    classified_regions.append(where_in_gene(closest_gene_start[0], sample))
                    
                elif closest_gene_end[0] is not None:
                    if sample["End"] >= closest_gene_end[0].start and closest_gene_start[0].name != closest_gene_end[0].name: # across multiple genes
                        classified_regions.append(where_in_gene(closest_gene_start[0], sample) + ";" + where_in_gene(closest_gene_end[0], sample)) # add areas in first and last gene (there may be genes inbetween which sample spans
                        
                    elif sample["Sense"] == "-": # in the promotor area and in gene
                        classified_regions.append("Promotor region," + where_in_gene(closest_gene_start[0], sample))
                        
                    else: # in the gene and the 3' area (not promotor)
                        classified_regions.append(where_in_gene(closest_gene_start[0], sample))
            
            elif closest_gene_end[0] is not None:
                if sample["End"] < closest_gene_end[0].end: # sample beg not contained in gene
                    if sample["Sense"] == "+":
                        classified_regions.append("Promotor region," + where_in_gene(closest_gene_end[0], sample)) # end contained in gene and beg in promotor region
                    else:
                        classified_regions.append(where_in_gene(closest_gene_end[0], sample)) # end contained in gene and beg not in promotor region
                        
                elif (abs(sample["Beg"] - closest_gene_start[0].start) < params["promotor"] or abs(sample["End"] - closest_gene_start[0].start) < params["promotor"]) and sample["Sense"] == "+": # beg and end to the left of a gene in the promotor region
                    classified_regions.append("Promotor region," + closest_gene_start[0].name)
                    
                elif (abs(sample["End"] - closest_gene_end[0].end) < params["promotor"] or abs(sample["Beg"] - closest_gene_end[0].end) < params["promotor"]) and sample["Sense"] == "-": # beg and end to the right of a gene and in the promotor region
                    classified_regions.append("Promotor region," + closest_gene_end[0].name)
                    
                else: # intergenic
                    classified_regions.append("intergenic")
                    
        else: # intergenic
             classified_regions.append("intergenic")
    return classified_regions
        
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
    
    genes = read_gff(gff_in, params)
    locations = pd.read_csv(csv_file_name, sep=params["sep"])
    
    # remove rows that don't begin with "###"
    locations = locations[locations['Score'].str.contains("###")]
    
    # rename columns so they work with python
    locations.columns = ['Score', '%_Div', '%_Del', '%_Ins', 'Query', 'Beg', 'End',
       'Length', 'Sense', 'Element', 'Family', 'Pos_Repeat_Beg',
       'Pos_Repeat_End', 'Pos_Repeat_Left', 'ID', 'Num_Assembled',
       '%_of_Ref']
    
    classified_regions = classify_regions(genes, locations, params)
    
    locations.insert(loc=1, column='gene_region', value=classified_regions)
    
    locations.to_csv(params["out_file"], sep='\t')
    
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
        "chrom" : str(chrom),
        "sep" : sep,
        "promotor" : args.p,
        "out_file" : f_out
    }
    
    main(gff_input, csv_file_name, params)
    gff_input.close()