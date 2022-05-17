#!/usr/bin/env python3
import argparse
import os
from pyfaidx import Fasta
import tarfile


def compress_dir(source_dir, name):
    """
    Compress directory into a tar.gz
    """
    with tarfile.open(name, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))

def parse_args():
    parser = argparse.ArgumentParser(description = 'Prepares the data needed by Armadillo by providing just the coordinates of the regions of interest.')
    parser.add_argument(
    '-g', '--reference', type = str, required = True,
    help = 'Reference genome')
    parser.add_argument(
    '-i', '--rois', type = str, required = True,
    help = 'Input file with regions of interest (BED or chr:st-end formatted)')
    parser.add_argument(
    '-I', '--identity', type = int, default = 90, required = False,
    help = 'Minimum identity for a hit to be considered a copy of the ROI (default: %(default)s).')
    parser.add_argument(
    '-L', '--lendiff', type = int, required = False, default = 15,
    help = 'Maximum length difference allowed between each hit and the input sequence (default: %(default)s%%).')
    parser.add_argument(
    '-m', '--mlen', type = int, default = 100, required = False,
    help = 'Minimum length (bp) allowed to each gene (default: %(default)s).')
    parser.add_argument(
    '-o', '--outputName', type = str, required = False, default = 'armadillo_data',
    help = 'Set name for the output directory (default: %(default)s).')
    return parser.parse_args()

def blat_parser(blat, filename):
    print_fa = False
    printed = 0
    for hit in blat:
        column=hit.strip().split()
        chr_stend=column[0].split(":")
        st_end=chr_stend[1].split("-")
        endd=st_end[1].split("_")
        querylength=int(endd[0])-int(st_end[0])
        if float(column[2])>args.identity and float(column[3])/querylength*100 >= 100-args.lendiff and float(column[7])/querylength*100 <=100+args.lendiff and float(column[5])<=3: #90% identity, alignment length +/- 15%, gaps 1
            print_fa = True
            printed += 1
            log=open("rois_copies_coords/"+filename, "a+") #We'll have a log file where we'll append any useful coordinate, so we don't use it twice
            if int(column[8])<int(column[9]): #We consider the start as the shortest coord, so if the it's in the negative strand, flip the coords
                log.write(column[1]+":"+column[8]+"-"+column[9]+"\n")
            else:
                log.write(column[1]+":"+column[9]+"-"+column[8]+"\n")
            log.close()
    if printed == 1:
        os.remove("rois_copies_coords/"+filename)
        print_fa = False
    return(print_fa)

args = parse_args()
href = Fasta(args.reference, rebuild=False) #Open the reference genome

##PREPARE THE DIRECTORIES##
input_file = os.path.abspath(args.rois)

os.mkdir(args.outputName)
d = os.getcwd()
os.chdir(d+"/"+args.outputName)
os.mkdir("rois_copies_coords")
os.mkdir("miniFASTA")

##READ THE INPUT ##
rois = open(input_file)
coord = None
hits = []
for line in rois:
    col = line.strip().split()
    chrom, start, end = col[0].replace(":", "-").split("-")

    if int(end) - int(start) < args.mlen: #Remove too short genes
        continue

    if col[0] == coord:
        hits.append(line.strip()) # Append lines if still the same query name
        continue

    if coord is None: # The first time it'll be blank, so we dont wont it to print anything
        coord = col[0] # Set the name of the new query
        continue
    
    coord = col[0] # Set the name of the new query
    extended_coord = f'{chrom}:{int(start)-100}-{int(end)+100}'
    
    if len(hits) > 1:
        print_fasta = blat_parser(hits, coord) # Parse the hits
        if print_fasta:
            miniFASTA = open("miniFASTA/" + coord + ".fa", "w+") #Create a FASTA for each region
            miniFASTA.write(">" + extended_coord + "\n" + href[chrom][int(start)-101:int(end)+100].seq+"\n") #We'll align the reads against a quite larger region so that the ones that overlap only in the flanks can align too
            miniFASTA.close()
    else:
        pass
    
    hits = [line.strip()] # Reset the list when changing the query name

## Remove duplicates
files = [file for file in os.listdir("rois_copies_coords")]
for file in files:
    for line in open("rois_copies_coords/"+file):
        line = line.strip()
        if line != file and line in files and file in files:
            os.remove("rois_copies_coords/"+line)
            os.remove("miniFASTA/"+line+".fa")
            files.remove(line)
            print("Removed ", line, ". Duplicate of: ", file, sep = "", end = "\n")
        else:
            continue

## Index the fasta files ##
fastas = [file for file in os.listdir("miniFASTA") if file.endswith(".fa")]
concat = ''.join([open("miniFASTA/"+f).read() for f in fastas])
refFASTA = open("armadillo_reference_genome.fa", "w+") #We'll write a merged reference genome
refFASTA.write(concat)
refFASTA.close()
os.system("samtools faidx armadillo_reference_genome.fa")

for fasta in fastas:
    if os.path.isfile("miniFASTA/"+fasta+".fai"): #Don't reindex if already done
        continue
    else:
        os.system("bwa index miniFASTA/" + fasta + ">/dev/null 2>&1")
        os.system("samtools faidx miniFASTA/" + fasta + ">/dev/null 2>&1")

## Create a list of all final rois
roisfile = open("rois", "w+")
for roi in [file for file in os.listdir("rois_copies_coords")]:
    roisfile.write(roi+"\n")
roisfile.close()

os.chdir(d)
compress_dir(os.path.abspath(args.outputName), f'/batchx/output/{args.outputName}.tar.gz')
