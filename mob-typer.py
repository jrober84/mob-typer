#!/usr/bin/env python
#
# Plasmid in silico typing based on MOB categories and replicon typing
#
# Authors - James Robertson
#
# Dependencies:
#
#
from __future__ import division
from argparse import (ArgumentParser, FileType)
import os, sys, re, collections, operator, math
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import SeqIO
from Bio.SeqUtils import GC
import pprint
from subprocess import Popen, PIPE
import bestBlastHit


def run_makeblastdb(fasta_file):

    p = Popen(['makeblastdb',
               '-in', '{}'.format(fasta_file),
               '-dbtype', 'nucl'],
              stdout=PIPE,
              stderr=PIPE)
    p.wait()




def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(description='Classify a plasmid sequence in fasta format')
    parser.add_argument('--outdir', type=str, required=True, help='Directory to put results')
    parser.add_argument('--infile', type=str, required=True, help='Input file to process')
    parser.add_argument('--evalue', type=str, required=False, help='evalue threshold',default=10)
    parser.add_argument('--repdb', type=str, required=False, help='Replicon protien typing database',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),'databases/rep.proteins.faa'))
    parser.add_argument('--pfinder', type=str, required=False, help='PLasmid finder DNA typing database',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),'databases/rep.elements.fas'))
    parser.add_argument('--mobdb', type=str, required=False, help='Relaxase typing database',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'databases/mob.proteins.faa'))
    parser.add_argument('--mpfdb', type=str, required=False, help='Mating pair formation typing database',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'databases/mpf.proteins.faa'))
    parser.add_argument('--repdb_type', type=str, required=False, help='Replicon typing database blast method',default='blastp')
    parser.add_argument('--pfinderdb_type', type=str, required=False, help='Replicon typing database blast method',
                        default='blastn')
    parser.add_argument('--mobdb_type', type=str, required=False, help='Relaxase typing database blast method', default='blastp')
    parser.add_argument('--mpfdb_type', type=str, required=False, help='Mating pair formation database blast method',
                        default='blastp')
    return parser.parse_args()

def run_blast(query_fasta,blast_db,method,evalue,outfile,wordsize=4):
    comline = NcbiblastxCommandline(query=query_fasta, db=blast_db, cmd=method, evalue=evalue, outfmt=5,
                                    out=outfile, word_size=wordsize)
    os.system(str(comline))

def run_prokka(outdir,locustag,prefix,fasta_file):
    comline = "prokka --compliant --force -o {} --compliant --addgenes --locustag {} --prefix {}  {} ;".format(outdir,locustag,prefix,fasta_file)
    os.system(str(comline))

def report(rephit,mobhit,mpfhit,prokka_stats):
    return

def calcFastaStats(fasta):
    num_seqs = 0;
    seq = ''
    for record in SeqIO.parse(fasta, "fasta"):
        num_seqs += 1
        seq = seq + record.seq
    genome_size = len(seq)
    gc = GC(seq)

    return {
        'num_seq': num_seqs,
        'size': genome_size,
        'gc_content': gc
    }

def getGenomeStats(nt_gene_fasta,full_fasta_file):
    genome_stat = calcFastaStats(full_fasta_file)
    gene_stat = calcFastaStats(nt_gene_fasta)

    return {
        'noncoding_bases': genome_stat['size'] - gene_stat['size'],
        'coding_percentage': gene_stat['size']/genome_stat['size']*100,
        'total_bases': genome_stat['size'],
        'num_genes': gene_stat['num_seq'],
        'gc_content': genome_stat['gc_content']
    }



def getOverlaps(positions):
    pStart = -1
    pEnd = -1
    regions = []
    size = len(positions)
    count = 0
    for start,end in positions:
        if pStart == -1:
            pStart = start
            pEnd = end
        elif(start > pEnd):
            regions.append([pStart, pEnd])
            pStart = start
            pEnd = end
        elif(end > pEnd):
            pEnd = end
        elif(size == count):
            if (start > pEnd):
                regions.append([pStart, pEnd])
            elif (end > pEnd):
                regions.append([pStart, end])
        count+=1
    return regions

def findNonCodingRegions(genbank_file):
    positions = list()
    for rec in SeqIO.parse(genbank_file, "genbank"):
        for feature in rec.features:
            if feature.type == 'Source':
                continue
            start = feature.location.start.position
            end = feature.location.end.position
            if(start > end):
                temp = end
                start = end
                end = temp
            positions.append([start,end])
    positions.sort(key=lambda tup: tup[1])
    regions = getOverlaps(positions)
    pStart = 0
    pEnd =0
    noncoding_regions = []
    for start,end in regions:
        if(start == 0):
            noncoding_regions.append([pStart,start ])
            pEnd = end
        else:
            noncoding_regions.append([pEnd , start ])
            pEnd = end
    return noncoding_regions



def extract_seqs(fasta,regions):
    sequences ={}
    for record in SeqIO.parse(fasta, "fasta"):
        id = record.id
        seq = record.seq
        for start,end in regions:
            sequences[str(id)+'|'+str(start)+':'+str(end)] = seq[start:end]
            print str(id)+'|'+str(start)+':'+str(end) + "\t" + seq[start:end]
    return sequences


def write_results(infile,out_file,rephit,mobhit,mpfhit,genome_stats):
    header = [
        'fasta_file',
        'genome_size',
        'gc_content',
        'number_of_genes',
        'genome_coding_fraction',
        'replication_type',
        'mob_type',
        'mpf_type',
        'predicted_transferability',
    ]
    report = "\t".join(header)
    score = 0
    if len(rephit) > 0:
        #pprint.pprint(rephit)
        reptype = ','.join(rephit)
        #pprint.pprint(reptype)
    else:
        reptype = '-'
    if len(mobhit) > 0:
        mobtype = mobhit['top_hit_id'].split('|')
    else:
        mobtype = ['-']
    if len(mpfhit) > 0:
        mpftype = mpfhit['top_hit_id'].split('|')
    else:
        mpftype = ['-']

    if mobtype[len(mobtype)-1] != '-' and mpftype[len(mpftype)-1] != '-':
        predicted_transferability = 'conjugative'
    elif mobtype[len(mobtype)-1] != '-' :
        predicted_transferability = 'mobilizable'
    else:
        predicted_transferability = 'non-mobilizable'

    report = report + "\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(infile,genome_stats['total_bases'],genome_stats['gc_content'],genome_stats['num_genes'], \
                                                      genome_stats['coding_percentage'],reptype,mobtype[len(mobtype)-1],mpftype[len(mpftype)-1],predicted_transferability)
    target = open(out_file, 'w')
    target.write(report)
    target.close()
    print report


'''
    Accepts a blast hsps object and filtering criteria: identity, coverage, evalue and
    returns an hsp object if it passes or None
'''

def filterBlastHit(hspObj,minIdent,minCoverage,minEvalue,query_len ) :
    if hspObj == None or hspObj.expect > minEvalue \
            or (hspObj.positives/hspObj.align_length*100) < minIdent \
            or (hspObj.align_length/query_len*100) < minCoverage:
        hspObj = None
    return hspObj

'''
    Accepts an array of blast hsps object and filtering criteria: identity, coverage, evalue and
    returns a list of replicon families which were found to have valid blast hits
'''

def repBlastFilter(hsps,minIdent,minCoverage,minEvalue ) :
    valid = []
    queries = {}
    for hit_id in hsps['all_hits']:

        for query in hsps['all_hits'][hit_id]:

            obj = filterBlastHit(hsps['all_hits'][hit_id][query]['hsp'], minIdent, minCoverage, minEvalue, hsps['all_hits'][hit_id][query]['query_len'])
            if obj != None:
                if query in queries:
                    continue
                id = hit_id.split('|')
                valid.append(id[len(id)-1])
                queries[query ] = ''
    return valid




def main():
    args = parse_args()

    if not args.outdir:
        print "Error, no ouput directory specified, please specify one"
        sys.exit()
    if not args.infile:
        print "Error, no fasta specified, please specify one"
        sys.exit()
    print "Processing .... "+args.infile
    sample_id = os.path.splitext(os.path.basename(args.infile))[0]
    prokka_dir = os.path.join(args.outdir,'PROKKA_'+sample_id)
    #run_prokka(prokka_dir, sample_id, sample_id, args.infile)
    nt_fasta_file = os.path.join(prokka_dir, sample_id+".fna")
    nt_gene_fasta =  os.path.join(prokka_dir, sample_id+".ffn")
    aa_gene_fasta = os.path.join(prokka_dir, sample_id + ".faa")
    genome_stats = getGenomeStats(nt_gene_fasta,args.infile)
    #run_makeblastdb(nt_fasta_file)

    pfinderblastfile = os.path.join(prokka_dir, 'replicon_blast_results_dna.xml')
    #run_blast(aa_gene_fasta, args.repdb, args.repdb_type, args.evalue, repblastfile)
    #Run Plasmid finder database
    #run_blast(args.pfinder, nt_fasta_file, args.pfinderdb_type, args.evalue, pfinderblastfile, 8)
    if (os.path.getsize(pfinderblastfile) > 0):
        reptypesDNA = repBlastFilter(bestBlastHit.bestBlastHit(pfinderblastfile).getBestHit(), 80, 60, 0.00001)
    else:
        reptypesDNA = {}

    repblastfile = os.path.join(prokka_dir, 'replicon_blast_results_aa.xml')
    #run_blast(aa_gene_fasta, args.repdb, args.repdb_type, args.evalue, repblastfile)
    if (os.path.getsize(repblastfile) > 0):
        reptypesAA = repBlastFilter(bestBlastHit.bestBlastHit( repblastfile).getBestHit(), 80, 80, 0.00001)
    else:
        reptypesAA = {}

    reptypes = list(set(sorted(reptypesDNA  + reptypesAA)))

    mobblastfile =  os.path.join(prokka_dir, 'mob_blast_results.xml')
    #run_blast(aa_gene_fasta, args.mobdb, args.mobdb_type, args.evalue, mobblastfile)
    if (os.path.getsize(mobblastfile) > 0):
        mobhit = bestBlastHit.bestBlastHit( mobblastfile).getBestHit()
        mobhit['top_hsp'] = filterBlastHit(mobhit['top_hsp'], 80, 80, 0.00001, mobhit['query_len'])
    else:
        mobhit = {}

    mpfblastfile =  os.path.join(prokka_dir, 'mpf_blast_results.xml')
    #run_blast(aa_gene_fasta, args.mpfdb, args.mpfdb_type, args.evalue, mpfblastfile)
    if(os.path.getsize(mpfblastfile) > 0):
        mpfhit = bestBlastHit.bestBlastHit( mpfblastfile).getBestHit()
        mpfhit['top_hsp'] = filterBlastHit(mpfhit['top_hsp'], 80, 80, 0.00001, mpfhit['query_len'])
    else:
        mpfhit = {}

    if (len(mobhit) > 0 ):
        query_id = mobhit['query_id']
        if (len(mpfhit) > 0 and mpfhit['query_id'] == query_id):
            mpfhit = {}

    if not 'top_hsp' in mobhit or mobhit['top_hsp'] == None:
        mobhit = {}
    if not 'top_hsp' in mpfhit or mpfhit['top_hsp'] == None:
        mpfhit = {}


    write_results(args.infile, os.path.join(args.outdir, sample_id+'_report.txt'), reptypes, mobhit, mpfhit, getGenomeStats(nt_gene_fasta,args.infile))






# call main function
if __name__ == '__main__':
	main()