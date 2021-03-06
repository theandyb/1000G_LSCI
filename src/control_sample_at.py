from pyfaidx import Fasta
import regex as re
import random
import pandas as pd
import argparse
from Bio.Seq import Seq

def full_motif(motif, ref):
    motif_seq = Seq(motif)
    motif_rc = motif_seq.reverse_complement().__str__()
    if ref in ['A','C']:
        final = "{}({})".format(motif, motif_rc)
    else:
        final = "{}({})".format(motif_rc, motif)
    return final

def main(ref_prefix = "chr"):
    parser = argparse.ArgumentParser(description="Sample control distribution.")
    parser.add_argument("-s", "--singleton", help="Location of singleton file", required=True)
    parser.add_argument("-f", "--fasta", help="FASTA file to grab sequence from", required=True)
    parser.add_argument("-o", "--output", help="Path to output", required=True)
    parser.add_argument("-n", "--nSample", help="Number of controls per singleton", type=int, default = 1)
    parser.add_argument("chrom", help="Chromosome we are sampling from")
    args = parser.parse_args()

    ref_file = args.fasta 
    singleton_file = args.singleton
    chrom = args.chrom
    nSample = args.nSample
    output_list = []
    # Create fasta object
    print("Reading fasta...")
    fasta_obj = Fasta(ref_file)
    seq = fasta_obj["{}{}".format(ref_prefix, chrom)] # hg37 does not have chr prefix, but hg38 does
    seqstr = seq[0:len(seq)].seq
    print("FASTA read!")
    # Iterate over singletons file
    print("Sampling control observations for singletons...")
    counter = 1
    with open(singleton_file) as fp:
        line = fp.readline()
        while line:
            # columns are: 0:chrom, 1:pos, 2:motif, 3:subtype, 4:alt, 5:id, 6:REF, 7:motif2, 8:subtype2
            content = line.strip().split(",")
            pos = int(content[1])
            ref = content[6]
            motif = content[7][:21]
            cat = content[8]
            if(not cat.startswith("A")):
                line = fp.readline()
                continue
            output_list.extend(sample_control(chrom, pos, ref, cat, seqstr, nSample))
            line = fp.readline()
            counter += 1
            if counter % 10000 == 0:
                print(counter)
                df = pd.DataFrame(output_list)
                df.to_csv(args.output, index = None, header=False, mode='a')
                output_list = []
    print("Done sampling...")
    if output_list:
        pd.DataFrame(output_list).to_csv(args.output, index = None, header=False, mode='a')
    print("Done!")

def sample_control(chrom, pos, ref, cat, seq, nSample, window=150, bp=10):
    sites1 = []
    sites2 = []
    newlist = []
    while(len(sites1) + len(sites2) < nSample + 5):
        lseg_lb = max((pos-1-window-bp), 0)
        lseg_ub = pos - bp - 1
        useg_lb = pos + bp
        useg_ub = upBound = min(len(seq), pos + window + bp)
        subseq1 = seq[lseg_lb:lseg_ub]
        subseq2 = seq[useg_lb:useg_ub]
        subseq1 = re.sub(r'^N+', '', subseq1)
        subseq2 = re.sub(r'N+$', '', subseq2)
        sites1 = [m.start() for m in re.finditer(ref, subseq1)]
        sites2 = [m.start() for m in re.finditer(ref, subseq2)]
        sites1 = [s for s in sites1 if (s >= bp+1 and s < (len(subseq1)-bp - 1))]
        sites2 = [s for s in sites2 if (s >= bp+1 and s < (len(subseq2)-bp -1))]
        window += 50 #expand window in edge case where mut_site is only ref_allele in window
    window -= 50
    while len(newlist) < nSample:
        flip = random.randint(0, 1)
        if ((flip == 0 and len(sites1)>0) or (len(sites2)==0)):
            subseq = subseq1
            sites = sites1
            c_direction = -1
        else:
            subseq = subseq2
            sites = sites2
            c_direction = 1
        if(len(sites)==0):
            print("Bad pos: {}".format(pos))
        ix = random.choice(sites)
        newSeq = subseq[(ix - bp):(ix+bp+1)]
        if not re.search("[NYRATCG]{21}", newSeq): # THIS CHANGED!
            #print(pos)
            sites.remove(ix)
            try:
                if c_direction == -1:
                    sites1.remove(ix)
                else:
                    sites2.remove(ix)
            except ValueError: 
                print("EYORE: Could not remove index {} from list of sites".format(ix))
                print("Newseq = {}".format(newSeq))
            continue
        if c_direction == -1:
            distance = window + bp - ix
        else:
            distance = bp + ix
        motif2 = full_motif(newSeq, ref)
        entry = {
            'chrom' : chrom,
            'pos' : pos,
            'motif' : newSeq,
            'cat': cat,
            'ref': ref,
            'window': window,
            'distance': distance,
            'motif2':motif2
        }
        newlist.append(entry)
        if c_direction == -1:
            sites1.remove(ix)
        else:
            sites2.remove(ix)
    return newlist



if __name__ == "__main__":
    #random.seed( 8675 ) # threeeeee ohhhhh niiiiii-eee-iiiiine
    random.seed( 1776 ) # round two
    main()
