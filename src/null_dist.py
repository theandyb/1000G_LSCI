from pyfaidx import Fasta
import regex as re
import random
import pandas as pd
import argparse
from Bio.Seq import Seq
from collections import defaultdict

def sample_control(pos, ref, seqstr, window = 150):
  cand_sites = []
  control_sites = []
  while(len(cand_sites) < 5):
    seg_lb = (pos - 1) - window
    seg_ub = (pos - 1) + window + 1
    subseq = seqstr[seg_lb:seg_ub]
    cand_sites = [m.start() for m in re.finditer(ref, subseq)]
    if window in cand_sites:
      cand_sites.remove(window) # remove singleton site
    window = window + 50
  # revert window to acctual size once candidates identified
  window -= 50
  while(len(control_sites) < 5):
    ix = random.choice(cand_sites)
    new_pos = ix - window + pos
    control_sites.append(new_pos)
    cand_sites.remove(ix)
  return(control_sites)
  

def main(ref_prefix = "chr"):
  parser = argparse.ArgumentParser(description="Sample control distribution.")
  parser.add_argument("-s", "--singleton", help="Location of singleton file", required=True)
  parser.add_argument("-f", "--fasta", help="FASTA file to grab sequence from", required=True)
  parser.add_argument("-o", "--output", help="Path to output", required=True)
  parser.add_argument("chrom", help="Chromosome we are sampling from")
  args = parser.parse_args()
  
  ref_file = args.fasta 
  singleton_file = args.singleton
  chrom = args.chrom
  out_dir = args.output
  
  singleton_table = {}
  control_table = {}
  
  print("Reading fasta...")
  fasta_obj = Fasta(ref_file)
  seq = fasta_obj["{}{}".format(ref_prefix, chrom)] # hg37 does not have chr prefix, but hg38 does
  seqstr = seq[0:len(seq)].seq
  print("FASTA read!")
  with open(singleton_file) as fp:
    line = fp.readline()
    while line:
      content = line.strip().split(",")
      pos = int(content[1])
      ref = content[6]
      motif = content[7][:21]
      cat = content[8]
      if not cat in singleton_table:
        singleton_table[cat] = defaultdict(int)
        control_table[cat] = defaultdict(int)
      # Sample control sites
      control_sites = sample_control(pos, ref, seqstr)
      # Add counts to corresponding tables
      singleton_table[cat][seqstr[(pos - 1) + 1000]] += 1
      for x in control_sites:
        control_table[cat][seqstr[(x - 1) + 1000]] += 1
      line = fp.readline()
  for cat_type in singleton_table:
    control_dict = dict(control_table[cat_type])
    singleton_dict = dict(singleton_table[cat_type])
    df_c = pd.DataFrame.from_dict(control_dict, orient='index', columns = ["Controls"])
    df_s = pd.DataFrame.from_dict(singleton_dict, orient='index', columns = ["Singletons"])
    
    file_s = out_dir + "/singleton_" + str(chrom) + "_" + cat_type + ".csv"
    file_c = out_dir + "/control_" + str(chrom) + "_" + cat_type + ".csv"
  
    df_c.to_csv(file_c, index = True, index_label = "Nuc", header=True)
    df_s.to_csv(file_s, index = True, index_label = "Nuc", header=True)
  return(0)

if __name__ == "__main__":
    #random.seed( 8675 ) # threeeeee ohhhhh niiiiii-eee-iiiiine
    random.seed( 1776 ) # round two
    main()
