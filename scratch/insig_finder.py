import pandas as pd
import statsmodels as sm
from pyfaidx import Fasta
from patsy import dmatrices
import statsmodels.formula.api as smf

ref_file = "/net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
fasta_obj = Fasta(ref_file)

def c_pos(subtype, chromosome):
  input_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/controls/ALL/pos_files/"
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt"
  pos_list = pd.read_csv(f_name, header=0, names = ['pos'], usecols=['pos'], squeeze = True)
  return pos_list

def s_pos(subtype, chromosome):
  input_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/singletons/ALL/pos_files/"
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt"
  pos_list = pd.read_csv(f_name, header=0, names = ['pos'], usecols=['pos'], squeeze = True)
  return pos_list

def get_count_table(pos_list, fasta_obj, chromosome, offset):
  seq = fasta_obj["{}{}".format("chr", chromosome)]
  seqstr = seq[0:len(seq)].seq
  results = {"A":0, "C":0, "G":0, "T":0}
  for index, value in pos_list.items():
    ix = value - 1 + offset
    nuc = seqstr[ix]
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

### TEST THIS!
def get_count_all(subtype, pos_list, fasta_obj, offset):
  results = {"A":0, "C":0, "G":0, "T":0}
  for chrom in range(1,23):
    seq = fasta_obj["{}{}".format("chr", chrom)]
    seqstr = seq[0:len(seq)].seq
    for index, value in pos_list.items():
      ix = value - 1 + offset
      nuc = seqstr[ix]
      if nuc in results.keys():
        results[nuc] += 1
  return(results)

def fit_model(subtype, chromosome, fasta_obj, offset):
  s_tab = get_count_table(s_pos(subtype, chromosome), fasta_obj, chromosome, offset)
  c_tab = get_count_table(c_pos(subtype, chromosome), fasta_obj, chromosome, offset)
  df = pd.DataFrame.merge(pd.DataFrame.from_dict(s_tab, orient = 'index', columns = ['singletons']), 
  pd.DataFrame.from_dict(c_tab, orient = 'index', columns = ['controls']), left_index= True, right_index=True)
  df.reset_index(inplace=True)
  df = df.rename(columns = {'index':'nuc'})
  df = pd.melt(df, id_vars = 'nuc', var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ status + nuc", df, family = sm.families.Poisson()).fit()
  return mod.deviance

def find_nonsig(subtype, chromosome, fasta_obj, max_d = 1000, sig_stat = 7.815):
  sig = True
  offset = 1
  sign = 1
  while sig and offset <= max_d:
    if fit_model(subtype, chromosome, fasta_obj, offset) < sig_stat:
      sig = False
      sign = 1
    elif fit_model(subtype, chromosome, fasta_obj, -1 * offset) < sig_stat:
      sig = False
      sign = -1
    else:
      offset += 1
  return offset*sign
