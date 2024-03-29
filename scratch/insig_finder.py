import pandas as pd
import statsmodels.api as sm
from pyfaidx import Fasta
from patsy import dmatrices
import statsmodels.formula.api as smf
import ray

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

## TESTING RAY
ray.init(num_cpus=26)

@ray.remote
def get_count_table_control(chromosome, subtype, offset):
  ref_file = "/net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  fasta_obj = Fasta(ref_file)
  control_pos = c_pos(subtype, chromosome)
  seq = fasta_obj["{}{}".format("chr", chromosome)]
  seqstr = seq[0:len(seq)].seq
  results = {"A":0, "C":0, "G":0, "T":0}
  for index, value in control_pos.items():
    ix = value - 1 + offset
    nuc = seqstr[ix]
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

@ray.remote
def get_count_table_singletons(chromosome, subtype, offset):
  singleton_pos = s_pos(subtype, chromosome)
  ref_file = "/net/snowwhite/home/beckandy/research/1000G_LSCI/reference_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  fasta_obj = Fasta(ref_file)
  seq = fasta_obj["{}{}".format("chr", chromosome)]
  seqstr = seq[0:len(seq)].seq
  results = {"A":0, "C":0, "G":0, "T":0}
  for index, value in singleton_pos.items():
    ix = value - 1 + offset
    nuc = seqstr[ix]
    if nuc in results.keys():
      results[nuc] += 1
  return(results)

# futures = [get_count_table_control.remote(i, "AT_TA", -2) for i in range(1,23)]
# df = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)

### TEST THIS!
def get_count_all(subtype, offset, status = "singleton"):
  if status == "singleton":
    futures = [get_count_table_singletons.remote(i, subtype, offset) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  else:
    futures = [get_count_table_control.remote(i, subtype, offset) for i in range(1,23)]
    results = pd.DataFrame.from_dict(ray.get(futures)).sum(axis=0)
  return(results)

def fit_model_all(subtype, offset):
  s_tab = get_count_all(subtype, offset, status = "singleton")
  s_tab = s_tab.reset_index(level=0)
  s_tab.columns = ['nuc', 'singletons']
  c_tab = get_count_all(subtype, offset, status = "control")
  c_tab = c_tab.reset_index(level=0)
  c_tab.columns = ['nuc', 'controls']
  df = pd.DataFrame.merge(s_tab, c_tab, on='nuc')
  df = pd.melt(df, id_vars = 'nuc', var_name = 'status', value_name = "n")
  mod = smf.glm("n ~ status + nuc", df, family = sm.families.Poisson()).fit()
  return mod.deviance


def find_nonsig_all(subtype, max_d = 1000, sig_stat = 7.815):
  sig = True
  offset = 1
  sign = 1
  while sig and offset <= max_d:
    if fit_model_all(subtype, offset) < sig_stat:
      sig = False
      sign = 1
    elif fit_model_all(subtype, -1 * offset) < sig_stat:
      sig = False
      sign = -1
    else:
      offset += 1
  return offset*sign
