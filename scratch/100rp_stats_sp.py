"""
Code for getting the deviance statistics at all positions within +/- 100bp window
"""
import pandas as pd
import statsmodels.api as sm
from pyfaidx import Fasta
from patsy import dmatrices
import statsmodels.formula.api as smf
import ray

def c_pos(subtype, chromosome):
  """Get the positions for controls for a given subtype"""
  input_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/controls/ALL/pos_files/"
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt"
  pos_list = pd.read_csv(f_name, header=0, names = ['pos'], usecols=['pos'], squeeze = True)
  return pos_list

def s_pos(subtype, chromosome):
  """Get the positions for singletons for a given subtype"""
  input_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/output/singletons/ALL/pos_files/"
  f_name = input_dir + subtype + "_" + str(chromosome) + ".txt"
  pos_list = pd.read_csv(f_name, header=0, names = ['pos'], usecols=['pos'], squeeze = True)
  return pos_list

ray.init(num_cpus=22)

@ray.remote
def get_count_table_control(chromosome, subtype, offset):
  """
  Get the count table for controls at a given relative position for a subtype.
  Note: this version has been made so as to be called in parallel via ray, but
  that it sucks in that each worker has to load the reference genome. Need to
  check if there's a python library that has the reference genome as on object
  that can maybe be passed around?
  """
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
  """
  Get the count table for singletons at a given relative position for a subtype.
  Note: this version has been made so as to be called in parallel via ray, but
  that it sucks in that each worker has to load the reference genome. Need to
  check if there's a python library that has the reference genome as on object
  that can maybe be passed around?
  """
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

def get_count_all(subtype, offset, status = "singleton"):
  """
  Get the singleton and control counts for a given relative position
  across all 22 autosomes. This version is parallelized via ray
  """
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
  n_s = sum(s_tab.singletons)
  n_c = sum(c_tab.controls)
  return {"dev":mod.deviance, "singletons":n_s, "controls":n_c, "offset":offset}

results = []
subtype = "cpg_GC_CG"
for offset in range(1, 1001):
  print(offset)
  results.append(fit_model_all(subtype, offset))
  results.append(fit_model_all(subtype, offset * -1))
  
final = pd.DataFrame.from_dict(results)
out_dir = "/net/snowwhite/home/beckandy/research/1000G_LSCI/scratch/100bp_sp/"
file_name = out_dir + subtype + ".csv"
final.to_csv(file_name, index = False)
