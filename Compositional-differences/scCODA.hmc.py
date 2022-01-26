#!/usr/bin/env python3
#-*- coding:utf-8 -*-


import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import pickle as pkl

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz

import sys

def main():

    #caste1, caste2, celltype = sys.argv[1:4]
    #caste1, caste2 = "Queen", "Gyne"
    #number_burnin, number_results = "10000", "800000"
    #caste1, caste2, num = sys.argv[1:4]

    #cell_counts = pd.read_table("num_batch.matrix.2.txt")
    #data_all = dat.from_pandas(cell_counts, covariate_columns=["Ant"])
    #data_all.obs["Condition"] = data_all.obs["Ant"].str.replace(r"_[0-9]+", "")

    #data_salm = data_all[data_all.obs["Condition"].isin([caste1, caste2])]
    #model_salm = mod.CompositionalAnalysis(data_salm, formula="Condition", reference_cell_type="automatic")
    #hmc_results = model_salm.sample_hmc(num_burnin=int(number_burnin), num_results=int(number_results))
    #da_results = model_salm.sample_hmc_da(num_burnin=int(number_burnin), num_results=int(number_results))
    #nuts_results = model_salm.sample_nuts()

    #hmc_f = caste1 + '_' + caste2 + '_auto_' + number_burnin + '_' + number_results + '.' + num + '.hmc'
    #hmc_results.save(hmc_f)
    #da_results.save(caste1 + '_' + caste2 + '_' + celltype + '_' + number_burnin + '_' + number_results + '.hmc_da')
    #nuts_results.save(caste1 + '_' + caste2 + '_' + celltype + '.nuts')

    with open(sys.argv[1],'rb') as f:
        hmc_results = pkl.load(f)
    hmc_results.set_fdr(est_fdr=0.2)
    alpha_df, beta_df = hmc_results.summary_prepare(est_fdr=0.2)
    beta_df.index = alpha_df.index
    all_df = pd.concat([alpha_df, beta_df], axis=1).rename(index = lambda s: sys.argv[1] + '_' + s)
    all_df.to_csv(sys.argv[1] + '.summary', sep='\t')


if __name__ == "__main__":
    main()
