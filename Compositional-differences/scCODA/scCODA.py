#!/usr/bin/env python3
# -*- coding:utf-8 -*-


# Setup
import warnings

warnings.filterwarnings("ignore")

import pandas as pd
import pickle as pkl

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz


def scCODA_analysis(caste1, caste2, number_burnin, number_results, cellcount_f, reference):
    """Data preparation : choose a dataset used for scCODA analysis"""
    cell_counts = pd.read_table(cellcount_f)

    # Use M.pha. dataset as example
    data_all = dat.from_pandas(cell_counts, covariate_columns=["Ant"])
    data_all.obs["Condition"] = data_all.obs["Ant"].str.replace(r"_[0-9]+", "")
    # Use Dmel age dataset as example
    #data_all = dat.from_pandas(cell_counts, covariate_columns=["Dmel"])
    #data_all.obs["Condition"] = data_all.obs["Dmel"].str.replace(r"_DGRP_r[1-3]|_w1118_r[1-3]", "")
    # Use Dmel head dataset as example
    #data_all = dat.from_pandas(cell_counts, covariate_columns=["Dmelhead"])
    #data_all.obs["Condition"] = data_all.obs["Dmelhead"].str.replace(r"_FCA[0-9]+", "")

    # Choose any pairwise in the dataset to assess the cell compositional differences
    data_salm = data_all[data_all.obs["Condition"].isin([caste1, caste2])]

    # Model setup and automatically select an appropriate cell type as the reference category
    # or use a prespecified reference cell type to identify compositional changes for the rest cell type
    model_salm = mod.CompositionalAnalysis(data_salm, formula="Condition", reference_cell_type=reference)

    # Choose a suitable MCMC sampling method and run
    # There are three different MCMC sampling methods avalible for scCODA:
    # 1) Hamiltonian Monte Carlo (HMC) sampling: sample_hmc()
    # 2) HMC sampling with Dual-averaging step size adaptation (Nesterovv, 2009): sample_hmc_da()
    # 3) No-U-Turn sampling (Hoffman and Gelman, 2014): sample_nuts
    hmc_results = model_salm.sample_hmc(num_burnin=int(number_burnin), num_results=int(number_results))
    #da_results = model_salm.sample_hmc_da(num_burnin=int(number_burnin), num_results=int(number_results))
    #nuts_results = model_salm.sample_nuts()

    out_name = caste1 + '_' + caste2 + '_auto_' + number_burnin + '_' + number_results
    hmc_results.save(out_name + '.hmc')
    #da_results.save(out_name + '.hmc_da')
    #nuts_results.save(out_name + '.nuts')

    # Saving MCMC results
    with open(out_name + '.hmc', 'rb') as f:
        hmc_results = pkl.load(f)

    # Set a desired False discovery rate (FDR) of 0.2 to reveal effects on caste1 and caste2 cells.
    hmc_results.set_fdr(est_fdr=0.2)
    alpha_df, beta_df = hmc_results.summary_prepare(est_fdr=0.2)
    beta_df.index = alpha_df.index
    all_df = pd.concat([alpha_df, beta_df], axis=1).rename(index=lambda s: out_name + '.hmc' + '_' + s)
    all_df.to_csv(out_name + '.hmc' + '.summary', sep='\t')


def main():

    # Choose Mpha as example
    cellcount_f = "Mpha.rep.matrix.txt"
    caste_list = ["Worker", "Queen", "Gyne", "Male"]
    for caste1 in caste_list:
        for caste2 in caste_list:
            scCODA_analysis(caste1, caste2, number_burnin="10000", number_results="800000",
                            cellcount_f, reference="automatic")
    # Choose Dmel as example
    # cellcount_f = "Dmel.age.rep.matrix.txt"
    # cellcount_f, reference = "Dmel.head.rep.matrix.txt", "c0"
    #scCODA_analysis("Male", "Female", number_burnin="10000", number_results="800000", cellcount_f, reference="automatic")

if __name__ == "__main__":
    main()
