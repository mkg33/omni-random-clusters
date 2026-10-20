#!/usr/bin/env python3
# Random clustering baseline: each cell gets a random cluster from 1..k (k = true number of clusters)

import argparse
import json

import numpy as np
import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument("--output_dir", required=True)
parser.add_argument("--name", required=True)
parser.add_argument("--data.h5ad", dest="h5ad", required=True)
parser.add_argument("--data.clusters_truth_num", dest="k", required=True)
parser.add_argument("--seed", type=int, required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.h5ad)
k = int(open(args.k).read())
labels = np.random.default_rng(args.seed).integers(1, k + 1, size=adata.n_obs)

# the metrics module also needs a PCA (for the structure metrics)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=30)

out = f"{args.output_dir}/{args.name}"
pd.DataFrame({"cell_id": adata.obs_names, "louvain": labels, "leiden": labels}).to_csv(
    f"{out}.clusters.tsv", sep="\t", index=False
)
pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names, columns=[f"PC{i}" for i in range(1, 31)]).to_csv(
    f"{out}.pca.tsv", sep="\t", index_label="cell_id"
)
# nothing to time or tune here
for f in ["timings", "clustering_info"]:
    json.dump({}, open(f"{out}.{f}.json", "w"))
