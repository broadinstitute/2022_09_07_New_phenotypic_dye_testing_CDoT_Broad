# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: base
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Copairs
# * **Details of the analysis in this notebook:**
# * **Data from :** CDoT
# * **Plates compared:**
#     * BR00122250 - stained with standard cellpaitning dyes
#     * BR00122246 - stained with Tocris Mitobrilliant dye along with the other CP dyes
#     * BR00122247 - stained with Phalloidin400LS dye along with other CP dyes
# * **Objective:** To understand the mAP of the plates stained with the new set of dyes.
# * **Normalization:** Negcon normalization
# * **mAP calculation:** mAP is calculated as difference to other treatments.

# %% [markdown]
# ### Modules import
# %%


import logging
from pathlib import Path

import numpy as np
import pandas as pd
from copairs.map import aggregate, run_pipeline
from tqdm import tqdm

logging.basicConfig(format="%(levelname)s:%(asctime)s:%(name)s:%(message)s")
logging.getLogger("copairs").setLevel(logging.INFO)

# %% [markdown]
# ### Reading the dataframe
# Batch 1 consists of plate with standard CP dye and other plate stained with Tocris Mitobrilliant dye

# %%
copairs_path = Path("copairs_csv")
gct_path = Path("gct")
batch_file = {
    str(v).split("/")[1].split("_")[-1].lower(): v
    for v in sorted(gct_path.rglob("*Batch*normalized_feature_select_batch.csv.gz"))
}
# %% [markdown]
# ### Analysis - Plate wise with respect to other treatments

# %% [markdown]
# #### Defining parameters to compute map

# %%
pert_col = "Metadata_broad_sample"

# %%
pos_sameby = [pert_col]
pos_diffby = []

neg_sameby = []
neg_diffby = [pert_col]
null_size = 10000

# %% [markdown]
# Rename batch names to their dye set
# %%
batch_dye = {
    "batch2": "Phalloidin400LS",
    "batch3": "Saguaro_then_CP",
    "batch5": "Saguaro",
}


dfs = {k: pd.read_csv(v) for k, v in batch_file.items()}
dye_data = {v: dfs[k] for k, v in batch_dye.items() if k in batch_dye}
for new_name, plate in (("Standard", "BR00122250"), ("MitoBrilliant", "BR00122246")):
    sub_df = dfs["batch1"]
    dye_data[new_name] = sub_df.loc[sub_df["Metadata_Plate"] == "BR00122250"]

# %% Analysis

agg_results = []
for dye, profile in dye_data.items():
    print(dye)
    metadata_names = [c for c in profile.columns if c.startswith("Metadata")]
    feature_names = [c for c in profile.columns if not c.startswith("Metadata")]
    feats = profile[feature_names].values
    dframe = profile[metadata_names]

    # To get MoA stuff remove negative controls
    feats = feats[dframe["Metadata_control_type"].isna(), :]
    dframe = dframe.loc(axis=0)[dframe["Metadata_control_type"].isna()]
    result = run_pipeline(
        dframe, feats, pos_sameby, pos_diffby, neg_sameby, neg_diffby, null_size
    )
    result.to_csv(copairs_path / f"stats_trt_{dye}.csv")

    agg_result = aggregate(result, sameby=pos_sameby, threshold=0.05)
    agg_result["Dye set"] = dye
    agg_results.append(agg_result)
    # agg_result.to_csv(copairs_path / f"agg_stats_trt_{dye}.csv")

# %% Combine all data sets
combined_df = pd.concat(agg_results, axis=0, ignore_index=True)

# %% [markdown]
# #### Adding metadata information to the combined_df
#

# %%
moa_metadata = pd.read_csv(copairs_path / "LC00009948_MoA_Common_Names.csv")
moa_metadata = moa_metadata.rename(columns={"BRD with batch": "Metadata_broad_sample"})


# %% [markdown]
# ##### Extracting BRD ID from BROAD sample name


# %%
def BRD_ID(i):
    if type(i) != float:
        ID = i.split("-")
        return ID[1]


combined_df["BRD ID"] = combined_df["Metadata_broad_sample"].map(BRD_ID)
moa_stats = pd.merge(combined_df, moa_metadata, on="BRD ID")

moa_stats.to_csv(copairs_path / f"moa_stats_5plates.csv", index=False)
