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
# Copairs analysis which includes the comparison of standard CP data with the Phalloidin 400LS data.
#
# * **Normalization:** Negcon normalization
# * **mAP calculation:** The following notebook includes the data from mAP calculated as difference to controls as well as the mAP calculated as difference to other treatments.

from pathlib import Path

# %%
import pandas as pd

# %% [markdown]
# #### Reading the dataframes

# %% [markdown]
# #### mAP calculated as difference to controls

# %% [markdown]
# ##### Standard CP data

# %%
# Reading the results dataframe of standard plate
copairs_dir = Path("copairs_csv")
name_csv = {
    "CellPainting": "Result_Negcon_wrt_Controls_StandardCP",
    "Saguaro": "Result_Negcon_wrt_Controls_batch5",
    # "Saguaro": "PrecisionValues_with_MoA_allplates_Negcon_wrt_Controls_48and49",
    "CP_Batch3": "Result_Negcon_wrt_Controls_batch3",
    "Phalloidin": "Result_Negcon_wrt_Controls_Phalloidin400LS",
    "MitoBrilliant": "Result_Negcon_wrt_Controls_Tocris_mitobrilliant",
}
dfs = {
    name: pd.read_csv(copairs_dir / f"{csv_file}.csv")
    for name, csv_file in name_csv.items()
}


def find_col(df: pd.DataFrame, names):
    for name in names:
        if name in df.columns:
            return name


# Homogeneise BRD ID column
def rename_id_col(df: pd.DataFrame) -> pd.DataFrame:
    if not "BRD ID" in df.columns:
        name = find_col(df, ("Metadata_BRD ID", "Metadata_broad_sample.1"))
        assert name is not None, f"Broad ID col not found"
        df = df.rename(columns={name: "BRD ID"})
    return df


def average_cell_count(
    # df: pd.DataFrame, reference: str = "BRD ID", suffix: str = ""
    df: pd.DataFrame,
    suffix: str = "",
) -> pd.DataFrame:
    # id_col = find_col(df, ("BRD ID", "Metadata_BRD ID"))
    id_col = "BRD ID"
    # id_col = "Metadata_broad_sample"
    # cc_col = find_col(df, ("Metadata_Count_Cells",))
    mean_cell_count = (
        rename_id_col(df).groupby(id_col)["Metadata_Count_Cells"].mean().reset_index()
    )
    # df[f"Metadata_Count_Cells_{suffix}"] = df["Metadata_Count_Cells"] / 100
    mean_cell_count[f"Metadata_Count_Cells_{suffix}_norm"] = (
        mean_cell_count["Metadata_Count_Cells"] / 100
    )
    del mean_cell_count["Metadata_Count_Cells"]
    # df["Metadata_Batch_Name"] = suffix
    return mean_cell_count


# %% [markdown]
# ##### Phalloidin 400 LS - long stoke shifted actin

# %%


moa_name_csv = {
    "std_mito_pha": "PrecisionValues_with_MoA_allplates_cellcount_Negcon_wrt_Controls",
    "cp_saguaro": "PrecisionValues_with_MoA_allplates_Negcon_wrt_Controls_48and49",
}

moa_dfs = {
    name: pd.read_csv(copairs_dir / f"{csv_file}.csv")
    for name, csv_file in moa_name_csv.items()
}

# df = rename_id_col(dfs[name])
# df["Metadata_Count_Cells_norm"] = df["Metadata_Count_Cells"] / 100
for name in ("CellPainting", "Saguaro", "MitoBrilliant"):
    moa_dfs["cp_saguaro"] = pd.merge(
        rename_id_col(moa_dfs["cp_saguaro"]),
        average_cell_count(dfs[name], suffix=name),
        on="BRD ID",
    )

combined_df = pd.merge(*moa_dfs.values(), on="BRD ID")
combined_df["MoA"] = combined_df["MoA_x"]
combined_df["Common Name"] = combined_df["Common Name_x"]

# %% [markdown]
# Save only the relevant fields for plotting
# %%

combined_df[
    [
        x
        for x in combined_df.columns
        if "average_precision" in x
        or "Count_Cells" in x
        or "_vs_" in x
        or x == "MoA"
        or x == "Common Name"
    ]
].to_csv(
    copairs_dir / "PrecisionValues_with_MoA_Negcon_wrt_Control_allplates_cellcount.csv",
)
