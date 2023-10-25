# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.1
#   kernelspec:
#     display_name: base
#     language: python
#     name: python3
# ---

# %% [markdown]
# ### figures - standard cp plate and the phalloidin400ls
#
# data were negcon normalized
### modules import
# %%

from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

# %% [markdown]
# ## map calculated as difference to controls
# %% [markdown]
# ##### reading the dataframe
# %%
copairs_dir = Path("copairs_csv")

combined_moa_cellcount_df = pd.read_csv(
    copairs_dir / "PrecisionValues_with_MoA_Negcon_wrt_Control_allplates_cellcount.csv"
)

# %% [code]
# set output folder
figs_dir = Path("figures") / "poster"
save_figures_externally = True

figs_dir.mkdir(exist_ok=True)


def quick_save(fig, fig_name: str, save_externally: bool = save_figures_externally):
    if save_externally:
        fig.write_html(str(figs_dir / f"{fig_name}"))
    else:
        fig.show("notebook")


# %% [markdown]
# ### comparison of mean average precision
#

# %% [markdown]
# #####  standard cellpainting data vs phalloidin 400ls
#

# %%
actin_fig = px.scatter(
    combined_moa_cellcount_df,
    x=combined_moa_cellcount_df["average_precision_std"],
    y=combined_moa_cellcount_df["average_precision_act"],
    labels={
        "average_precision_std": "mean average precision - standard cellpainting dyes",
        "average_precision_act": "mean average preicison - <br> phalloidin 400ls (long-stoke shifted)",
    },
    color=combined_moa_cellcount_df["MoA"],
)
actin_fig.update_layout(legend=dict(orientation="h"), height=800, width=1000)
quick_save(actin_fig, "actin_fig.html")

# %% [markdown]
# ### mean average precision values of all compounds
#

# %%
combined_box_plot = go.Figure()

field_name = (
    ("average_precision_std", "Standard CellPainting dyes"),
    ("average_precision_act", "Phalloidin 400LS"),
    ("average_precision_mito", "MitoBrilliant"),
    ("mean_average_precision_batch3", "Saguaro + CellPainting"),
    ("mean_average_precision_batch5", "Saguaro"),
)
for field, name in field_name:
    combined_box_plot.add_trace(
        go.Box(
            y=combined_moa_cellcount_df[field],
            name=name,
            boxpoints="all",
            hovertext=combined_moa_cellcount_df["MoA"]
            + "-"
            + combined_moa_cellcount_df["Common Name"],
        )
    )
combined_box_plot.update_layout(
    height=800,
    width=1000,
    font_family="arial",
    font=dict(size=14, color="black"),
    boxmode="group",
    yaxis_title="mean average precision",
)
quick_save(combined_box_plot, "combined_box_plot.html")

# %%  [markdown]
# Plot cell count vs mAP for each dye
# %%

import seaborn as sns

# %% [markdown]
# ### mean average precision values -  moa
#
# the size of the markers represent the average number of cells present in the replicates. the number of cells were normalized by dividing the actual number by 100 for easier plotting.

# %%
scatter_plot = go.Figure()

map_field_name = (
    (
        "average_precision_std",
        "Metadata_Count_Cells_Std_norm",
        "Standard CellPainting dyes",
    ),
    ("average_precision_act", "Metadata_Count_Cells_act_norm", "Phalloidin 400LS"),
    (
        "mean_average_precision_batch3",
        "Metadata_Count_Cells_CellPainting_norm",
        "Saguaro + CellPainting",
    ),
    ("mean_average_precision_batch5", "Metadata_Count_Cells_Saguaro_norm", "Saguaro"),
    # ("average_precision_mito", "MitoBrilliant"),
)
for map_field, cc_field, name in map_field_name:
    scatter_plot.add_trace(
        go.Scatter(
            # x=combined_moa_cellcount_df["MoA"],
            x=combined_moa_cellcount_df[cc_field],
            y=combined_moa_cellcount_df[map_field],
            # hovertext=[combined_moa_cellcount_df[map_field]],
            hovertext=[combined_moa_cellcount_df["MoA"]],
            mode="markers",
            name=name,
            marker_size=combined_moa_cellcount_df[cc_field],
        )
    )

scatter_plot.update_layout(
    height=1000,
    width=1500,
    font_family="arial",
    font=dict(size=14, color="black"),
    boxmode="group",
    yaxis_title="mean average precision",
    legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
)
scatter_plot.update_xaxes(tickangle=90, categoryorder="total ascending")
# quick_save(scatter_plot, "scatter_plot.html")
quick_save(scatter_plot, "scatter_cc_vs_map.html")

# %% [markdown]
# ### difference in mean average precision values
# the negative values indicate the better performance of phalloidin 400ls

# %%
dfs = []
for map_field, cc_field, name in map_field_name:
    tmp_df = pd.DataFrame()
    tmp_df["Cell Count"] = df[cc_field]
    tmp_df["Mean Average Precision"] = df[map_field]
    tmp_df["Dye Set"] = name
    dfs.append(tmp_df)
dye_sets_df = pd.concat(dfs, axis=0)

ax = sns.histplot(dye_sets_df, x="Cell Count", hue="Dye Set", element="step")
ax.spines[["right", "top"]].set_visible(False)
sns.move_legend(ax, loc="upper left")
import matplotlib.pyplot as plt

plt.savefig(figs_dir / "cell_count_distribution.png", dpi=200)
plt.close()
# %%
fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x=combined_moa_cellcount_df["MoA"],
        y=combined_moa_cellcount_df["std_vs_act"],
        mode="markers",
        hovertext=combined_moa_cellcount_df["Common Name"],
    )
)
fig.update_layout(
    height=1000, width=1700, font_family="arial", font=dict(size=14, color="black")
)

fig.update_yaxes(title="difference in <br> mean average precision")
fig.update_xaxes(categoryorder="total ascending")
quick_save(fig, "dmAP_performance.html")
