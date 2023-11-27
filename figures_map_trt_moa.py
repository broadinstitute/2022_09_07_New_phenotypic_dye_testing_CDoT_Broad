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
# ### Figures
# * **Normalization:** Negcon normalization
# * **mAP calculation:** mAP is calculated as difference to other treatments.

# %%
### Modules import
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt

# %% [markdown]
# #### Reading the dataframes

# %%
copairs_path = Path("copairs_csv")
stats_moa = pd.read_csv(
    # copairs_path / "PrecisionValues_with_MoA_allplates_Negcon_Trmts.csv"
    copairs_path
    / "moa_stats_5plates.csv",
)

# %% [markdown]
# ### Comparison of Mean average precision
#
# %% Compare all dyes against each other
# fig = sns.pairplot(data=stats_moa, vars=["Dye set"],hue="Dye set")
figs_path = Path("figures/")
sns.set_style("whitegrid")
fig = sns.boxenplot(
    data=stats_moa.sort_values("mean_average_precision"),
    y="Dye set",
    hue="Dye set",
    x="mean_average_precision",
)
plt.title("Dye biological retrievability, based on other treatments")
plt.tight_layout()
plt.savefig(figs_path / "boxenplot_dye_set_map.png", dpi=200)
plt.close()
# %% make each dye a column and compare each one of them separately
# %% List retrievability of each MoA group

moa_map_table = stats_moa.pivot_table(
    values="mean_average_precision", index="MoA", columns="Dye set"
)
sns.set_style("whitegrid")
ax = sns.pairplot(data=moa_map_table)
# ax = sns.catplot(
#     data=stats_moa.sort_values("mean_average_precision"),
#     y="MoA",
#     hue="Dye set",
#     x="mean_average_precision",
#     alpha=0.5,
#     col="Dye set",
# )
plt.title("Biological retrievability")
# sns.move_legend(ax, loc="lower left", bbox_to_anchor=(0, 1))
plt.tight_layout()
plt.savefig(figs_path / "moa_dye_set_map_pairplot.png", dpi=200)
plt.close()

# %% Difference between all vs saguaro


dsaguaro = moa_map_table.subtract(moa_map_table["Saguaro"], axis=0).drop(
    "Saguaro", axis=1
)


melted_dsaguaro = dsaguaro.melt(ignore_index=False).reset_index().sort_values("value")
melted_dsaguaro["Mechanism of Action (MoA)"] = melted_dsaguaro.MoA.str.replace(
    "inhibitor", "inh."
)
ax = sns.catplot(
    data=melted_dsaguaro,
    x="value",
    y="Mechanism of Action (MoA)",
    hue="Dye set",
    alpha=0.7,
)
ax.set_yticklabels(size=7)
# sns.move_legend(ax, loc="lower left", bbox_to_anchor=(-0.8, 0))
sns.move_legend(ax, loc="lower left")

plt.xlabel("mAP(X) - mAP(Saguaro)")
plt.tight_layout()
plt.savefig(figs_path / "catplot_moa_dsaguaro.png", dpi=200)
plt.title("mAP of Dyes - Saguaro ")
plt.close()

# %% Above p threshold
# %% Cell size?
# Calculate the MoAs where Saguaro over and underperforms
