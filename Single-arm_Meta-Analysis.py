# needed modules: "pip install pandas statsmodels matplotlib numpy" if error raised
import pandas as pd
import numpy as np
from statsmodels.stats.meta_analysis import combine_effects
import matplotlib.pyplot as plt

# create a dataframe from a labelled (CAP and HAP) csv file
df = pd.read_csv("S-A M-A.csv")

def logit_transform(x, n):
    """Compute logit and variance for meta-analysis."""
    valid = (x > 0) & (n > x)
    yi = np.log(x / (n - x))
    vi = 1 / x + 1 / (n - x)
    return yi[valid], vi[valid], valid

def run_single_arm_meta(x, n, label="CAP"):
    yi, vi, mask = logit_transform(x, n)
    res = combine_effects(yi, vi, method_re="dl")
    return res, df.loc[mask]

# run separately for HAP and CAP
res_cap, df_cap = run_single_arm_meta(df["x_CAP"], df["n_CAP"], label="CAP")
res_hap, df_hap = run_single_arm_meta(df["x_HAP"], df["n_HAP"], label="HAP")

def fused_forest_plot(df_cap, res_cap, df_hap, res_hap):
    # combine CAP and HAP
    df_cap = df_cap.copy()
    df_cap["label"] = "CAP"
    df_cap["x"] = df_cap["x_CAP"]
    df_cap["n"] = df_cap["n_CAP"]

    df_hap = df_hap.copy()
    df_hap["label"] = "HAP"
    df_hap["x"] = df_hap["x_HAP"]
    df_hap["n"] = df_hap["n_HAP"]

    df_all = pd.concat([df_cap, df_hap], ignore_index=True)

    # compute logits and CIs
    df_all["yi"] = np.log(df_all["x"] / (df_all["n"] - df_all["x"]))
    df_all["vi"] = 1 / df_all["x"] + 1 / (df_all["n"] - df_all["x"])
    df_all["prop"] = np.exp(df_all["yi"]) / (1 + np.exp(df_all["yi"]))
    df_all["ci_low"] = np.exp(df_all["yi"] - 1.96 * np.sqrt(df_all["vi"])) / (1 + np.exp(df_all["yi"] - 1.96 * np.sqrt(df_all["vi"])))
    df_all["ci_upp"] = np.exp(df_all["yi"] + 1.96 * np.sqrt(df_all["vi"])) / (1 + np.exp(df_all["yi"] + 1.96 * np.sqrt(df_all["vi"])))

    # for plotting, create labels
    df_all["StudyLabel"] = df_all["Study"] + " (" + df_all["label"] + ")"
    df_all = df_all.sort_values(by="label", ascending=False).reset_index(drop=True)

    # get pooled values
    pooled = {
        "CAP": np.exp(res_cap.summary_frame().loc["random effect", "eff"]) / (1 + np.exp(res_cap.summary_frame().loc["random effect", "eff"])),
        "HAP": np.exp(res_hap.summary_frame().loc["random effect", "eff"]) / (1 + np.exp(res_hap.summary_frame().loc["random effect", "eff"])),
    }
    ci_low = {
        "CAP": np.exp(res_cap.summary_frame().loc["random effect", "ci_low"]) / (1 + np.exp(res_cap.summary_frame().loc["random effect", "ci_low"])),
        "HAP": np.exp(res_hap.summary_frame().loc["random effect", "ci_low"]) / (1 + np.exp(res_hap.summary_frame().loc["random effect", "ci_low"])),
    }
    ci_upp = {
        "CAP": np.exp(res_cap.summary_frame().loc["random effect", "ci_upp"]) / (1 + np.exp(res_cap.summary_frame().loc["random effect", "ci_upp"])),
        "HAP": np.exp(res_hap.summary_frame().loc["random effect", "ci_upp"]) / (1 + np.exp(res_hap.summary_frame().loc["random effect", "ci_upp"])),
    }

    # plot
    fig, ax = plt.subplots(figsize=(7, 0.5 * len(df_all) + 2))

    color_map = {"CAP": "blue", "HAP": "green"}
    for i, row in df_all.iterrows():
        ax.errorbar(row["prop"], i,
                    xerr=[[row["prop"] - row["ci_low"]], [row["ci_upp"] - row["prop"]]],
                    fmt='o', capsize=4, color=color_map[row["label"]])

    # add pooled lines
    for label, y in zip(["CAP", "HAP"], [-1, -2]):
        ax.axvline(pooled[label], linestyle="--", color=color_map[label], label=f"{label} pooled")
        ax.fill_betweenx([y - 0.3, y + 0.3], ci_low[label], ci_upp[label], color=color_map[label], alpha=0.3)

    ax.set_yticks(range(len(df_all)))
    ax.set_yticklabels(df_all["StudyLabel"])
    ax.set_xlabel("Prevalence")
    ax.set_title("HAP and CAP Prevalence")
    ax.invert_yaxis()
    ax.legend()
    plt.tight_layout()
    # optional save in svg format for further editting
    plt.savefig('S-A_M-A.svg', format='svg', dpi=300, bbox_inches='tight')
    plt.show()

fused_forest_plot(df_cap, res_cap, df_hap, res_hap)
