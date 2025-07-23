# needed modules: "pip install pandas statsmodels matplotlib numpy" if error raised
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.stats.proportion as smp
import statsmodels.api as sm

# load and prepare data
file_path = "FILENAME.csv"
data = pd.read_csv(file_path)

# check the column names
ma_columns = ["Study", "Species", "X", "N"]
data = data[ma_columns]

# calculating prevalence and 95% CI per study
data["Prevalence"] = data["X"] / data["N"]
data["CI_Lower"], data["CI_Upper"] = zip(*[smp.proportion_confint(x, n, alpha=0.05, method="wilson") for x, n in zip(data["X"], data["N"])])

# proper presentation with percent
data["Prevalence"] *= 100
data["CI_Lower"] *= 100
data["CI_Upper"] *= 100

# group by species
species_groups = data.groupby("Species")

# create containers for meta-analysis summaries
summary_species = []
summary_prev = []
summary_CI_lower = []
summary_CI_upper = []
# between-study variance
summary_tau2 = []
# Cochran's Q
summary_Q = []
# I²
summary_I2 = []

def meta_analysis_random_effects(x_events, n_total):
    """
    x_events: Number of positive cases (X)
    n_total: Number of total individuals tested (N)
        
    returns: 
        summary_p: pooled proportion
        ci_lower: lower bound of 95% CI
        ci_upper: upper bound of 95% CI
        tau2: between-study variance estimate
        Q: Cochran's Q statistic
        I2: I-squared statistic
    """

    # declaring numpy arrays
    x_events = np.array(x_events, dtype=float)
    n_total = np.array(n_total, dtype=float)
    # prevalence calculation 0-1
    p_i = x_events / n_total
    # Logit transform (theta_i) and approximate variance
    p_i = np.clip(p_i, 1e-8, 1-1e-8)
    # logit
    theta_i = np.log(p_i / (1 - p_i))
    # Var(logit(p)) approximation
    var_i = 1/x_events + 1/(n_total - x_events)
    # inverse-variance
    w_i = 1 / var_i
    # logit scale
    theta_fixed = np.sum(w_i * theta_i) / np.sum(w_i)
    # Cochran's Q
    Q = np.sum(w_i * (theta_i - theta_fixed)**2)
    df = len(theta_i) - 1
    # DerSimonian–Laird between-study variance
    tau2 = max(0, (Q - df) / (np.sum(w_i) - np.sum(w_i**2)/np.sum(w_i)))
    # random-effect weights
    w_star = 1 / (var_i + tau2)
    theta_star = np.sum(w_star * theta_i) / np.sum(w_star)
    # weighted SE under random-effects
    se_star = np.sqrt(1 / np.sum(w_star))
    # 95% CI on the logit scale
    ci_low_logit = theta_star - 1.96 * se_star
    ci_high_logit = theta_star + 1.96 * se_star
    # back-transform to proportion
    summary_p = 1 / (1 + np.exp(-theta_star))
    ci_lower_p = 1 / (1 + np.exp(-ci_low_logit))
    ci_upper_p = 1 / (1 + np.exp(-ci_high_logit))
    # I² statistic
    I2 = max(0, (Q - df) / Q) * 100 if Q > df else 0.0

    return summary_p, ci_lower_p, ci_upper_p, tau2, Q, I2


# meta-analysis for each species
for species, group in species_groups:
    x_events = group["X"].values
    n_total = group["N"].values
    
    # random-effects meta-analysis
    pooled_p, ciL, ciU, tau2, Qval, i2val = meta_analysis_random_effects(x_events, n_total)
    
    # percentage
    pooled_p *= 100
    ciL *= 100
    ciU *= 100
    
    # append results to the lists
    summary_species.append(species)
    summary_prev.append(pooled_p)
    summary_CI_lower.append(ciL)
    summary_CI_upper.append(ciU)
    summary_tau2.append(tau2)
    summary_Q.append(Qval)
    summary_I2.append(i2val)

# dataframe from dictionary for summary data
summary_df = pd.DataFrame({
    "Study": ["Summary"] * len(summary_species),
    "Species": summary_species,
    "Prevalence": summary_prev,
    "CI_Lower": summary_CI_lower,
    "CI_Upper": summary_CI_upper,
    "tau2": summary_tau2,
    "Q": summary_Q,
    "I2(%)": summary_I2
})

# order microorganisms by prevalence (largest first)
summary_df = summary_df.sort_values(by="Prevalence", ascending=True).reset_index(drop=True)

# forest plot
fig, ax = plt.subplots(figsize=(7, max(9, len(summary_df) * 0.5)))
y_positions = range(len(summary_df))

# summary plot with error bars
ax.errorbar(
    summary_df["Prevalence"], 
    y_positions,
    xerr=[
        summary_df["Prevalence"] - summary_df["CI_Lower"], 
        summary_df["CI_Upper"] - summary_df["Prevalence"]
    ],
    fmt='o', color='blue', label="Prevalence"
)

# make sure to set the max and min limits for the plot by changing ax.set_xlin
ax.set_yticks(y_positions)
ax.set_yticklabels(summary_df["Study"] + " - " + summary_df["Species"])
ax.set_xlabel("Prevalence (%)")
ax.set_title("Forest Plot of Bacterial Prevalence (Random-Effects)")
ax.set_xlim(-5, 60)

# This part adds pooled_prev summary effect vertical line
for i, row in summary_df.iterrows():
    pooled_prev = row["Prevalence"]
    ax.axvline(pooled_prev, linestyle="--", color="red", alpha=0.3)

plt.grid(True, linestyle="--", alpha=0.5)
plt.savefig('FILENAME.svg', format='svg', dpi=300, bbox_inches='tight')
plt.show()

# 6) Save summary table
summary_df.to_csv("FILENAME.csv", index=False)
