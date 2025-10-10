import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file = "/Users/JB/Rotation_bkslab/250203_alphafold3/LSD_AF3.xlsx" #Change to your xlsx file, must have 'L-pLDDT score' and 'Active' columns
df = pd.read_excel(file, sheet_name="D4")

def rolling_hit_rate(
    df,
    score_col="L-pLDDT score",
    label_col="Active",
    window=50,
    higher_is_better=True,
    min_periods=None
):
    
    dd = df[[score_col, label_col]].dropna().copy()
    dd[label_col] = dd[label_col].astype(int)

    # sort so that "best" scores are on the left of the plot
    dd = dd.sort_values(score_col, ascending=not higher_is_better).reset_index(drop=True)

    if min_periods is None:
        min_periods = max(25, window // 5)

    # rolling counts and means over the ranked list
    rcount = dd[label_col].rolling(window, min_periods=min_periods).count()
    rmean  = dd[label_col].rolling(window, min_periods=min_periods).mean()

    # Wilson 95% confidence interval for a proportion
    z = 1.96
    n = rcount.astype(float)
    p = rmean.fillna(0.0)
    denom  = 1.0 + (z**2)/n.replace(0, np.nan)
    center = (p + (z**2)/(2*n.replace(0, np.nan))) / denom
    halfw  = (z*np.sqrt((p*(1.0-p)/n.replace(0, np.nan)) + (z**2)/(4*n.replace(0, np.nan)**2))) / denom
    lo95 = (center - halfw).clip(lower=0.0)
    hi95 = (center + halfw).clip(upper=1.0)

    out = pd.DataFrame({
        "score": dd[score_col].to_numpy(),
        "hit_rate": rmean.to_numpy(),
        "lo95": lo95.to_numpy(),
        "hi95": hi95.to_numpy(),
        "n": n.to_numpy()
    })

    # Drop positions before the first complete-ish window
    out = out[~np.isnan(out["hit_rate"])].reset_index(drop=True)
    return out


def plot_rolling_hit_rate(
    roll_df,
    baseline=None,
    method_name="AF3",
    color="tab:red",
    score_label="Score",
    title="Rolling hit-rate vs. score"
):
    """
    Plot rolling hit-rate with Wilson 95% CI and optional baseline line.
    """
    plt.rcParams.update({
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 18,
        "legend.fontsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12
    })

    fig = plt.figure(figsize=(8,6))
    x = roll_df["score"].to_numpy()
    y = roll_df["hit_rate"].to_numpy()
    lo = roll_df["lo95"].to_numpy()
    hi = roll_df["hi95"].to_numpy()

    plt.plot(x, y, color=color, linewidth=2, label=f"{method_name} (rolling)")
    plt.fill_between(x, lo, hi, color=color, alpha=0.2, linewidth=0)
    baseline = df["Active"].mean()

    if baseline is not None:
        plt.axhline(baseline, linestyle="--", color="gray", linewidth=2,
                    label=f"Random baseline = {baseline:.2f}")

    plt.xlabel(score_label)
    plt.ylabel("Hit rate (actives fraction)")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()