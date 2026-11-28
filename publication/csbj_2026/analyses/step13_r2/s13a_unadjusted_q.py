"""Editor R2, point 1: q-values with and without the signal adjustment.
Re-runs of the five lists with code that also reports the unadjusted q (step13_r2/qvalues/*); checks that the
adjusted q and the kept terms are identical to the published runs (step09, step11), then summarizes.
Output: step13_r2/s13a_q_comparison.tsv, s13a_terms_HIV.tsv"""
import glob, os
import pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(H, "..")
OLD = {"HIV": "step09_end_to_end/HIV_dependency_factors.symbols_2026*", "König": "step11_three_screens/Konig_2008_*_2026*",
       "Zhou": "step11_three_screens/Zhou_2008_*_2026*", "Hypoxia": "step09_end_to_end/hypoxia-genes.symbols_2026*",
       "Vaccination": "step09_end_to_end/vaccination_NAKAYA_*_2026*"}
NEW = {"HIV": "HIV_dependency_factors", "König": "Konig_2008", "Zhou": "Zhou_2008", "Hypoxia": "hypoxia-genes", "Vaccination": "vaccination_NAKAYA"}
rows = []
for k in OLD:
    o = pd.read_csv(glob.glob(os.path.join(A, OLD[k], "*_term_qvalues.tsv"))[0], sep="\t")
    n = pd.read_csv(glob.glob(os.path.join(H, "qvalues", NEW[k] + "*_2026*", "*_term_qvalues.tsv"))[0], sep="\t")
    m = o.merge(n, on=["Library", "Term", "Iteration"], suffixes=("_old", ""))
    assert len(m) == len(o) == len(n), k
    assert (abs(m["q-value_old"] - m["q-value"]) < 1e-12).all() and (m.Kept_old == m.Kept).all(), k
    qr = m["q-value (no signal adjustment)"]
    # the unadjusted q is slightly lower for a few strongly enriched terms (Monte Carlo: the adjusted q of smaller sizes
    # comes from 1,000 lists per size step, the unadjusted q from the 10,000 lists of the input size); never above 0.05
    low = (qr < m["q-value"] - 1e-12); assert (m.loc[low, "q-value"] <= 0.05).all(), k
    rows.append(dict(list=k, selected=len(m), kept_adjusted_q_lt_0_10=int(m.Kept.sum()),
                     unadjusted_q_lt_0_10=int((qr < 0.10).sum()), unadjusted_q_lt_0_20=int((qr < 0.20).sum()),
                     kept_but_unadjusted_ge_0_10=int((m.Kept & (qr >= 0.10)).sum()),
                     unadjusted_lower_mc=int(low.sum()), max_unadjusted_q_of_kept=round(qr[m.Kept].max(), 3), median_ratio_unadj_adj=round((qr / m["q-value"])[m.Kept].median(), 2)))
    if k == "HIV":
        m[["Library", "Iteration", "Term", "p-value", "q-value", "q-value (no signal adjustment)", "Kept"]].to_csv(os.path.join(H, "s13a_terms_HIV.tsv"), sep="\t", index=False)
d = pd.DataFrame(rows); d.to_csv(os.path.join(H, "s13a_q_comparison.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250); print(d.to_string(index=False))
