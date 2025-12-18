# Explanation: Stratum and CDF

## What is a Stratum?

A **stratum** (plural: strata) is a subgroup or category that you divide your data into. Think of it as a "bin" or "bucket" that groups similar observations together.

### In the Context of iGEA P-Value Correction

In our case, a stratum is a **combination of conditions** that affect p-values. We group permutation results that have the same:
- Library (e.g., "GO BP")
- Iteration bin (e.g., "2-3")
- Term size bin (e.g., "50-100")
- Overlap size bin (e.g., "4")

### Example

Let's say you have a p-value from:
- Library: "GO BP"
- Iteration: 2
- Term size: 75 genes
- Overlap: 4 genes

This would belong to the stratum:
```
(GO BP, iteration_bin=2, term_size_bin="50-100", overlap_bin="4")
```

All permutation results with these same characteristics are grouped together in this stratum.

### Why Use Strata?

Because p-values behave differently under different conditions:
- P-values in iteration 1 are typically smaller than in iteration 5
- P-values in GO BP library are different from Reactome
- P-values for large terms (300+ genes) are different from small terms (<50 genes)

By creating separate strata, we can build **separate null distributions** for each combination, which gives us more accurate corrections.

### Visual Example

```
All Permutation Results (1.25M)
│
├── Stratum 1: (GO BP, Iteration 1, Term <50, Overlap 3)
│   └── 15,234 p-values → Build null distribution #1
│
├── Stratum 2: (GO BP, Iteration 1, Term <50, Overlap 4)
│   └── 12,456 p-values → Build null distribution #2
│
├── Stratum 3: (GO BP, Iteration 2, Term 50-100, Overlap 4)
│   └── 8,923 p-values → Build null distribution #3
│
└── ... (527 total strata with enough data)
```

---

## What is a CDF?

**CDF** stands for **Cumulative Distribution Function**. It tells you: "What fraction of values are less than or equal to a given value?"

### Simple Example

Imagine you have 100 p-values from permutations:
```
[0.001, 0.003, 0.005, 0.008, 0.012, 0.015, 0.020, 0.025, ...]
```

The CDF answers questions like:
- "What fraction of p-values are ≤ 0.01?" → Answer: 30% (30 out of 100)
- "What fraction of p-values are ≤ 0.05?" → Answer: 100% (all of them)

### Visual Example

```
P-value distribution:     CDF (Cumulative):
                          
0.05 |                   1.0 |                    ████████████
     |                        |                  █
0.04 |                        |                █
     |                        |              █
0.03 |                        |            █
     |                        |          █
0.02 |                        |        █
     |                        |      █
0.01 |                        |    █
     |                        |  █
0.00 |________________       0.0 |________________
     0    20   40   60           0.00  0.01  0.02  0.03  0.04  0.05
     Count of p-values              P-value
```

The CDF curve shows:
- At p=0.01, CDF ≈ 0.3 (30% of values ≤ 0.01)
- At p=0.02, CDF ≈ 0.6 (60% of values ≤ 0.02)
- At p=0.05, CDF = 1.0 (100% of values ≤ 0.05)

### How We Use CDF for Correction

When we observe a raw p-value (e.g., 0.01), we ask:
> "In the null distribution for this stratum, what fraction of p-values are ≤ 0.01?"

The CDF gives us the answer. This fraction becomes our **corrected p-value**.

### Example Calculation

**Stratum:** (GO BP, Iteration 2, Term 50-100, Overlap 4)

**Null distribution** (from 8,923 permutations):
- 10% of p-values are ≤ 0.005
- 30% of p-values are ≤ 0.01  ← Your observed p-value
- 60% of p-values are ≤ 0.02
- 100% of p-values are ≤ 0.05

**Your observed p-value:** 0.01

**Correction:**
- Look up in CDF: 30% of null p-values are ≤ 0.01
- **Corrected p-value = 0.30**

This means your p-value of 0.01 is actually at the 30th percentile of the null distribution, not as significant as it first appeared!

### Why This Works

Under the null hypothesis, p-values should be **uniformly distributed** (all values equally likely). But in iGEA, they're not uniform because of dependencies.

By using the CDF to transform p-values, we're saying:
> "Instead of comparing to a uniform distribution, compare to the actual null distribution for this specific condition."

After correction, the p-values **become uniform** under the null, which is what we need for proper statistical testing.

---

## Putting It Together

### Step-by-Step Process

1. **Create Strata**: Group permutation results by (Library, Iteration, Term size, Overlap)

2. **Build CDFs**: For each stratum with enough data, create a CDF from the p-values

3. **Correct P-values**: When you see a new p-value:
   - Find its stratum
   - Look up the p-value in that stratum's CDF
   - The CDF value is the corrected p-value

### Example

**Input:**
- Raw p-value: 0.01
- Library: GO BP
- Iteration: 2
- Term size: 75
- Overlap: 4

**Process:**
1. Identify stratum: (GO BP, Iteration 2, Term 50-100, Overlap 4)
2. Find CDF for this stratum
3. Query CDF: "What fraction of null p-values are ≤ 0.01?"
4. Answer: 0.30 (30%)

**Output:**
- Corrected p-value: 0.30

**Interpretation:**
- The raw p-value of 0.01 seemed significant
- But in this stratum, 30% of null p-values are ≤ 0.01
- So it's not as significant as it appeared
- The corrected p-value of 0.30 reflects this

---

## Key Takeaways

1. **Stratum** = A group of similar conditions (Library × Iteration × Term size × Overlap)

2. **CDF** = A function that tells you what fraction of values are below a threshold

3. **Correction** = Transform raw p-value using the CDF to account for dependencies

4. **Result** = P-values that are uniform under the null, ready for proper statistical testing

---

## Technical Note

In the code, we use `scipy.stats.rv_histogram` to create a CDF from a histogram:

```python
# Create histogram of p-values
counts, bins = np.histogram(p_values, bins=100, range=(0, 0.05))

# Normalize to get probabilities
counts = counts / counts.sum()

# Create CDF object
cdf = rv_histogram((counts, bins))

# Use it to get corrected p-value
corrected = cdf.cdf(0.01)  # Returns fraction ≤ 0.01
```

This gives us a smooth, continuous CDF even though we started with discrete p-values.
