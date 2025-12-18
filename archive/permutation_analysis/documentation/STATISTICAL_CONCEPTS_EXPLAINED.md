# Statistical Concepts Explained: Stratum and CDF

## What is a Stratum?

A **stratum** (plural: strata) is a subgroup or category that shares similar characteristics. In statistics, stratification means dividing your data into homogeneous groups based on specific criteria.

### In Our Context

For iGEA p-value correction, we create strata by grouping permutation results that have the **same combination** of:
- Library (e.g., "GO BP", "Reactome")
- Iteration bin (e.g., iteration 1, iterations 2-3, etc.)
- Term size bin (e.g., 50-100 genes)
- Overlap size bin (e.g., 4 overlapping genes)

### Example

Imagine you have these permutation results:

| Library | Iteration | Term Size | Overlap | P-value |
|---------|-----------|-----------|---------|---------|
| GO BP   | 1         | 150       | 4       | 0.012   |
| GO BP   | 1         | 150       | 4       | 0.008   |
| GO BP   | 1         | 150       | 4       | 0.015   |
| GO BP   | 2         | 150       | 4       | 0.020   |
| Reactome| 1         | 150       | 4       | 0.025   |

The first three rows belong to the **same stratum** because they share:
- Same library: "GO BP"
- Same iteration bin: 1
- Same term size bin: 100-200 (since 150 falls in this range)
- Same overlap bin: 4

The fourth row is in a **different stratum** (different iteration).
The fifth row is in yet another **different stratum** (different library).

### Why Stratify?

P-values behave differently depending on these factors. For example:
- P-values in iteration 1 tend to be smaller than in iteration 5
- P-values for GO BP terms tend to be smaller than for Reactome terms
- P-values for large terms (300+ genes) behave differently than small terms (<50 genes)

By creating separate strata, we can model the **null distribution** (what p-values look like under the null hypothesis) for each specific combination of conditions.

### Visual Example

```
All Permutation Results (1.25 million rows)
│
├── Stratum 1: (GO BP, Iteration 1, Term 50-100, Overlap 3)
│   └── 15,000 p-values → Build null distribution for this group
│
├── Stratum 2: (GO BP, Iteration 1, Term 50-100, Overlap 4)
│   └── 12,000 p-values → Build null distribution for this group
│
├── Stratum 3: (GO BP, Iteration 2, Term 50-100, Overlap 3)
│   └── 8,000 p-values → Build null distribution for this group
│
└── ... (527 total strata)
```

---

## What is a CDF?

**CDF** stands for **Cumulative Distribution Function**. It tells you the probability that a random variable is less than or equal to a given value.

### Simple Explanation

Think of it as answering the question: "What percentage of values are ≤ X?"

### Example with Heights

Imagine you measured the heights of 100 people:

| Height (cm) | Number of People | Percentage | CDF Value |
|-------------|------------------|------------|-----------|
| ≤ 150       | 5                | 5%         | 0.05      |
| ≤ 160       | 20               | 20%        | 0.20      |
| ≤ 170       | 50               | 50%        | 0.50      |
| ≤ 180       | 80               | 80%        | 0.80      |
| ≤ 190       | 95               | 95%        | 0.95      |
| ≤ 200       | 100              | 100%       | 1.00      |

The CDF at height 170 is 0.50, meaning 50% of people are ≤ 170 cm tall.

### Visual Representation

```
CDF
1.0 |                                    ┌─────
    |                                ┌───┘
0.8 |                            ┌───┘
    |                        ┌───┘
0.5 |                    ┌───┘
    |                ┌───┘
0.2 |            ┌───┘
    |        ┌───┘
0.0 |───────┘
    └─────────────────────────────────────
    150  160  170  180  190  200
         Height (cm)
```

### In Our Context: P-Value CDF

For each stratum, we build a CDF from the permutation p-values:

**Example Stratum: (GO BP, Iteration 1, Term 50-100, Overlap 4)**

We have 12,000 p-values from permutations. We create a histogram:

| P-value Range | Count | Percentage | CDF |
|---------------|-------|------------|-----|
| 0.000 - 0.001 | 200   | 1.7%       | 0.017 |
| 0.001 - 0.002 | 500   | 4.2%       | 0.059 |
| 0.002 - 0.005 | 2000  | 16.7%      | 0.226 |
| 0.005 - 0.010 | 3000  | 25.0%      | 0.476 |
| 0.010 - 0.020 | 4000  | 33.3%      | 0.809 |
| 0.020 - 0.050 | 2300  | 19.2%      | 1.000 |

**What this tells us:**
- 1.7% of null p-values are ≤ 0.001
- 5.9% of null p-values are ≤ 0.002
- 22.6% of null p-values are ≤ 0.005
- 47.6% of null p-values are ≤ 0.010
- 80.9% of null p-values are ≤ 0.020

### How We Use CDF for Correction

When we observe a p-value of 0.01 in this stratum:

1. **Look up the CDF**: The CDF at 0.01 is 0.476
2. **Interpret**: This means 47.6% of null p-values in this stratum are ≤ 0.01
3. **Correct**: We transform the p-value to 0.476 (the percentile rank)

**Why?** Under the null hypothesis, p-values should be uniformly distributed (each value equally likely). If we observe a p-value that's at the 47.6th percentile of the null distribution, we correct it to 0.476 to make it uniform.

### Visual Example of Correction

```
Raw P-value Distribution (skewed, not uniform)
│
│     ████
│   ████████
│  ████████████
│ ████████████████
└─────────────────────
 0.00             0.05
     (Many small p-values)
     
After CDF Transformation (uniform)
│
│ ████████████████████
│ ████████████████████
│ ████████████████████
│ ████████████████████
└─────────────────────
 0.00             1.00
     (Uniform distribution)
```

---

## Putting It Together: How Stratum + CDF Work

### Step-by-Step Example

**Scenario:** You observe a p-value of 0.008 in:
- Library: GO BP
- Iteration: 2
- Term size: 150 genes
- Overlap: 4 genes

**Step 1: Identify the Stratum**
- Library: GO BP
- Iteration bin: 2 (iterations 2-3)
- Term size bin: 100-200 (150 falls here)
- Overlap bin: 4

**Step 2: Find the CDF for This Stratum**
From permutation data, we built a CDF. Let's say:
- CDF(0.005) = 0.30 (30% of null p-values ≤ 0.005)
- CDF(0.008) = 0.45 (45% of null p-values ≤ 0.008)
- CDF(0.010) = 0.55 (55% of null p-values ≤ 0.010)

**Step 3: Correct the P-value**
- Raw p-value: 0.008
- CDF(0.008) = 0.45
- **Corrected p-value: 0.45**

**Interpretation:**
- The raw p-value of 0.008 is at the 45th percentile of the null distribution
- After correction, it becomes 0.45, which is uniform under the null
- This corrected p-value can now be used for multiple testing correction (like FDR)

### Why This Works

1. **Stratification** ensures we compare "apples to apples" - we only compare p-values from similar conditions
2. **CDF transformation** converts the skewed null distribution to uniform, making standard statistical tests valid
3. **Result**: Corrected p-values are now independent of library, iteration, term size, and overlap

---

## Summary

- **Stratum**: A subgroup of data with the same characteristics (library, iteration bin, term size bin, overlap bin)
- **CDF**: A function that tells you what percentage of values are ≤ a given value
- **Together**: We build separate CDFs for each stratum, then use them to transform p-values to be uniform under the null

This allows us to correct for the dependencies and make valid statistical inferences!
