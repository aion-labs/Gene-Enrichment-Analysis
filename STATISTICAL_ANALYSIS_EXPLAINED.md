# Statistical Analysis in iGEA: A Comprehensive Guide

## Table of Contents

1. [Introduction: Why Do We Need Statistics?](#introduction)
2. [The Core Question](#core-question)
3. [Permutation Testing: The Big Idea](#permutation-testing)
4. [Step 1: Generating Random Data (Permutations)](#generating-permutations)
5. [Step 2: Analyzing Random Networks](#analyzing-random-networks)
6. [Step 3: Building the Reference Distribution](#building-reference)
7. [Step 4: Runtime Analysis - Comparing Your Data](#runtime-analysis)
8. [Understanding the Results](#understanding-results)
9. [Key Concepts Summary](#summary)

---

## Introduction: Why Do We Need Statistics? {#introduction}

When you analyze a list of genes with iGEA, you get results showing which biological pathways, processes, or functions are enriched. But here's the critical question: **Are these results meaningful, or could they have happened by chance?**

### A Simple Analogy

Imagine you're a teacher and you notice that students who sit in the front row tend to get better grades. Is this meaningful, or just random chance? To find out, you could:

1. **Randomly assign students to seats** many times (say, 1,000 times)
2. **Count how often** the front-row students get better grades just by chance
3. **Compare** your actual observation to these random results

If your actual observation is much better than what happened in 999 out of 1,000 random assignments, you can be confident it's not just chance.

**This is exactly what iGEA does**, but instead of students and grades, we're looking at:
- **Genes** and **biological pathways**
- **Real gene lists** vs. **random gene lists**
- **Network connectivity** (how well genes and pathways connect together)

---

## The Core Question {#core-question}

When you analyze a gene list with iGEA, you want to know:

> **"Is the network of connections between my genes and enriched pathways better than what we'd expect from a random list of genes?"**

### Why This Matters

- **Real biological gene lists** (e.g., genes involved in a disease) should form **well-connected networks** because they work together in biological processes
- **Random gene lists** should form **loosely connected, fragmented networks** because there's no biological reason for them to connect

If your gene list forms a better-connected network than random lists, it suggests your genes are **biologically meaningful** and work together in real biological processes.

---

## Permutation Testing: The Big Idea {#permutation-testing}

### What is Permutation Testing?

**Permutation testing** is a statistical method that answers: "Could this result have happened by chance?"

The process:
1. **Generate many random versions** of your data (permutations)
2. **Analyze each random version** the same way you analyzed your real data
3. **Build a distribution** of what random results look like
4. **Compare your real result** to this distribution

### The iGEA Approach

In iGEA, we use permutation testing to answer:

> **"Is the network connectivity of my gene list better than random gene lists of the same size?"**

**Network connectivity** means: How well do genes and pathways connect together? Do they form large, dense clusters, or are they fragmented?

---

## Step 1: Generating Random Data (Permutations) {#generating-permutations}

### What Are Permutations?

**Permutations** are random gene lists created by randomly selecting genes from the background (all possible genes in the organism).

### The Permutation Generation Process

#### 1. **Define the Background**
- Start with all known genes in the organism (e.g., ~20,000 human genes)
- This is your "pool" of possible genes

#### 2. **Create Random Gene Lists**
For each permutation:
- **Randomly select** a specific number of genes (e.g., 50, 100, 150, 200, etc.)
- **No biological meaning** - just random selection
- **Same size** as the gene list sizes you want to test

#### 3. **Run iGEA on Each Random List**
- Process each random gene list through iGEA **exactly like a real gene list**
- Use the same libraries (pathway databases)
- Use the same parameters (p-value threshold, etc.)
- Get enrichment results

#### 4. **Store the Results**
- Save enrichment results for each permutation
- Each permutation file contains: which pathways were enriched, which genes were in each pathway

### Example: Generating 1,000 Permutations

```
Permutation 1: Random 200 genes → iGEA → Results saved
Permutation 2: Random 200 genes → iGEA → Results saved
Permutation 3: Random 200 genes → iGEA → Results saved
...
Permutation 1,000: Random 200 genes → iGEA → Results saved
```

### Why Multiple Sizes?

We generate permutations for **different gene list sizes** (50, 100, 150, 200, 250, etc.) because:
- Network connectivity depends on the number of genes
- A list of 50 genes will naturally have different connectivity than 200 genes
- We need reference data for each size

### Current Setup

- **~20,000 permutations** total
- **Multiple gene list sizes**: 50, 100, 150, 200, 250, 300, 400, 500, 750, 1000
- **Multiple libraries**: Reactome, KEGG, GO Biological Process, GO Molecular Function, GO Cellular Component
- **Strict p-value threshold**: 0.01 (only highly significant enrichments)

---

## Step 2: Analyzing Random Networks {#analyzing-random-networks}

### What is a Network?

When iGEA finds enriched pathways, it creates a **network**:
- **Nodes**: Genes and pathways (terms)
- **Edges**: Connections between genes and pathways (a gene is "in" a pathway)
- **Clusters**: Groups of pathways connected through shared genes

### Visual Example

```
Real Gene List Network:
    Gene1 ──┐
    Gene2 ─┼── Pathway A
    Gene3 ──┘
    
    Gene2 ──┐
    Gene4 ─┼── Pathway B
    Gene5 ──┘
    
    Gene3 ──┐
    Gene5 ─┼── Pathway C
    Gene6 ──┘

This forms ONE LARGE CLUSTER because pathways share genes.
```

```
Random Gene List Network:
    GeneX ── Pathway D
    
    GeneY ── Pathway E
    
    GeneZ ── Pathway F

This forms THREE SMALL CLUSTERS (no connections).
```

### How We Analyze Random Networks

For each permutation result, we compute **network metrics**:

#### 1. **Cluster Metrics** (Most Important)

**What is a cluster?**
- A cluster is a group of pathways connected through shared genes
- If Pathway A and Pathway B both contain Gene X, they're in the same cluster

**Metrics we measure:**
- **Largest cluster size**: How many genes + pathways in the biggest cluster?
- **Largest cluster genes**: How many genes in the biggest cluster?
- **Largest cluster terms**: How many pathways in the biggest cluster?
- **Largest cluster edges**: How many gene-pathway connections in the biggest cluster?
- **Largest cluster density**: How "tightly connected" is the biggest cluster?

**Why density matters:**
- Density = (actual connections) / (maximum possible connections)
- In iGEA, each gene can connect to at most 1 pathway per library
- So maximum connections = (number of genes) × (number of libraries)
- Higher density = more connections = better connectivity

#### 2. **Network-Wide Metrics**

- **Number of clusters**: How many separate groups?
  - Fewer clusters = better (more connected)
  - More clusters = worse (more fragmented)
- **Average cluster size**: Average size across all clusters
- **Fraction in largest cluster**: What percentage of everything is in the biggest cluster?

### The Analysis Process

For each permutation file:
1. **Load the enrichment results**
2. **Build the network** (genes ↔ pathways)
3. **Find clusters** (connected components)
4. **Compute metrics** for each cluster
5. **Store the metrics**

### Example: Analyzing One Permutation

```
Permutation File: random_200_genes_001.json

Enrichment Results:
- Pathway A: contains genes [Gene1, Gene2, Gene3]
- Pathway B: contains genes [Gene2, Gene4, Gene5]
- Pathway C: contains genes [Gene3, Gene5, Gene6]

Network Analysis:
- 6 genes, 3 pathways
- 1 cluster (all pathways connected through shared genes)
- Largest cluster: 6 genes, 3 pathways, 9 edges
- Density: 9 / (6 × 3) = 0.50

Metrics Saved:
- largest_cluster_genes: 6
- largest_cluster_terms: 3
- largest_cluster_edges: 9
- largest_cluster_density: 0.50
- n_connected_components: 1
```

---

## Step 3: Building the Reference Distribution {#building-reference}

### What is a Reference Distribution?

After analyzing all permutations, we have **thousands of measurements** for each metric. We organize these into a **reference distribution** (also called a "null distribution").

### The Distribution

Think of it like this: If you measured the height of 1,000 random people, you'd get:
- **Average height**: 170 cm
- **Standard deviation**: 10 cm
- **Range**: 150-190 cm

Similarly, for random gene lists of size 200:
- **Average largest cluster size**: 89 genes + pathways
- **Standard deviation**: 29
- **Range**: 30-150

### How We Build It

For each gene list size (50, 100, 200, etc.):

1. **Collect all metrics** from all permutations of that size
2. **Compute statistics**:
   - **Mean** (average)
   - **Standard deviation** (spread)
   - **Minimum/Maximum**
   - **Median** (middle value)

3. **Store in Parquet files** (fast, efficient format)
   - One file per gene list size
   - Organized by library combinations
   - Easy to query and filter

### Example: Reference Distribution for Size 200

```
Gene List Size: 200
Number of Permutations: 1,000

Metric: largest_cluster_genes
  Mean: 64.7 genes
  Std Dev: 20.8
  Min: 10
  Max: 114
  Median: 67

Metric: largest_cluster_terms
  Mean: 24.8 pathways
  Std Dev: 9.0
  Min: 2
  Max: 45
  Median: 25

Metric: largest_cluster_density
  Mean: 0.35
  Std Dev: 0.04
  Min: 0.15
  Max: 0.55
  Median: 0.36
```

### Why This Matters

This reference distribution tells us: **"What does a random gene list of size 200 typically look like?"**

When we analyze your real gene list, we can compare it to this reference.

---

## Step 4: Runtime Analysis - Comparing Your Data {#runtime-analysis}

### What Happens When You Run iGEA?

When you submit a gene list for analysis, here's what happens:

#### Phase 1: Enrichment (Main Thread)

1. **Load your gene list**
2. **Run iGEA enrichment** for each library
   - Find enriched pathways
   - Iteratively refine the gene set
   - Generate results

#### Phase 2: Statistical Benchmarking (Parallel Thread)

**While enrichment is running**, a separate process:

1. **Checks which libraries have permutation data**
   - Not all libraries may have been included in permutation generation
   - Only libraries with data can be used for statistics

2. **Loads the appropriate reference distribution**
   - Finds the Parquet file for your gene list size (or nearest size)
   - Filters to only include libraries you selected that have permutation data

3. **Computes null distribution statistics**
   - Loads cluster statistics from Parquet
   - Filters by your selected libraries
   - Computes mean, std dev, etc. for each metric

#### Phase 3: Network Analysis

After enrichment completes:

1. **Build your network**
   - Create the gene ↔ pathway network from your results
   - Find clusters (connected components)

2. **Compute your metrics**
   - Same metrics as for random data:
     - Largest cluster size, genes, terms, edges, density
     - Number of clusters
     - Average cluster size
     - etc.

#### Phase 4: Statistical Comparison

Compare your metrics to the reference distribution:

For each metric:
1. **Get your value** (e.g., largest_cluster_genes = 174)
2. **Get reference statistics** (mean = 64.7, std = 20.8)
3. **Compute Z-score**: 
   ```
   Z = (your_value - mean) / std_dev
   Z = (174 - 64.7) / 20.8 = 5.25
   ```
4. **Compute percentile**: What percentage of random lists had values lower than yours?
   - Z = 5.25 means you're **much better** than random
   - Percentile ≈ 100% (top 0.001%)

### Example: Full Runtime Process

```
Your Gene List: 200 genes related to HIV infection

[Parallel Thread] Loading reference distribution...
  ✓ Found Parquet file for size 200
  ✓ Filtered to your libraries: Reactome, KEGG, GO BP, GO MF, GO CC
  ✓ Computed statistics in 1.5 seconds

[Main Thread] Running iGEA enrichment...
  ✓ Reactome: 15 enriched pathways
  ✓ KEGG: 8 enriched pathways
  ✓ GO BP: 22 enriched pathways
  ✓ GO MF: 12 enriched pathways
  ✓ GO CC: 9 enriched pathways

[After Enrichment] Building network...
  ✓ 200 genes, 66 pathways
  ✓ 1 large cluster (174 genes, 66 pathways, 330 edges)
  ✓ Density: 0.38

[Statistical Comparison]
  Your largest_cluster_genes: 174
  Random average: 64.7
  Z-score: 5.25
  Percentile: 100%
  Status: ✓✓ SIGNIFICANTLY BETTER
```

---

## Understanding the Results {#understanding-results}

### What Do the Numbers Mean?

When you get benchmark results, you'll see metrics like:

#### Z-Score

**What it is**: How many standard deviations above (or below) the random average you are.

**Interpretation**:
- **Z > 2.0**: Much better than random (top ~2.5%)
- **Z > 1.0**: Better than random (top ~16%)
- **Z ≈ 0**: Similar to random (middle 68%)
- **Z < -1.0**: Worse than random (bottom ~16%)
- **Z < -2.0**: Much worse than random (bottom ~2.5%)

**Example**:
- Your Z-score = 5.25
- This means you're **5.25 standard deviations** above the random average
- This is **extremely rare** in random data (top 0.0001%)
- **Conclusion**: Your result is highly significant

#### Percentile

**What it is**: What percentage of random gene lists had values lower than yours.

**Interpretation**:
- **Percentile > 95%**: Top 5% (very good)
- **Percentile > 84%**: Top 16% (good)
- **Percentile ≈ 50%**: Average (similar to random)
- **Percentile < 16%**: Bottom 16% (worse than random)
- **Percentile < 5%**: Bottom 5% (much worse than random)

**Example**:
- Your percentile = 100%
- This means **all** random lists had lower values than yours
- **Conclusion**: Your result is exceptional

#### Status Indicators

- **✓✓ SIGNIFICANTLY BETTER**: Z > 2.0 and Percentile > 95%
  - Your gene list forms a much better-connected network than random
  - Strong evidence of biological coherence
  
- **✓ BETTER**: Z > 1.0 and Percentile > 84%
  - Your gene list forms a better-connected network than random
  - Evidence of biological coherence
  
- **~ SIMILAR**: -1.0 < Z < 1.0 and 16% < Percentile < 84%
  - Your gene list is similar to random
  - May indicate weak biological coherence or need for more genes
  
- **✗ WORSE**: Z < -1.0 and Percentile < 16%
  - Your gene list forms a less-connected network than random
  - Unusual - may indicate data quality issues

### Key Metrics Explained

#### Largest Cluster Genes

**What it measures**: How many genes participate in the main connected component.

**Why it matters**: 
- Real biological gene lists should have many genes working together
- Random lists typically have fewer genes in the largest cluster

**Good value**: Much higher than random average

#### Largest Cluster Terms

**What it measures**: How many pathways are connected together in the main cluster.

**Why it matters**:
- Real biological processes involve multiple related pathways
- Random lists typically have fewer pathways connected

**Good value**: Much higher than random average

#### Largest Cluster Density

**What it measures**: How tightly connected the main cluster is (connections / maximum possible).

**Why it matters**:
- Real biological networks are densely connected
- Random networks are more sparse

**Good value**: Higher than random average (typically > 0.35)

#### Number of Clusters

**What it measures**: How many separate groups of pathways exist.

**Why it matters**:
- Real biological gene lists should form fewer, larger clusters
- Random lists form many small, disconnected clusters

**Good value**: **Lower** than random average (fewer is better)

### Example: Interpreting Results

```
Your Gene List: 200 genes, HIV-related

Metric: largest_cluster_genes
  Your value: 174 genes
  Random average: 64.7 genes
  Z-score: 5.25
  Percentile: 100%
  Status: ✓✓ SIGNIFICANTLY BETTER

Interpretation:
- Your 200 genes form a network where 174 genes are in one large cluster
- Random 200-gene lists typically have only 65 genes in the largest cluster
- Your result is in the top 0.0001% of random lists
- This strongly suggests your genes work together biologically
```

---

## Key Concepts Summary {#summary}

### 1. **Permutation Testing**
- Generate many random versions of your data
- Analyze them the same way
- Compare your real result to the random results

### 2. **Network Connectivity**
- Genes and pathways form networks
- Real biological gene lists form well-connected networks
- Random gene lists form fragmented networks

### 3. **Clusters**
- Groups of pathways connected through shared genes
- Larger, denser clusters = better connectivity
- Real lists should have fewer, larger clusters than random lists

### 4. **Reference Distribution**
- Built from thousands of random gene list analyses
- Tells us what "typical" random results look like
- Organized by gene list size and library combinations

### 5. **Statistical Comparison**
- Compare your metrics to the reference distribution
- Z-score: How many standard deviations above/below average
- Percentile: What percentage of random lists were worse
- Status: Simple indicator (SIGNIFICANTLY BETTER, BETTER, SIMILAR, WORSE)

### 6. **Runtime Process**
- Enrichment runs in main thread
- Statistical benchmarking runs in parallel thread
- Only libraries with permutation data are used for statistics
- Full network results include all selected libraries

### 7. **Interpretation**
- **SIGNIFICANTLY BETTER**: Strong evidence of biological coherence
- **BETTER**: Evidence of biological coherence
- **SIMILAR**: May indicate weak coherence or need for more genes
- **WORSE**: Unusual - may indicate data quality issues

---

## Frequently Asked Questions

### Q: Why do we need so many permutations?

**A**: More permutations = more accurate reference distribution. With 1,000 permutations, we can reliably estimate what random results look like. With only 10 permutations, the estimate would be unreliable.

### Q: What if my gene list size doesn't match exactly?

**A**: The system uses the nearest available size or interpolates between two sizes. For example, if you have 175 genes and we have data for 150 and 200, we interpolate.

### Q: Why are some libraries excluded from statistics?

**A**: Only libraries that were included in permutation generation can be used for statistics. If a library wasn't included in permutations, we can't compare against random data for that library. However, you still get full enrichment results for all libraries.

### Q: What does "density" mean?

**A**: Density measures how "tightly connected" a cluster is. It's the ratio of actual connections to maximum possible connections. Higher density = more connections = better connectivity.

### Q: Can I trust results that are "SIMILAR" to random?

**A**: "SIMILAR" doesn't mean your genes aren't biologically related. It might mean:
- Your gene list is too small to show strong connectivity
- Your genes are related but through pathways not in the libraries used
- The connectivity pattern is different but not necessarily worse

### Q: How long does statistical benchmarking take?

**A**: Very fast! The parallel computation typically takes 1-2 seconds, running simultaneously with enrichment (which takes much longer). There's essentially no time cost.

---

## Conclusion

The statistical analysis in iGEA provides a rigorous way to assess whether your gene list results are biologically meaningful or could have occurred by chance. By comparing your network connectivity to thousands of random gene lists, we can quantify how exceptional your results are and provide confidence in their biological significance.

The system is designed to be:
- **Fast**: Parallel computation with no performance penalty
- **Accurate**: Based on thousands of permutations
- **Flexible**: Works with different gene list sizes and library combinations
- **Transparent**: Clear reporting of which libraries were used and how results compare to random

If you have questions about interpreting your specific results, refer to the status indicators and percentiles - they provide clear, quantitative measures of how your gene list compares to random expectations.
