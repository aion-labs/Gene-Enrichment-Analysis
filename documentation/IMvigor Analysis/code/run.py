#!/usr/bin/env python3
"""
IMvigor210 Biomarker Analysis Pipeline
Main orchestrator — runs all modules in sequence.
"""

import os
import sys
import time

# Ensure code directory is on path
sys.path.insert(0, '/code')

from module1_data import run_module1
from module2_clinical import run_module2
from module3_genomic import run_module3
from module4_expression import run_module4
from module5_signature import run_module5
from module6_integrated import run_module6
from module7_report import run_module7

RESULTS_DIR = '/results'


def main():
    start = time.time()

    os.makedirs(RESULTS_DIR, exist_ok=True)

    print("=" * 60)
    print("IMvigor210 BIOMARKER ANALYSIS PIPELINE")
    print("=" * 60)

    # Module 1: Data loading, QC, preprocessing
    data = run_module1()

    # Module 2: Clinical feature analysis
    m2_results = run_module2(data)

    # Module 3: Genomic alteration analysis
    m3_results = run_module3(data)

    # Module 4: Differential expression & pathway analysis
    m4_results = run_module4(data)

    # Module 5: Predictive gene expression signature (depends on M4)
    m5_results = run_module5(data, m4_results)

    # Module 6: Integrated model & summary (depends on M2-M5)
    m6_results = run_module6(data, m2_results, m3_results, m4_results, m5_results)

    # Module 7: PDF report
    run_module7()

    elapsed = time.time() - start
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print(f"Elapsed: {elapsed:.1f} seconds")
    print(f"Results saved to: {RESULTS_DIR}")
    print("=" * 60)

    # List outputs
    outputs = sorted(os.listdir(RESULTS_DIR))
    print(f"\nGenerated {len(outputs)} output files:")
    for f in outputs:
        size = os.path.getsize(os.path.join(RESULTS_DIR, f))
        print(f"  {f} ({size:,} bytes)")


if __name__ == '__main__':
    main()
