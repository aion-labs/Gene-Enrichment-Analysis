#!/usr/bin/env python3
"""
Check progress of running permutation distribution job.

This script analyzes the output directory to determine:
- How many permutations have been completed
- Estimated time remaining
- Whether the job is progressing normally
"""

import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, Set

# Detect Code Ocean environment
if Path("/results").exists() and Path("/results").is_dir():
    RESULTS_DIR = Path("/results")
else:
    RESULTS_DIR = Path(__file__).resolve().parent / "results"

OUTPUT_DIR = RESULTS_DIR / "permutation_results"


def get_completed_permutations(output_dir: Path, size: int) -> Set[int]:
    """Get set of completed permutation indices for a given size."""
    completed = set()
    size_dir = output_dir / f"size_{size}"
    
    if not size_dir.exists():
        return completed
    
    for file in size_dir.glob("permutation_*.tsv"):
        try:
            # Extract permutation index from filename
            # Format: permutation_0001.tsv -> 1
            idx_str = file.stem.replace("permutation_", "")
            perm_idx = int(idx_str)
            completed.add(perm_idx)
        except (ValueError, AttributeError):
            continue
    
    return completed


def check_progress():
    """Check and report progress of permutation job."""
    if not OUTPUT_DIR.exists():
        print(f"❌ Output directory not found: {OUTPUT_DIR}")
        print("   The job may not have started yet.")
        return
    
    # Load configuration
    config_file = OUTPUT_DIR / "config.json"
    if not config_file.exists():
        print(f"❌ Configuration file not found: {config_file}")
        print("   Cannot determine expected parameters.")
        return
    
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    n_permutations = config.get("n_permutations", 1000)
    sizes = config.get("gene_list_sizes", [50])
    timestamp = config.get("timestamp", "")
    
    print("=" * 70)
    print("PERMUTATION JOB PROGRESS REPORT")
    print("=" * 70)
    print(f"Configuration timestamp: {timestamp}")
    print(f"Permutations per size: {n_permutations:,}")
    print(f"Gene list sizes: {sizes}")
    print(f"Total sizes to process: {len(sizes)}")
    print("-" * 70)
    
    # Check progress for each size
    total_completed = 0
    total_expected = len(sizes) * n_permutations
    size_progress = {}
    
    for size in sizes:
        completed = get_completed_permutations(OUTPUT_DIR, size)
        completed_count = len(completed)
        remaining = n_permutations - completed_count
        pct = (completed_count / n_permutations * 100) if n_permutations > 0 else 0
        
        size_progress[size] = {
            "completed": completed_count,
            "remaining": remaining,
            "percent": pct
        }
        total_completed += completed_count
        
        status = "✅" if completed_count == n_permutations else "⏳" if completed_count > 0 else "⏸️"
        print(f"{status} Size {size:4d}: {completed_count:4d}/{n_permutations:4d} ({pct:5.1f}%) - {remaining:4d} remaining")
    
    print("-" * 70)
    overall_pct = (total_completed / total_expected * 100) if total_expected > 0 else 0
    print(f"OVERALL: {total_completed:,}/{total_expected:,} ({overall_pct:.1f}%) completed")
    print(f"         {total_expected - total_completed:,} remaining")
    print("=" * 70)
    
    # Estimate time remaining (very rough estimate)
    # Based on AWS recommendations: ~0.05-0.1 seconds per permutation
    # This is a conservative estimate
    avg_time_per_perm = 0.1  # seconds
    remaining_permutations = total_expected - total_completed
    estimated_seconds = remaining_permutations * avg_time_per_perm
    
    if estimated_seconds > 0:
        hours = estimated_seconds / 3600
        print(f"\n⏱️  ROUGH ESTIMATE (assuming 0.1s/permutation):")
        print(f"   Remaining time: ~{hours:.1f} hours")
        print(f"   (This is a rough estimate - actual time depends on many factors)")
    
    # Check for recent activity
    print("\n📊 RECENT ACTIVITY:")
    recent_files = []
    for size in sizes:
        size_dir = OUTPUT_DIR / f"size_{size}"
        if size_dir.exists():
            for file in size_dir.glob("permutation_*.tsv"):
                try:
                    mtime = file.stat().st_mtime
                    recent_files.append((mtime, file))
                except:
                    pass
    
    if recent_files:
        recent_files.sort(reverse=True)
        latest_file, latest_mtime = recent_files[0][1], recent_files[0][0]
        latest_time = datetime.fromtimestamp(latest_mtime)
        time_since = datetime.now().timestamp() - latest_mtime
        
        print(f"   Latest file: {latest_file.name}")
        print(f"   Last modified: {latest_time.strftime('%Y-%m-%d %H:%M:%S')}")
        
        if time_since < 3600:
            print(f"   ✅ Active in last hour ({int(time_since/60)} minutes ago)")
        elif time_since < 86400:
            print(f"   ⚠️  Last activity {int(time_since/3600)} hours ago")
        else:
            print(f"   ❌ No activity in last {int(time_since/86400)} days - job may be stuck!")
    else:
        print("   ⚠️  No output files found")
    
    # Recommendations
    print("\n💡 RECOMMENDATIONS:")
    if overall_pct < 10:
        print("   - Job is still in early stages")
        print("   - Consider checking if parallelization is working (n_jobs parameter)")
    elif overall_pct < 50:
        print("   - Job is progressing but still has significant work remaining")
        print("   - Monitor for continued progress")
    elif overall_pct < 90:
        print("   - Job is more than halfway complete")
        print("   - Should complete in reasonable time if progress continues")
    else:
        print("   - Job is nearly complete!")
        print("   - Should finish soon")
    
    if time_since > 3600 and overall_pct < 100:
        print("\n⚠️  WARNING: No recent activity detected!")
        print("   - Job may be stuck or waiting for resources")
        print("   - Check Code Ocean logs for errors")
        print("   - Consider canceling and restarting with fewer permutations for testing")


if __name__ == "__main__":
    check_progress()




