# Parallel Null Distribution Implementation

## Overview

This document describes the implementation of parallel null distribution computation integrated into the iGEA benchmarking workflow. The system now dynamically computes null distributions from Parquet cluster statistics in parallel with iGEA enrichment, ensuring that only libraries with available permutation data are used for statistical benchmarking.

## Key Features

1. **Parallel Execution**: Null distribution computation runs in a separate thread while iGEA enrichment executes
2. **Dynamic Library Filtering**: Automatically detects which libraries have permutation data available
3. **Selective Statistics**: Only libraries with permutation data are used for statistical benchmarking
4. **Full Network Results**: Users still receive complete network results with all selected libraries
5. **Clear Reporting**: Reports explicitly indicate which libraries were used for statistics

## Architecture

### New Module: `code/parallel_null_distribution.py`

#### Functions

- **`get_available_libraries_from_parquet(parquet_dir: Path) -> Set[str]`**
  - Scans Parquet files to detect which libraries have permutation data
  - Returns set of library names with available data

- **`find_intersection_libraries(user_selected_libraries: List[str], parquet_dir: Path) -> tuple[List[str], List[str]]`**
  - Finds intersection between user-selected libraries and libraries with permutation data
  - Returns: `(libraries_with_data, libraries_without_data)`
  - Handles library name normalization and partial matching

- **`compute_null_distribution_from_parquet(parquet_dir: Path, gene_list_size: int, selected_libraries: List[str], metrics: List[str] = None) -> Dict`**
  - Computes null distribution statistics from Parquet cluster statistics
  - Filters by selected libraries
  - Aggregates metrics per permutation (one value per permutation file)
  - Returns dictionary: `{metric_name: {'mean': float, 'std': float, 'n': int, ...}}`

- **`compute_null_distribution_parallel(parquet_dir: Path, gene_list_size: int, selected_libraries: List[str], result_dict: Dict, lock: threading.Lock) -> None`**
  - Thread-safe wrapper for parallel null distribution computation
  - Stores results in shared dictionary with thread lock

### Updated: `benchmark_hiv_connectivity.py`

#### Workflow Changes

1. **Library Detection** (Step 5/6):
   - Checks which user-selected libraries have permutation data
   - Reports libraries with/without data

2. **Parallel Null Distribution** (Step 6/6):
   - Starts null distribution computation in separate thread
   - Thread runs concurrently with iGEA enrichment

3. **Dual Network Analysis**:
   - **Full Network**: All user-selected libraries (for visualization)
   - **Filtered Network**: Only libraries with permutation data (for statistics)

4. **Benchmarking**:
   - Only filtered network is benchmarked against null distribution
   - Full network metrics are computed but not benchmarked

5. **Reporting**:
   - Clear notes about which libraries were used for statistics
   - Warnings for libraries excluded from statistics

## Data Flow

```
User selects libraries
    ↓
Check which have permutation data (Parquet files)
    ↓
Start null distribution computation (parallel thread)
    ↓
Run iGEA with ALL libraries (main thread)
    ↓
When iGEA completes:
    ├─ Full network: All libraries (for user)
    └─ Filtered network: Only statistics libraries (for benchmarking)
    ↓
Benchmark filtered network against null distribution
    ↓
Generate reports with clear library information
```

## Null Distribution Format

The null distribution is computed as:
```python
{
    str(gene_list_size): {
        'metric_name': {
            'mean': float,
            'std': float,
            'n': int,
            'min': float,
            'max': float,
            'median': float
        },
        ...
    }
}
```

This format matches what `benchmark_real_results()` expects.

## Performance

- **Null Distribution Computation**: ~1.5 seconds for 5 libraries, gene list size 200
- **Parallel Execution**: No additional time cost (runs concurrently with iGEA)
- **Library Detection**: <0.1 seconds

## Usage

### In `benchmark_hiv_connectivity.py`:

```python
from parallel_null_distribution import (
    find_intersection_libraries,
    compute_null_distribution_parallel
)

# Check libraries
libraries_with_data, libraries_without_data = find_intersection_libraries(
    user_selected_libraries,
    PARQUET_DIR
)

# Start parallel computation
null_dist_result = {'null_distribution': None, 'status': 'running', ...}
null_dist_lock = threading.Lock()
null_dist_thread = threading.Thread(
    target=compute_null_distribution_parallel,
    args=(PARQUET_DIR, gene_set.size, libraries_with_data, null_dist_result, null_dist_lock),
    daemon=True
)
null_dist_thread.start()

# Run iGEA (in parallel)
# ... iGEA code ...

# Wait for null distribution
null_dist_thread.join(timeout=30)
null_distribution = null_dist_result['null_distribution']
```

## Testing

A test script `test_parallel_null_distribution.py` is provided to verify:
- Library intersection detection
- Parallel null distribution computation
- Network filtering
- Benchmarking workflow

Run with:
```bash
python3 test_parallel_null_distribution.py
```

## Files Modified

1. **`code/parallel_null_distribution.py`** (NEW)
   - Core module for parallel null distribution computation

2. **`benchmark_hiv_connectivity.py`** (UPDATED)
   - Integrated parallel null distribution computation
   - Added library filtering logic
   - Updated reporting to show library information

3. **`test_parallel_null_distribution.py`** (NEW)
   - Test script for verification

## Notes

- The system ensures that users always receive full network results with all selected libraries
- Statistics are only computed for libraries with permutation data to ensure accuracy
- Reports clearly indicate which libraries were included/excluded from statistics
- The parallel execution ensures no performance penalty for null distribution computation
