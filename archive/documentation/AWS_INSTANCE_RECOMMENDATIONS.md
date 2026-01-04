# AWS Instance Recommendations for Permutation Distribution Script

## Workload Characteristics

- **CPU-intensive**: Statistical calculations (Fisher's exact test) for 20,000 permutations
- **Parallel processing**: Script uses multiprocessing (configurable with `--n-jobs`)
- **Memory**: Moderate (2-4 GB estimated, libraries loaded once)
- **Storage I/O**: High (writing 20,000 TSV files)
- **Network**: Low (no external data transfer needed)
- **Expected runtime**: 
  - Single-threaded: ~110 hours
  - 8 cores: ~14 hours
  - 16 cores: ~7 hours
  - 32 cores: ~3.5 hours

## Recommended Instance Types

### 🏆 **Best Overall: c6i.xlarge or c6i.2xlarge**

**c6i.xlarge** (Recommended for cost-effectiveness):
- **vCPUs**: 4 cores
- **Memory**: 8 GB
- **Network**: Up to 12.5 Gbps
- **EBS Bandwidth**: Up to 10 Gbps
- **Estimated Runtime**: ~28 hours (with 4 parallel jobs)
- **Cost**: ~$0.17/hour = ~$4.76 for full run
- **Best for**: Budget-conscious runs, good balance of speed and cost

**c6i.2xlarge** (Recommended for speed):
- **vCPUs**: 8 cores
- **Memory**: 16 GB
- **Network**: Up to 12.5 Gbps
- **EBS Bandwidth**: Up to 10 Gbps
- **Estimated Runtime**: ~14 hours (with 8 parallel jobs)
- **Cost**: ~$0.34/hour = ~$4.76 for full run
- **Best for**: Faster completion, still cost-effective

### ⚡ **Fastest: c6i.4xlarge or c6i.8xlarge**

**c6i.4xlarge**:
- **vCPUs**: 16 cores
- **Memory**: 32 GB
- **Network**: Up to 25 Gbps
- **EBS Bandwidth**: Up to 20 Gbps
- **Estimated Runtime**: ~7 hours (with 16 parallel jobs)
- **Cost**: ~$0.68/hour = ~$4.76 for full run
- **Best for**: Time-sensitive runs, same cost as 2xlarge but faster

**c6i.8xlarge**:
- **vCPUs**: 32 cores
- **Memory**: 64 GB
- **Network**: Up to 50 Gbps
- **EBS Bandwidth**: Up to 40 Gbps
- **Estimated Runtime**: ~3.5 hours (with 32 parallel jobs)
- **Cost**: ~$1.36/hour = ~$4.76 for full run
- **Best for**: Maximum speed, same total cost

### 💰 **Budget Option: c6i.large**

**c6i.large**:
- **vCPUs**: 2 cores
- **Memory**: 4 GB
- **Network**: Up to 12.5 Gbps
- **EBS Bandwidth**: Up to 10 Gbps
- **Estimated Runtime**: ~55 hours (with 2 parallel jobs)
- **Cost**: ~$0.085/hour = ~$4.68 for full run
- **Best for**: Minimal cost, but very slow

## Why C6i (Intel) Instances?

1. **CPU-optimized**: Designed for compute-intensive workloads
2. **Latest generation**: Better price/performance than C5 instances
3. **Good EBS bandwidth**: Important for writing 20,000 TSV files
4. **Cost-effective**: Best price per vCPU for CPU-bound tasks

## Alternative: C6a (AMD) Instances

**c6a.xlarge** (AMD alternative):
- **vCPUs**: 4 cores
- **Memory**: 8 GB
- **Cost**: ~10-15% cheaper than c6i
- **Performance**: Similar to c6i for this workload
- **Best for**: If cost is primary concern

## Storage Recommendations

### EBS Volume Configuration:
- **Type**: gp3 (General Purpose SSD)
- **Size**: 50-100 GB (sufficient for 20,000 TSV files + libraries)
- **IOPS**: 3,000 (default, sufficient for TSV writes)
- **Throughput**: 125 MB/s (default, sufficient)

### Estimated Storage Needs:
- Libraries: ~500 MB - 2 GB (depends on library sizes)
- Output TSV files: ~100-500 MB (20,000 files, ~5-25 KB each)
- Logs: ~10-50 MB
- **Total**: ~2-3 GB minimum, 50 GB provides comfortable margin

## Code Ocean Configuration

### Recommended Settings:

```yaml
# For c6i.2xlarge (recommended)
instance_type: c6i.2xlarge
cpu_cores: 8
memory_gb: 16
storage_gb: 50
```

### Command to Run:

```bash
python generate_permutation_distribution.py --n-jobs 8
```

Adjust `--n-jobs` to match the number of vCPUs:
- c6i.large: `--n-jobs 2`
- c6i.xlarge: `--n-jobs 4`
- c6i.2xlarge: `--n-jobs 8`
- c6i.4xlarge: `--n-jobs 16`
- c6i.8xlarge: `--n-jobs 32`

## Cost Comparison

| Instance Type | vCPUs | Runtime | Cost/Hour | Total Cost |
|--------------|-------|---------|-----------|------------|
| c6i.large    | 2     | ~55h    | $0.085    | ~$4.68     |
| c6i.xlarge   | 4     | ~28h    | $0.17     | ~$4.76     |
| c6i.2xlarge  | 8     | ~14h    | $0.34     | ~$4.76     |
| c6i.4xlarge  | 16    | ~7h     | $0.68     | ~$4.76     |
| c6i.8xlarge  | 32    | ~3.5h   | $1.36     | ~$4.76     |

**Key Insight**: All options cost approximately the same (~$4.76) because:
- More cores = faster completion
- Cost scales linearly with cores
- Runtime scales inversely with cores
- **Total cost ≈ constant**

## Final Recommendation

### 🎯 **Primary Recommendation: c6i.2xlarge**
- **Why**: Best balance of speed (14 hours) and cost (~$4.76)
- **Use**: `--n-jobs 8`
- **Best for**: Most users, good compromise

### ⚡ **If Time-Critical: c6i.4xlarge**
- **Why**: Completes in ~7 hours (half the time)
- **Use**: `--n-jobs 16`
- **Best for**: When you need results quickly

### 💰 **If Budget-Conscious: c6i.xlarge**
- **Why**: Still reasonable speed (~28 hours) at same cost
- **Use**: `--n-jobs 4`
- **Best for**: When you can wait longer

## Important Notes

1. **Resumability**: The script can resume, so you can stop/start without losing progress
2. **Spot Instances**: Consider using spot instances for 50-70% cost savings (but risk of interruption)
3. **Monitoring**: Monitor CPU utilization - should be near 100% for all cores
4. **Storage**: Ensure sufficient EBS storage for output files
5. **Logs**: Check `permutation_distribution.log` for progress and errors

## Code Ocean Specific Tips

1. **Set timeout**: Ensure Code Ocean timeout is longer than expected runtime
2. **Output directory**: Results will be in `permutation_results/` directory
3. **Checkpoint**: Script automatically saves progress, can resume if interrupted
4. **Resource limits**: Ensure Code Ocean project has sufficient compute credits

## Troubleshooting

### If script is slow:
- Check CPU utilization (should be near 100%)
- Verify `--n-jobs` matches instance vCPUs
- Check EBS IOPS if disk writes are slow

### If out of memory:
- Reduce `--n-jobs` (fewer parallel processes)
- Use larger instance (more memory)

### If storage full:
- Increase EBS volume size
- Clean up old results if resuming

