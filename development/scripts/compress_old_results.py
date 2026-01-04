#!/usr/bin/env python3
"""
Compress results folders and files that are older than a specified age.

This script finds all folders and files in the results/ directory that are older than
the specified threshold (default: 10 hours) and compresses them:
- Folders: compressed using tar.gz
- Files: compressed using gzip (.gz)
The original items are removed after successful compression.

Usage:
    python3 compress_old_results.py [--hours 10] [--dry-run]
"""

import os
import subprocess
import shutil
import argparse
from pathlib import Path
from datetime import datetime, timedelta


def get_folder_size(folder_path):
    """Calculate total size of a folder in bytes."""
    total_size = 0
    try:
        for dirpath, dirnames, filenames in os.walk(folder_path):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                try:
                    total_size += os.path.getsize(fp)
                except (OSError, FileNotFoundError):
                    pass
    except (OSError, PermissionError):
        pass
    return total_size


def compress_old_items(results_dir, hours_threshold=10, dry_run=False):
    """
    Compress folders and files older than the specified threshold.
    
    Args:
        results_dir: Path to results directory
        hours_threshold: Age threshold in hours (default: 10)
        dry_run: If True, only show what would be compressed without actually doing it
    """
    results_path = Path(results_dir)
    cutoff_time = datetime.now() - timedelta(hours=hours_threshold)
    current_time = datetime.now()
    
    print("=" * 80)
    print("COMPRESSING OLD RESULTS FOLDERS AND FILES")
    print("=" * 80)
    print(f"Current time: {current_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Cutoff time ({hours_threshold} hours ago): {cutoff_time.strftime('%Y-%m-%d %H:%M:%S')}")
    if dry_run:
        print("DRY RUN MODE - No files will be compressed or deleted")
    print()
    
    # Find all directories and files
    folders_to_compress = []
    files_to_compress = []
    
    for item in results_path.iterdir():
        # Skip if already compressed
        if item.name.endswith('.tar.gz') or item.name.endswith('.gz'):
            continue
            
        mtime = datetime.fromtimestamp(item.stat().st_mtime)
        age_hours = (current_time - mtime).total_seconds() / 3600
        
        if age_hours > hours_threshold:
            if item.is_dir():
                folder_size = get_folder_size(item)
                folders_to_compress.append((item, mtime, age_hours, folder_size))
            elif item.is_file():
                file_size = item.stat().st_size
                files_to_compress.append((item, mtime, age_hours, file_size))
    
    # Sort by modification time (oldest first)
    folders_to_compress.sort(key=lambda x: x[1])
    files_to_compress.sort(key=lambda x: x[1])
    
    total_items = len(folders_to_compress) + len(files_to_compress)
    
    if total_items == 0:
        print("No folders or files older than {} hours found.".format(hours_threshold))
        return
    
    # Show summary
    total_folders_size = sum(f[3] for f in folders_to_compress) / (1024**3)  # GB
    total_files_size = sum(f[3] for f in files_to_compress) / (1024**3)  # GB
    total_size = total_folders_size + total_files_size
    
    print(f"Found {len(folders_to_compress)} folders and {len(files_to_compress)} files to compress")
    print(f"Total size: {total_size:.2f} GB ({total_folders_size:.2f} GB folders, {total_files_size:.2f} GB files)\n")
    
    if folders_to_compress:
        print("Folders to compress:")
        for folder, mtime, age_hours, size_bytes in folders_to_compress:
            size_mb = size_bytes / (1024**2)
            print(f"  [FOLDER] {folder.name:50s}  {age_hours:5.1f}h old  {size_mb:7.2f} MB")
    
    if files_to_compress:
        print("\nFiles to compress:")
        for file, mtime, age_hours, size_bytes in files_to_compress:
            size_mb = size_bytes / (1024**2)
            print(f"  [FILE]   {file.name:50s}  {age_hours:5.1f}h old  {size_mb:7.2f} MB")
    
    if dry_run:
        print("\n[DRY RUN] Would compress the above items.")
        return
    
    print("\n" + "=" * 80)
    print("Compressing items...")
    print("=" * 80)
    
    compressed_folders = 0
    compressed_files = 0
    failed_folders = 0
    failed_files = 0
    total_compressed_size = 0
    
    # Compress folders
    if folders_to_compress:
        print("\n--- Compressing Folders ---")
        for folder, mtime, age_hours, size_bytes in folders_to_compress:
            archive_name = results_path / f"{folder.name}.tar.gz"
            
            # Skip if archive already exists
            if archive_name.exists():
                print(f"\n⚠️  Archive already exists: {archive_name.name}")
                continue
            
            print(f"\nCompressing folder: {folder.name} ({size_bytes / (1024**2):.2f} MB)...")
            try:
                # Create tar.gz archive
                result = subprocess.run(
                    ['tar', '-czf', str(archive_name), '-C', str(results_path), folder.name],
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                # Verify archive was created
                if archive_name.exists():
                    archive_size = archive_name.stat().st_size / (1024 * 1024)  # MB
                    compression_ratio = (1 - archive_size / (size_bytes / (1024**2))) * 100 if size_bytes > 0 else 0
                    print(f"  ✓ Created: {archive_name.name} ({archive_size:.2f} MB, {compression_ratio:.1f}% reduction)")
                    
                    # Remove original folder
                    shutil.rmtree(folder)
                    print(f"  ✓ Removed original folder: {folder.name}")
                    
                    compressed_folders += 1
                    total_compressed_size += archive_size
                else:
                    print(f"  ✗ Failed to create archive")
                    failed_folders += 1
                    
            except subprocess.CalledProcessError as e:
                print(f"  ✗ Error compressing {folder.name}: {e.stderr}")
                failed_folders += 1
            except Exception as e:
                print(f"  ✗ Error: {e}")
                failed_folders += 1
    
    # Compress files
    if files_to_compress:
        print("\n--- Compressing Files ---")
        for file, mtime, age_hours, size_bytes in files_to_compress:
            gz_name = results_path / f"{file.name}.gz"
            
            # Skip if already compressed
            if gz_name.exists():
                print(f"\n⚠️  Compressed file already exists: {gz_name.name}")
                continue
            
            print(f"\nCompressing file: {file.name} ({size_bytes / (1024**2):.2f} MB)...")
            try:
                # Create gzip archive (gzip automatically removes original file)
                result = subprocess.run(
                    ['gzip', str(file)],
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                # Verify compressed file was created
                if gz_name.exists():
                    gz_size = gz_name.stat().st_size / (1024 * 1024)  # MB
                    compression_ratio = (1 - gz_size / (size_bytes / (1024**2))) * 100 if size_bytes > 0 else 0
                    print(f"  ✓ Created: {gz_name.name} ({gz_size:.2f} MB, {compression_ratio:.1f}% reduction)")
                    print(f"  ✓ Original file removed: {file.name}")
                    
                    compressed_files += 1
                    total_compressed_size += gz_size
                else:
                    print(f"  ✗ Failed to create compressed file")
                    failed_files += 1
                    
            except subprocess.CalledProcessError as e:
                print(f"  ✗ Error compressing {file.name}: {e.stderr}")
                failed_files += 1
            except Exception as e:
                print(f"  ✗ Error: {e}")
                failed_files += 1
    
    print("\n" + "=" * 80)
    print("COMPRESSION SUMMARY")
    print("=" * 80)
    print(f"Successfully compressed: {compressed_folders} folders, {compressed_files} files")
    if failed_folders > 0 or failed_files > 0:
        print(f"Failed: {failed_folders} folders, {failed_files} files")
    if compressed_folders > 0 or compressed_files > 0:
        print(f"Total compressed size: {total_compressed_size:.2f} MB")
    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description='Compress old results folders using tar.gz'
    )
    parser.add_argument(
        '--results-dir',
        type=str,
        default='results',
        help='Path to results directory (default: results)'
    )
    parser.add_argument(
        '--hours',
        type=int,
        default=10,
        help='Age threshold in hours (default: 10)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be compressed without actually doing it'
    )
    
    args = parser.parse_args()
    
    compress_old_items(
        args.results_dir,
        hours_threshold=args.hours,
        dry_run=args.dry_run
    )


if __name__ == '__main__':
    main()

