#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge pi values and Fst values for two populations
Usage: python 00.merge_pi_Fst.py pop1 pop1_pi_file pop2 pop2_pi_file fst_file output_file
Example: python 00.merge_pi_Fst.py QH QH.win50Kstep5K.pi.windowed.pi QL QL.win50Kstep5K.pi.windowed.pi QH_QL.win50Kstep5k.fst.windowed.weir.fst QH_QL.win50Kstep5k.pi.fst.tab
"""

import sys
import pandas as pd
import numpy as np
from collections import defaultdict
import os

def read_pi_file(pi_file, pop_name):
    """
    Read pi file with format: CHROM  BIN_START  BIN_END  N_VARIANTS  PI
    """
    pi_data = {}
    try:
        print(f"Reading pi file for {pop_name}: {pi_file}")
        with open(pi_file, 'r') as f:
            # Skip header line
            header = next(f).strip()
            print(f"File header: {header}")
            
            for line_num, line in enumerate(f, 2):
                parts = line.strip().split()
                if len(parts) < 5:
                    print(f"Warning: Line {line_num} has incorrect format: {line.strip()}")
                    continue
                
                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    n_variants = int(parts[3])
                    pi_value = float(parts[4])
                except ValueError as e:
                    print(f"Warning: Line {line_num} value conversion error: {e}")
                    continue
                
                # Use chrom:start-end as key
                key = f"{chrom}:{start}-{end}"
                
                pi_data[key] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'n_variants': n_variants,
                    f'pi_{pop_name}': pi_value,
                    f'n_var_{pop_name}': n_variants
                }
        
        print(f"Successfully read pi file for {pop_name}: {len(pi_data)} windows")
        return pi_data
    except Exception as e:
        print(f"Error reading pi file {pi_file}: {e}")
        sys.exit(1)

def read_fst_file(fst_file):
    """
    Read Fst file with format: CHROM  BIN_START  BIN_END  N_VARIANTS  WEIGHTED_FST  MEAN_FST
    """
    fst_data = {}
    try:
        print(f"Reading Fst file: {fst_file}")
        with open(fst_file, 'r') as f:
            # Skip header line
            header = next(f).strip()
            print(f"Fst file header: {header}")
            
            for line_num, line in enumerate(f, 2):
                parts = line.strip().split()
                if len(parts) < 6:
                    print(f"Warning: Fst file line {line_num} has incorrect format: {line.strip()}")
                    continue
                
                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    n_variants = int(parts[3])
                    weighted_fst = float(parts[4])
                    mean_fst = float(parts[5])
                except ValueError as e:
                    print(f"Warning: Fst file line {line_num} value conversion error: {e}")
                    continue
                
                # Use chrom:start-end as key
                key = f"{chrom}:{start}-{end}"
                
                fst_data[key] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'fst_n_variants': n_variants,
                    'weighted_fst': weighted_fst,
                    'mean_fst': mean_fst
                }
        
        print(f"Successfully read Fst file: {len(fst_data)} windows")
        return fst_data
    except Exception as e:
        print(f"Error reading Fst file {fst_file}: {e}")
        sys.exit(1)

def merge_data(pi_data1, pop1_name, pi_data2, pop2_name, fst_data):
    """
    Merge pi and Fst data
    """
    print("\nStarting data merge...")
    
    # Collect all windows
    all_keys = set(list(pi_data1.keys()) + list(pi_data2.keys()) + list(fst_data.keys()))
    print(f"Total unique windows: {len(all_keys)}")
    
    merged_data = []
    
    # First sort by natural chromosome and position order
    # Create a list to store sorted keys
    sorted_keys = []
    
    for key in all_keys:
        if ':' in key and '-' in key:
            chrom, pos_range = key.split(':')
            start, end = map(int, pos_range.split('-'))
            sorted_keys.append((chrom, start, end, key))
        else:
            print(f"Warning: Incorrect key format: {key}")
    
    # Custom sorting function for chromosome numbers (chr1, chr2, ..., chr10, chrX, chrY, etc.)
    def chromosome_key(chrom):
        # Remove possible 'chr' prefix
        chrom_str = str(chrom).lower().replace('chr', '')
    
    # Sort by chromosome and start position
    sorted_keys.sort(key=lambda x: (chromosome_key(x[0]), x[1]))
    
    for chrom, start, end, key in sorted_keys:
        # Get data from each source
        pi1 = pi_data1.get(key, {})
        pi2 = pi_data2.get(key, {})
        fst = fst_data.get(key, {})
        
        # Prepare row data
        row = {
            'CHROM': chrom,
            'BIN_START': start,
            'BIN_END': end,
            'POS_ID': key
        }
        
        # Add pi values
        row[f'{pop1_name}_PI'] = pi1.get(f'pi_{pop1_name}', np.nan)
        row[f'{pop2_name}_PI'] = pi2.get(f'pi_{pop2_name}', np.nan)
        
        # Calculate PI_RATIO
        pi1_val = pi1.get(f'pi_{pop1_name}', 0)
        pi2_val = pi2.get(f'pi_{pop2_name}', 0)
        if pi2_val != 0:
            row['PI_RATIO'] = pi1_val / pi2_val
        else:
            row['PI_RATIO'] = np.nan
        
        # Add Fst values
        row['WEIGHTED_FST'] = fst.get('weighted_fst', np.nan)
        row['MEAN_FST'] = fst.get('mean_fst', np.nan)
        
        # Add variant counts (optional)
        row[f'{pop1_name}_N_VAR'] = pi1.get(f'n_var_{pop1_name}', 0)
        row[f'{pop2_name}_N_VAR'] = pi2.get(f'n_var_{pop2_name}', 0)
        row['FST_N_VAR'] = fst.get('fst_n_variants', 0)
        
        merged_data.append(row)
    
    # Convert to DataFrame
    df = pd.DataFrame(merged_data)
    
    # Reorder columns
    base_cols = ['CHROM', 'BIN_START', 'BIN_END', 'POS_ID', 
                 f'{pop1_name}_PI', f'{pop2_name}_PI', 'PI_RATIO',
                 'WEIGHTED_FST', 'MEAN_FST']
    
    # Keep only existing columns
    final_cols = [col for col in base_cols if col in df.columns]
    df = df[final_cols]
    
    return df


def save_output(df, output_file):
    """
    Save results to file
    """
    try:
        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Save file
        df.to_csv(output_file, sep='\t', index=False, na_rep='NA')
        print(f"\nResults saved to: {output_file}")
        print(f"Total rows: {len(df)}")
        
        # Display data preview
        print("\nData preview (first 5 rows):")
        print(df.head().to_string())
        
        # Display basic statistics
        print("\nBasic statistics:")
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            stats = df[numeric_cols].describe()
            print(stats.round(4).to_string())
        
        # Check for missing values
        print("\nMissing value statistics:")
        missing = df.isnull().sum()
        missing_pct = (missing / len(df) * 100).round(2)
        for col in df.columns:
            if missing[col] > 0:
                print(f"  {col}: {missing[col]} rows ({missing_pct[col]}%)")
    
    except Exception as e:
        print(f"Error saving results file: {e}")
        sys.exit(1)

def main():
    # Check arguments
    if len(sys.argv) != 7:
        print("Error: Incorrect number of arguments!")
        print("Usage: python 00.merge_pi_Fst.py pop1 pop1_pi_file pop2 pop2_pi_file fst_file output_file")
        print("Example: python 00.merge_pi_Fst.py QH QH.pi.txt QL QL.pi.txt QH_QL.fst.txt QH_QL.merged.tab")
        sys.exit(1)
    
    # Parse arguments
    pop1_name = sys.argv[1]
    pop1_pi_file = sys.argv[2]
    pop2_name = sys.argv[3]
    pop2_pi_file = sys.argv[4]
    fst_file = sys.argv[5]
    output_file = sys.argv[6]
    
    print("=" * 60)
    print("Script parameters:")
    print(f"Population 1 name: {pop1_name}")
    print(f"Population 1 pi file: {pop1_pi_file}")
    print(f"Population 2 name: {pop2_name}")
    print(f"Population 2 pi file: {pop2_pi_file}")
    print(f"Fst file: {fst_file}")
    print(f"Output file: {output_file}")
    print("=" * 60)
    
    # Check if input files exist
    for file_path in [pop1_pi_file, pop2_pi_file, fst_file]:
        if not os.path.exists(file_path):
            print(f"Error: File does not exist: {file_path}")
            sys.exit(1)
    
    # Read data
    pi_data1 = read_pi_file(pop1_pi_file, pop1_name)
    pi_data2 = read_pi_file(pop2_pi_file, pop2_name)
    fst_data = read_fst_file(fst_file)
    
    # Merge data
    merged_df = merge_data(pi_data1, pop1_name, pi_data2, pop2_name, fst_data)
    
    # Save results
    save_output(merged_df, output_file)
    
    print("\n" + "=" * 60)
    print("Processing complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()