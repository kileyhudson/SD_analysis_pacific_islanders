#!/usr/bin/env python3
#Analyze Inspector error statistics 

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import argparse
from collections import defaultdict
import re

def parse_summary_statistics(file_path):
    """Parse the summary statistics file"""
    stats = {}
    
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Extract key statistics using regex
    patterns = {
        'total_length': r'Total length (\d+)',
        'n50': r'N50 (\d+)',
        'num_contigs': r'Number of contigs (\d+)',
        'longest_contig': r'Longest contig (\d+)',
        'qv_score': r'QV ([\d.]+)',
        'mapping_rate': r'Mapping rate /% ([\d.]+)',
        'depth': r'Depth ([\d.]+)',
        'structural_errors': r'Structural error (\d+)',
        'small_scale_errors_total': r'Total small-scale assembly error (\d+)',
        'small_scale_errors_per_mbp': r'Small-scale assembly error /per Mbp ([\d.]+)'
    }
    
    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            stats[key] = float(match.group(1)) if '.' in match.group(1) else int(match.group(1))
    
    # Extract error type counts
    error_patterns = {
        'expansion': r'Expansion (\d+)',
        'collapse': r'Collapse (\d+)',
        'haplotype_switch': r'Haplotype switch (\d+)',
        'inversion': r'Inversion (\d+)',
        'base_substitution': r'Base substitution (\d+)',
        'small_expansion': r'Small-scale expansion (\d+)',
        'small_collapse': r'Small-scale collapse (\d+)'
    }
    
    for key, pattern in error_patterns.items():
        match = re.search(pattern, content)
        if match:
            stats[key] = int(match.group(1))
    
    return stats

def parse_small_scale_errors(file_path):
    """Parse small-scale error BED file"""
    errors = []
    
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found")
        return pd.DataFrame()
    
    # Read the file
    df = pd.read_csv(file_path, sep='\t', comment='#', 
                     names=['contig', 'start', 'end', 'base_contig', 'base_read', 
                            'supporting_reads', 'depth', 'type', 'pvalue'])
    
    # Calculate error sizes
    df['size'] = df['end'] - df['start']
    
    # For substitutions, size is 1
    df.loc[df['type'] == 'BaseSubstitution', 'size'] = 1
    
    return df

def parse_structural_errors(file_path):
    """Parse structural error BED file"""
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found")
        return pd.DataFrame()
    
    # Read the file
    df = pd.read_csv(file_path, sep='\t', comment='#',
                     names=['contig', 'start', 'end', 'supporting_reads', 
                            'type', 'size_info', 'haplotype_info', 'depth_left', 
                            'depth_right', 'depth_min', 'read_names', 'hap_switch_info'])
    
    # Extract size from size_info column
    def extract_size(size_info):
        if pd.isna(size_info):
            return 0
        match = re.search(r'Size=(\d+)', str(size_info))
        if match:
            return int(match.group(1))
        # For HaplotypeSwitch, get the first size
        match = re.search(r'Size=(\d+);', str(size_info))
        if match:
            return int(match.group(1))
        return 0
    
    df['size'] = df['size_info'].apply(extract_size)
    
    # Handle HaplotypeSwitch positions (they have semicolon-separated values)
    def parse_position(pos):
        if ';' in str(pos):
            # Return the first position for simplicity
            return int(str(pos).split(';')[0])
        return int(pos)
    
    df['start'] = df['start'].apply(parse_position)
    df['end'] = df['end'].apply(parse_position)
    
    return df

def calculate_error_statistics(small_errors, struct_errors, summary_stats):
    """Calculate comprehensive error statistics"""
    stats = {}
    
    # Total assembly length
    total_length = summary_stats.get('total_length', 0)
    total_length_mbp = total_length / 1_000_000
    
    # Small-scale error statistics
    if not small_errors.empty:
        small_error_types = small_errors['type'].value_counts().to_dict()
        stats['small_scale_errors'] = {
            'total': len(small_errors),
            'types': small_error_types,
            'total_bp': small_errors['size'].sum(),
            'mean_size': small_errors['size'].mean(),
            'median_size': small_errors['size'].median(),
            'errors_per_mbp': len(small_errors) / total_length_mbp if total_length_mbp > 0 else 0
        }
    else:
        stats['small_scale_errors'] = {'total': 0, 'types': {}, 'total_bp': 0}
    
    # Structural error statistics
    if not struct_errors.empty:
        struct_error_types = struct_errors['type'].value_counts().to_dict()
        stats['structural_errors'] = {
            'total': len(struct_errors),
            'types': struct_error_types,
            'total_bp': struct_errors['size'].sum(),
            'mean_size': struct_errors['size'].mean(),
            'median_size': struct_errors['size'].median(),
            'errors_per_mbp': len(struct_errors) / total_length_mbp if total_length_mbp > 0 else 0
        }
    else:
        stats['structural_errors'] = {'total': 0, 'types': {}, 'total_bp': 0}
    
    # Combined statistics
    all_errors_bp = stats['small_scale_errors']['total_bp'] + stats['structural_errors']['total_bp']
    stats['combined'] = {
        'total_errors': stats['small_scale_errors']['total'] + stats['structural_errors']['total'],
        'total_error_bp': all_errors_bp,
        'error_fraction': all_errors_bp / total_length if total_length > 0 else 0,
        'assembly_length': total_length,
        'assembly_length_mbp': total_length_mbp
    }
    
    # Add summary statistics
    stats['assembly_metrics'] = {
        'n50': summary_stats.get('n50', 0),
        'num_contigs': summary_stats.get('num_contigs', 0),
        'qv_score': summary_stats.get('qv_score', 0),
        'mapping_rate': summary_stats.get('mapping_rate', 0),
        'depth': summary_stats.get('depth', 0)
    }
    
    return stats

def analyze_sample(sample_dir, sample_name, haplotype):
    """Analyze a single sample/haplotype combination"""
    hap_dir = os.path.join(sample_dir, haplotype)
    
    # Parse files
    summary_stats = parse_summary_statistics(os.path.join(hap_dir, 'summary_statistics'))
    small_errors = parse_small_scale_errors(os.path.join(hap_dir, 'small_scale_error.bed'))
    struct_errors = parse_structural_errors(os.path.join(hap_dir, 'structural_error.bed'))
    
    # Calculate statistics
    stats = calculate_error_statistics(small_errors, struct_errors, summary_stats)
    
    # Add sample information
    stats['sample_info'] = {
        'sample': sample_name,
        'haplotype': haplotype
    }
    
    return stats, small_errors, struct_errors

def create_summary_table(all_stats):
    """Create a summary table from all statistics"""
    rows = []
    
    for stats in all_stats:
        row = {
            'sample': stats['sample_info']['sample'],
            'haplotype': stats['sample_info']['haplotype'],
            'assembly_length': stats['combined']['assembly_length'],
            'n50': stats['assembly_metrics']['n50'],
            'num_contigs': stats['assembly_metrics']['num_contigs'],
            'qv_score': stats['assembly_metrics']['qv_score'],
            'mapping_rate': stats['assembly_metrics']['mapping_rate'],
            'depth': stats['assembly_metrics']['depth'],
            'small_errors_total': stats['small_scale_errors']['total'],
            'small_errors_bp': stats['small_scale_errors']['total_bp'],
            'small_errors_per_mbp': stats['small_scale_errors']['errors_per_mbp'],
            'struct_errors_total': stats['structural_errors']['total'],
            'struct_errors_bp': stats['structural_errors']['total_bp'],
            'struct_errors_per_mbp': stats['structural_errors']['errors_per_mbp'],
            'total_errors': stats['combined']['total_errors'],
            'total_error_bp': stats['combined']['total_error_bp'],
            'error_fraction': stats['combined']['error_fraction']
        }
        
        # Add error type counts
        for error_type, count in stats['small_scale_errors']['types'].items():
            row[f'small_{error_type}'] = count
        
        for error_type, count in stats['structural_errors']['types'].items():
            row[f'struct_{error_type}'] = count
        
        rows.append(row)
    
    return pd.DataFrame(rows)

def main():
    parser = argparse.ArgumentParser(description='Analyze Inspector error statistics')
    parser.add_argument('--input-dir', default='/projects/standard/hsiehph/shared/globus-incoming/assembly_qc_files',
                        help='Input directory containing sample folders')
    parser.add_argument('--output-dir', default='/projects/standard/hsiehph/shared/DIR_homes/hudso501/analysisPI/data',
                        help='Output directory for results')
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Get all sample directories
    sample_dirs = sorted([d for d in os.listdir(args.input_dir) if os.path.isdir(os.path.join(args.input_dir, d))])
    
    print(f"Found {len(sample_dirs)} samples to analyze")
    
    all_stats = []
    all_small_errors = []
    all_struct_errors = []
    
    # Analyze each sample
    for sample_dir_name in sample_dirs:
        sample_path = os.path.join(args.input_dir, sample_dir_name)
        sample_name = sample_dir_name.split('_', 1)[1]  # Remove number prefix
        
        print(f"Analyzing {sample_name}...")
        
        for haplotype in ['hap1', 'hap2']:
            try:
                stats, small_errors, struct_errors = analyze_sample(sample_path, sample_name, haplotype)
                all_stats.append(stats)
                
                # Add sample info to error dataframes
                if not small_errors.empty:
                    small_errors['sample'] = sample_name
                    small_errors['haplotype'] = haplotype
                    all_small_errors.append(small_errors)
                
                if not struct_errors.empty:
                    struct_errors['sample'] = sample_name
                    struct_errors['haplotype'] = haplotype
                    all_struct_errors.append(struct_errors)
                    
            except Exception as e:
                print(f"Error processing {sample_name} {haplotype}: {str(e)}")
    
    # Create summary table
    summary_df = create_summary_table(all_stats)
    summary_df = summary_df.sort_values(['sample', 'haplotype'])
    
    # Save summary table
    summary_file = os.path.join(args.output_dir, 'inspector_error_summary.tsv')
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"Saved summary to {summary_file}")
    
    # Combine all errors
    if all_small_errors:
        combined_small = pd.concat(all_small_errors, ignore_index=True)
        small_errors_file = os.path.join(args.output_dir, 'all_small_scale_errors.tsv')
        combined_small.to_csv(small_errors_file, sep='\t', index=False)
        print(f"Saved {len(combined_small)} small-scale errors to {small_errors_file}")
    
    if all_struct_errors:
        combined_struct = pd.concat(all_struct_errors, ignore_index=True)
        struct_errors_file = os.path.join(args.output_dir, 'all_structural_errors.tsv')
        combined_struct.to_csv(struct_errors_file, sep='\t', index=False)
        print(f"Saved {len(combined_struct)} structural errors to {struct_errors_file}")
    
    # Print summary statistics
    print("\n=== OVERALL STATISTICS ===")
    print(f"Total samples analyzed: {len(sample_dirs)}")
    print(f"Total haplotypes analyzed: {len(all_stats)}")
    print(f"Mean QV score: {summary_df['qv_score'].mean():.2f}")
    print(f"Mean assembly size: {summary_df['assembly_length'].mean() / 1e9:.2f} Gbp")
    print(f"Total small-scale errors: {summary_df['small_errors_total'].sum()}")
    print(f"Total structural errors: {summary_df['struct_errors_total'].sum()}")
    print(f"Mean error rate: {summary_df['error_fraction'].mean():.6f}")
    
    # Save detailed statistics
    stats_file = os.path.join(args.output_dir, 'inspector_detailed_stats.txt')
    with open(stats_file, 'w') as f:
        f.write("INSPECTOR ERROR ANALYSIS - DETAILED STATISTICS\n")
        f.write("=" * 50 + "\n\n")
        
        for stats in all_stats:
            f.write(f"Sample: {stats['sample_info']['sample']} - {stats['sample_info']['haplotype']}\n")
            f.write(f"Assembly length: {stats['combined']['assembly_length']:,} bp\n")
            f.write(f"Total errors: {stats['combined']['total_errors']}\n")
            f.write(f"Total error bp: {stats['combined']['total_error_bp']:,}\n")
            f.write(f"Error fraction: {stats['combined']['error_fraction']:.6f}\n")
            f.write("\nSmall-scale errors by type:\n")
            for error_type, count in stats['small_scale_errors']['types'].items():
                f.write(f"  {error_type}: {count}\n")
            f.write("\nStructural errors by type:\n")
            for error_type, count in stats['structural_errors']['types'].items():
                f.write(f"  {error_type}: {count}\n")
            f.write("\n" + "-" * 30 + "\n\n")
    
    print(f"\nDetailed statistics saved to {stats_file}")

if __name__ == "__main__":
    main()
