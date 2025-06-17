#!/usr/bin/env python3
#Visualize Inspector error statistics
#Creates plots and additional analyses

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse

def plot_error_distributions(df, output_dir):
    """Create various error distribution plots"""
    
    # Set style
    plt.style.use('seaborn-v0_8-darkgrid')
    sns.set_palette("husl")
    
    # 1. Error counts by sample
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Small-scale errors
    sample_data = df.groupby('sample')['small_errors_total'].sum().sort_values(ascending=False)
    sample_data.plot(kind='bar', ax=ax1)
    ax1.set_title('Small-scale Errors by Sample')
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('Total Small-scale Errors')
    ax1.tick_params(axis='x', rotation=45)
    
    # Structural errors
    sample_data = df.groupby('sample')['struct_errors_total'].sum().sort_values(ascending=False)
    sample_data.plot(kind='bar', ax=ax2)
    ax2.set_title('Structural Errors by Sample')
    ax2.set_xlabel('Sample')
    ax2.set_ylabel('Total Structural Errors')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'error_counts_by_sample.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Error rates scatter plot
    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(df['small_errors_per_mbp'], df['struct_errors_per_mbp'], 
                        s=df['assembly_length']/1e7, alpha=0.6)
    ax.set_xlabel('Small-scale Errors per Mbp')
    ax.set_ylabel('Structural Errors per Mbp')
    ax.set_title('Error Rates by Assembly')
    
    # Add sample labels for outliers
    for idx, row in df.iterrows():
        if row['small_errors_per_mbp'] > df['small_errors_per_mbp'].quantile(0.9) or \
           row['struct_errors_per_mbp'] > df['struct_errors_per_mbp'].quantile(0.9):
            ax.annotate(f"{row['sample']}-{row['haplotype']}", 
                       (row['small_errors_per_mbp'], row['struct_errors_per_mbp']),
                       fontsize=8, alpha=0.7)
    
    plt.savefig(os.path.join(output_dir, 'error_rates_scatter.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. QV score distribution
    fig, ax = plt.subplots(figsize=(10, 6))
    df['qv_score'].hist(bins=30, ax=ax)
    ax.axvline(df['qv_score'].mean(), color='red', linestyle='--', 
               label=f'Mean QV: {df["qv_score"].mean():.1f}')
    ax.set_xlabel('QV Score')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of QV Scores')
    ax.legend()
    plt.savefig(os.path.join(output_dir, 'qv_distribution.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Error type breakdown
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Small-scale error types
    small_cols = [col for col in df.columns if col.startswith('small_') and col != 'small_errors_total' 
                  and col != 'small_errors_bp' and col != 'small_errors_per_mbp']
    small_totals = df[small_cols].sum()
    small_totals = small_totals[small_totals > 0]
    
    if not small_totals.empty:
        small_totals.plot(kind='pie', ax=ax1, autopct='%1.1f%%')
        ax1.set_title('Small-scale Error Types')
        ax1.set_ylabel('')
    
    # Structural error types
    struct_cols = [col for col in df.columns if col.startswith('struct_') and col != 'struct_errors_total' 
                   and col != 'struct_errors_bp' and col != 'struct_errors_per_mbp']
    struct_totals = df[struct_cols].sum()
    struct_totals = struct_totals[struct_totals > 0]
    
    if not struct_totals.empty:
        struct_totals.plot(kind='pie', ax=ax2, autopct='%1.1f%%')
        ax2.set_title('Structural Error Types')
        ax2.set_ylabel('')
    
    plt.savefig(os.path.join(output_dir, 'error_type_breakdown.png'), dpi=300, bbox_inches='tight')
    plt.close()

def analyze_error_sizes(small_errors_file, struct_errors_file, output_dir):
    """Analyze error size distributions"""
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    # Small-scale error sizes
    if os.path.exists(small_errors_file):
        small_df = pd.read_csv(small_errors_file, sep='\t')
        if 'size' in small_df.columns:
            small_df['size'].hist(bins=50, ax=ax1)
            ax1.set_xlabel('Error Size (bp)')
            ax1.set_ylabel('Count')
            ax1.set_title('Small-scale Error Size Distribution')
            ax1.set_yscale('log')
            
            # Add statistics
            stats_text = f"Mean: {small_df['size'].mean():.1f} bp\n"
            stats_text += f"Median: {small_df['size'].median():.1f} bp\n"
            stats_text += f"Max: {small_df['size'].max()} bp"
            ax1.text(0.7, 0.9, stats_text, transform=ax1.transAxes, 
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Structural error sizes
    if os.path.exists(struct_errors_file):
        struct_df = pd.read_csv(struct_errors_file, sep='\t')
        if 'size' in struct_df.columns:
            struct_df[struct_df['size'] > 0]['size'].hist(bins=50, ax=ax2)
            ax2.set_xlabel('Error Size (bp)')
            ax2.set_ylabel('Count')
            ax2.set_title('Structural Error Size Distribution')
            ax2.set_yscale('log')
            ax2.set_xscale('log')
            
            # Add statistics
            sizes = struct_df[struct_df['size'] > 0]['size']
            stats_text = f"Mean: {sizes.mean():.1f} bp\n"
            stats_text += f"Median: {sizes.median():.1f} bp\n"
            stats_text += f"Max: {sizes.max()} bp"
            ax2.text(0.05, 0.9, stats_text, transform=ax2.transAxes,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'error_size_distributions.png'), dpi=300, bbox_inches='tight')
    plt.close()

def create_summary_report(df, output_dir):
    """Create a text summary report"""
    
    report_file = os.path.join(output_dir, 'inspector_analysis_report.txt')
    
    with open(report_file, 'w') as f:
        f.write("INSPECTOR ERROR ANALYSIS REPORT\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("OVERALL STATISTICS\n")
        f.write("-" * 30 + "\n")
        f.write(f"Total samples: {df['sample'].nunique()}\n")
        f.write(f"Total assemblies analyzed: {len(df)}\n")
        f.write(f"Mean assembly size: {df['assembly_length'].mean() / 1e9:.2f} Gbp\n")
        f.write(f"Mean N50: {df['n50'].mean() / 1e6:.1f} Mbp\n")
        f.write(f"Mean QV score: {df['qv_score'].mean():.2f} ± {df['qv_score'].std():.2f}\n")
        f.write(f"Mean mapping rate: {df['mapping_rate'].mean():.1f}%\n\n")
        
        f.write("ERROR STATISTICS\n")
        f.write("-" * 30 + "\n")
        f.write(f"Total small-scale errors: {df['small_errors_total'].sum():,}\n")
        f.write(f"Total structural errors: {df['struct_errors_total'].sum():,}\n")
        f.write(f"Mean small-scale error rate: {df['small_errors_per_mbp'].mean():.2f} per Mbp\n")
        f.write(f"Mean structural error rate: {df['struct_errors_per_mbp'].mean():.4f} per Mbp\n")
        f.write(f"Total error bp: {df['total_error_bp'].sum():,}\n")
        f.write(f"Mean error fraction: {df['error_fraction'].mean():.6f}\n\n")
        
        f.write("TOP 5 SAMPLES BY ERROR RATE\n")
        f.write("-" * 30 + "\n")
        top_error = df.nlargest(5, 'error_fraction')[['sample', 'haplotype', 'error_fraction', 'total_errors']]
        f.write(top_error.to_string(index=False))
        f.write("\n\n")
        
        f.write("TOP 5 SAMPLES BY QV SCORE\n")
        f.write("-" * 30 + "\n")
        top_qv = df.nlargest(5, 'qv_score')[['sample', 'haplotype', 'qv_score', 'total_errors']]
        f.write(top_qv.to_string(index=False))
        f.write("\n\n")
        
        f.write("ERROR TYPE BREAKDOWN\n")
        f.write("-" * 30 + "\n")
        
        # Small-scale errors
        small_cols = [col for col in df.columns if col.startswith('small_') and 
                      col not in ['small_errors_total', 'small_errors_bp', 'small_errors_per_mbp']]
        for col in small_cols:
            if col in df.columns and df[col].sum() > 0:
                f.write(f"{col}: {df[col].sum():,}\n")
        
        f.write("\n")
        
        # Structural errors
        struct_cols = [col for col in df.columns if col.startswith('struct_') and 
                       col not in ['struct_errors_total', 'struct_errors_bp', 'struct_errors_per_mbp']]
        for col in struct_cols:
            if col in df.columns and df[col].sum() > 0:
                f.write(f"{col}: {df[col].sum():,}\n")
    
    print(f"Report saved to {report_file}")

def main():
    parser = argparse.ArgumentParser(description='Visualize Inspector error statistics')
    parser.add_argument('--data-dir', default='/projects/standard/hsiehph/shared/DIR_homes/hudso501/analysisPI/data',
                        help='Directory containing the analysis output files')
    args = parser.parse_args()
    
    # Load summary data
    summary_file = os.path.join(args.data_dir, 'inspector_error_summary.tsv')
    if not os.path.exists(summary_file):
        print(f"Error: Could not find {summary_file}")
        print("Please run analyze_inspector_errors.py first")
        return
    
    df = pd.read_csv(summary_file, sep='\t')
    print(f"Loaded data for {len(df)} assemblies")
    
    # Create plots
    print("Creating visualizations...")
    plot_error_distributions(df, args.data_dir)
    
    # Analyze error sizes
    small_errors_file = os.path.join(args.data_dir, 'all_small_scale_errors.tsv')
    struct_errors_file = os.path.join(args.data_dir, 'all_structural_errors.tsv')
    analyze_error_sizes(small_errors_file, struct_errors_file, args.data_dir)
    
    # Create summary report
    create_summary_report(df, args.data_dir)
    
    print("Visualization complete!")

if __name__ == "__main__":
    main()
