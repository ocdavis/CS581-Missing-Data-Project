#!/usr/bin/env python3
"""
Convert tabular data (CSV, Excel, JSON, Parquet) to a pandas DataFrame 
and export as a LaTeX table.
"""

import argparse
import sys
import pandas as pd


def load_file(filepath):
    """Load data from various file formats into a pandas DataFrame."""
    file_ext = filepath.split('.')[-1].lower()
    
    try:
        if file_ext == 'csv':
            return pd.read_csv(filepath, index_col=0, dtype={"replicate": str})
        elif file_ext in ['xlsx', 'xls']:
            return pd.read_excel(filepath, index_col=0)
        elif file_ext == 'json':
            return pd.read_json(filepath, index_col=0)
        elif file_ext == 'parquet':
            return pd.read_parquet(filepath, index_col=0)
        else:
            raise ValueError(f"Unsupported file format: .{file_ext}")
    except Exception as e:
        print(f"Error loading file: {e}", file=sys.stderr)
        sys.exit(1)


def convert_to_latex(df, caption=None, label=None, decimals=2, output_file=None, index=True):
    """Convert DataFrame to LaTeX table and optionally save to file."""
    try:
        latex_table = df.to_latex(
            float_format=lambda x: f'{x:.{decimals}f}',
            caption=caption,
            label=label,
            column_format='l' + 'c' * (len(df.columns)),
            index=index
        )
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(latex_table)
            print(f"LaTeX table saved to: {output_file}")
        else:
            print(latex_table)
            
    except Exception as e:
        print(f"Error converting to LaTeX: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Convert tabular data (CSV, Excel, JSON, Parquet) to LaTeX table.',
        epilog='Supported formats: CSV, XLSX, XLS, JSON, Parquet'
    )
    
    parser.add_argument(
        'input_file',
        help='Path to the input file (CSV, XLSX, XLS, JSON, or Parquet)'
    )
    
    parser.add_argument(
        '-o', '--output',
        dest='output_file',
        default=None,
        help='Output file path for LaTeX table (default: print to stdout)'
    )
    
    parser.add_argument(
        '-c', '--caption',
        default="temp",
        help='LaTeX table caption'
    )
    
    parser.add_argument(
        '-l', '--label',
        default="tab:temp",
        help='LaTeX table label (for cross-referencing)'
    )
    
    parser.add_argument(
        '-d', '--decimals',
        type=int,
        default=2,
        help='Number of decimal places to display (default: 2)'
    )
    
    parser.add_argument(
        '-ni','--no-index',
        action='store_false',
        help='Remove the index'
    )
    
    args = parser.parse_args()
    
    # Load the file
    print(f"Loading {args.input_file}...")
    df = load_file(args.input_file)
    
    print(f"DataFrame shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    
    # Convert to LaTeX
    convert_to_latex(
        df,
        caption=args.caption,
        label=args.label,
        decimals=args.decimals,
        output_file=args.output_file,
        index=args.no_index
    )


if __name__ == '__main__':
    main()
