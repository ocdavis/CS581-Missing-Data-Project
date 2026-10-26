import pandas as pd
import argparse
import sys

COLOR_PALETTE = ['red', 'blue', 'green', 'yellow', 'cyan', 'magenta', 
                 'orange', 'purple', 'brown', 'pink', 'gray', 'violet']


def get_row_colors(row):
    """Map unique values in a row to colors. Single occurrences stay uncolored."""
    value_counts = row.value_counts()
    colors_map = {}
    color_idx = 0
    
    for value, count in value_counts.items():
        if count > 1:
            colors_map[value] = COLOR_PALETTE[color_idx % len(COLOR_PALETTE)]
            color_idx += 1
        else:
            colors_map[value] = 'white'
    
    return [colors_map[val] for val in row]


def load_data(input_file):
    """Load data from CSV, Excel, JSON, or Parquet file."""
    ext = input_file.split('.')[-1].lower()
    
    if ext == 'csv':
        return pd.read_csv(input_file, dtype={"replicate": str})
    elif ext in ('xlsx', 'xls'):
        return pd.read_excel(input_file)
    elif ext == 'json':
        return pd.read_json(input_file)
    elif ext == 'parquet':
        return pd.read_parquet(input_file)
    else:
        raise ValueError(f"Unsupported file format: {input_file}")


def convert_to_latex(df, caption=None, label=None, decimals=2, output_file=None, index=True):
    """Convert DataFrame to LaTeX table with row-based coloring."""
    latex_lines = []
    
    # Table opening
    latex_lines.append(r'\begin{table}[h]')
    latex_lines.append(r'\centering')
    if caption:
        latex_lines.append(f'\\caption{{{caption}}}')
    if label:
        latex_lines.append(f'\\label{{{label}}}')
    
    # Tabular environment
    col_format = 'l' + 'c' * len(df.columns)
    latex_lines.append(f'\\begin{{tabular}}{{{col_format}}}')
    latex_lines.append(r'\toprule')
    
    # Headers
    header = ('Index & ' if index else '') + ' & '.join(df.columns.astype(str)) + r' \\'
    latex_lines.append(header)
    latex_lines.append(r'\midrule')
    
    # Data rows
    for idx, row in df.iterrows():
        colors = get_row_colors(row)
        cells = [str(idx)] if index else []
        
        for val, color in zip(row, colors):
            formatted = f'{val:.{decimals}f}' if isinstance(val, (int, float)) else str(val)
            cells.append(f'\\cellcolor{{{color}}} {formatted}' if color != 'white' else formatted)
        
        latex_lines.append(' & '.join(cells) + r' \\')
    
    latex_lines.append(r'\bottomrule')
    latex_lines.append(r'\end{tabular}')
    latex_lines.append(r'\end{table}')
    
    output = '\n'.join(latex_lines)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(r'\documentclass{article}' + '\n')
            f.write(r'\usepackage[table]{xcolor}' + '\n')
            f.write(r'\usepackage{booktabs}' + '\n')
            f.write(r'\begin{document}' + '\n\n')
            f.write(output)
            f.write('\n\n' + r'\end{document}' + '\n')
        print(f"LaTeX table written to {output_file}")
    else:
        print(output)


def main():
    parser = argparse.ArgumentParser(description='Convert tabular data to LaTeX tables with colored cells.')
    parser.add_argument('input_file', help='Input file (CSV, XLSX, XLS, JSON, Parquet)')
    parser.add_argument('-o', '--output', help='Output LaTeX file (default: print to stdout)')
    parser.add_argument('-c', '--caption', default='temp', help='Table caption (default: temp)')
    parser.add_argument('-l', '--label', default='tab:temp', help='Table label (default: tab:temp)')
    parser.add_argument('-d', '--decimals', type=int, default=2, help='Number of decimal places (default: 2)')
    parser.add_argument('-ni', '--no-index', action='store_true', help='Do not include index column')
    
    args = parser.parse_args()
    
    try:
        df = load_data(args.input_file)
        convert_to_latex(
            df,
            caption=args.caption,
            label=args.label,
            decimals=args.decimals,
            output_file=args.output,
            index=not args.no_index
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
