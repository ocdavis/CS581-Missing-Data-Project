import pandas as pd
import sys
import argparse
from collections import defaultdict

def analyze_phylogeny_scores(ranking_file, scores_file, output_file=None, drop_columns=None):
    """Analyzes ranking consistency and returns results as a DataFrame.
    
    Args:
        ranking_file: Path to rankings CSV file
        scores_file: Path to scores CSV file
        output_file: Optional output file path (auto-detects format from extension)
        drop_columns: List of columns to drop from results (e.g., ['Ranking', 'Method Count'])
    """
    
    rankings_df = pd.read_csv(ranking_file, dtype={"replicate": str})
    scores_df = pd.read_csv(scores_file, dtype={"replicate": str})
    
    method_columns = [col for col in rankings_df.columns if col != 'replicate']
    
    # Collect results in a list of dicts
    results = []
    all_ranges = []
    
    for idx, row in rankings_df.iterrows():
        replicate = row['replicate']
        
        ranking_groups = defaultdict(list)
        for method in method_columns:
            ranking_value = row[method]
            ranking_groups[ranking_value].append(method)
        
        replicate_ranges = []
        
        for ranking_value, methods in sorted(ranking_groups.items()):
            if len(methods) > 1:
                scores = scores_df[
                    (scores_df['replicate'] == replicate) & 
                    (scores_df['method_label'].isin(methods))
                ]['mpl_score'].values
                
                if len(scores) > 0:
                    score_range = max(scores) - min(scores)
                    replicate_ranges.append(score_range)
                    
                    results.append({
                        'Replicate': replicate,
                        'Ranking': ranking_value,
                        'Methods': ', '.join(methods),
                        'Method Count': len(methods),
                        'Score Range': score_range,
                        'Min Score': min(scores),
                        'Max Score': max(scores)
                    })
                    all_ranges.append(score_range)
    
    # Create DataFrame
    results_df = pd.DataFrame(results)
    
    # Drop specified columns
    if drop_columns:
        results_df = results_df.drop(columns=[col for col in drop_columns if col in results_df.columns])
    
    # Print console output
    print(results_df.to_string(index=False))
    print()
    
    if all_ranges:
        print(f"Overall average range: {sum(all_ranges) / len(all_ranges):.4f}")
        print(f"Total ranking groups analyzed: {len(all_ranges)}")
    else:
        print("No ranking groups with multiple methods found.")
    
    # Optionally save to file
    if output_file:
        file_ext = output_file.split('.')[-1].lower()
        try:
            if file_ext == 'csv':
                results_df.to_csv(output_file, index=False)
            elif file_ext in ['xlsx', 'xls']:
                results_df.to_excel(output_file, index=False)
            elif file_ext == 'html':
                results_df.to_html(output_file, index=False)
            else:
                raise ValueError(f"Unsupported output format: .{file_ext}")
            print(f"\nResults saved to {output_file}")
        except Exception as e:
            print(f"Error saving file: {e}", file=sys.stderr)
            sys.exit(1)
    
    return results_df


def main():
    parser = argparse.ArgumentParser(
        description='Analyze phylogeny ranking consistency and score ranges.',
        epilog='Supported output formats: CSV, XLSX, HTML'
    )
    
    parser.add_argument(
        'ranking_file',
        help='Path to rankings CSV file'
    )
    
    parser.add_argument(
        'scores_file',
        help='Path to scores CSV file'
    )
    
    parser.add_argument(
        '-o', '--output',
        dest='output_file',
        default=None,
        help='Output file path (auto-detects format from extension: .csv, .xlsx, .html)'
    )
    
    parser.add_argument(
        '-dr', '--drop-ranking',
        action='store_true',
        help='Drop the Ranking column from output'
    )
    
    parser.add_argument(
        '-dm', '--drop-method-count',
        action='store_true',
        help='Drop the Method Count column from output'
    )
    
    parser.add_argument(
        '--drop-both',
        action='store_true',
        help='Drop both Ranking and Method Count columns'
    )
    
    args = parser.parse_args()
    
    # Determine which columns to drop
    drop_columns = []
    if args.drop_both or (args.drop_ranking and args.drop_method_count):
        drop_columns = ['Ranking', 'Method Count']
    else:
        if args.drop_ranking:
            drop_columns.append('Ranking')
        if args.drop_method_count:
            drop_columns.append('Method Count')
    
    analyze_phylogeny_scores(
        args.ranking_file,
        args.scores_file,
        args.output_file,
        drop_columns=drop_columns if drop_columns else None
    )


if __name__ == "__main__":
    main()
