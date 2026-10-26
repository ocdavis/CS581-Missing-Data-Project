import argparse
import pandas as pd
from utils import load_config

CONFIG_HELP = r"""
The TOML configuration file must contain:

networks_csv = "/path/to/networks.csv"
errors_csv = "/path/to/network_errors.csv"
mpl_csv = "/path/to/mpl_scores.csv"
output = "/path/to/output/directory"
method1 = "CAMUS"
method2 = "PhyloNET"

Example:

    networks_csv = "/data/networks.csv"
    errors_csv = "/data/network_errors.csv"
    mpl_csv = "/data/mpl_scores.csv"
    output = "/data/results"
    method1 = "CAMUS"
    method2 = "PhyloNET"
"""

REQUIRED_KEYS = {
    "networks_csv",
    "errors_csv",
    "mpl_csv",
    "method1",
    "method2"
}


def load_and_prepare_data(networks_csv, errors_csv, mpl_csv):
    """Load all CSV files and return DataFrames"""
    df_networks = pd.read_csv(networks_csv)
    df_errors = pd.read_csv(errors_csv)
    df_mpl = pd.read_csv(mpl_csv)
    return df_networks, df_errors, df_mpl


def find_network_differences(df_networks, method1, method2):
    """Identify replicates where networks differ between two methods"""
    df_networks['networks_differ'] = df_networks[method1] != df_networks[method2]
    num_different = df_networks['networks_differ'].sum()
    fraction_different = num_different / len(df_networks)
    different_replicates = df_networks[df_networks['networks_differ']]['replicate'].tolist()
    
    return num_different, fraction_different, different_replicates


def compare_metrics(df_errors, df_mpl, method1, method2, different_replicates):
    """Compare error and MPL metrics between two methods for differing networks"""
    
    # Prepare error data
    df_errors_method1 = df_errors[df_errors['method_label'] == method1][['replicate', 'false_average']].rename(
        columns={'false_average': f'error_{method1.lower()}'}
    )
    df_errors_method2 = df_errors[df_errors['method_label'] == method2][['replicate', 'false_average']].rename(
        columns={'false_average': f'error_{method2.lower()}'}
    )
    
    # Prepare MPL data
    df_mpl_method1 = df_mpl[df_mpl['method_label'] == method1][['replicate', 'mpl_score']].rename(
        columns={'mpl_score': f'mpl_{method1.lower()}'}
    )
    df_mpl_method2 = df_mpl[df_mpl['method_label'] == method2][['replicate', 'mpl_score']].rename(
        columns={'mpl_score': f'mpl_{method2.lower()}'}
    )
    
    # Merge error data
    error_comparison = df_errors_method1.merge(df_errors_method2, on='replicate', how='outer')
    error_comparison['error_improvement'] = error_comparison[f'error_{method1.lower()}'] - error_comparison[f'error_{method2.lower()}']
    
    # Merge MPL data
    mpl_comparison = df_mpl_method1.merge(df_mpl_method2, on='replicate', how='outer')
    mpl_comparison['mpl_improvement'] = mpl_comparison[f'mpl_{method2.lower()}'] - mpl_comparison[f'mpl_{method1.lower()}']
    
    # Filter for only replicates where networks differ
    error_comparison_diff = error_comparison[error_comparison['replicate'].isin(different_replicates)].dropna()
    mpl_comparison_diff = mpl_comparison[mpl_comparison['replicate'].isin(different_replicates)].dropna()
    
    return error_comparison_diff, mpl_comparison_diff


def print_results(num_different, fraction_different, error_comparison, mpl_comparison, method1, method2, output_path):
    """Print and save analysis results"""
    
    output_text = []
    output_text.append(f"Network Comparison: {method1} vs {method2}\n")
    output_text.append(f"{'='*60}\n")
    output_text.append(f"Networks differ: {num_different} out of {num_different + (len(error_comparison) - num_different)} replicates\n")
    output_text.append(f"Fraction different: {fraction_different:.2%}\n\n")
    
    # Error analysis
    if not error_comparison.empty:
        avg_error_improvement = error_comparison['error_improvement'].mean()
        output_text.append(f"Network Error Analysis (false_average - lower is better)\n")
        output_text.append(f"{'-'*60}\n")
        output_text.append(f"Average error improvement: {avg_error_improvement:.6f}\n")
        output_text.append(f"(Positive = {method2} is better)\n\n")
        output_text.append("Detailed Error Comparison:\n")
        output_text.append(error_comparison.to_string(index=False))
        output_text.append("\n\n")
    
    # MPL analysis
    if not mpl_comparison.empty:
        avg_mpl_improvement = mpl_comparison['mpl_improvement'].mean()
        output_text.append(f"MPL Score Analysis (higher is better)\n")
        output_text.append(f"{'-'*60}\n")
        output_text.append(f"Average MPL improvement: {avg_mpl_improvement:.2f}\n")
        output_text.append(f"(Positive = {method2} is better)\n\n")
        output_text.append("Detailed MPL Comparison:\n")
        output_text.append(mpl_comparison.to_string(index=False))
        output_text.append("\n")
    
    # Print to console
    result_text = "".join(output_text)
    print(result_text)
    
    # Save to file
    if output_path:
        with open(output_path, 'w') as f:
            f.write(result_text)
        print(f"\nResults saved to: {output_path}")


def analyze_networks(networks_csv, errors_csv, mpl_csv, method1, method2, output_path):
    """Main analysis function"""
    
    # Load data
    df_networks, df_errors, df_mpl = load_and_prepare_data(networks_csv, errors_csv, mpl_csv)
    
    # Find differences
    num_different, fraction_different, different_replicates = find_network_differences(df_networks, method1, method2)
    
    # Compare metrics
    error_comparison, mpl_comparison = compare_metrics(df_errors, df_mpl, method1, method2, different_replicates)
    
    # Print results
    print_results(num_different, fraction_different, error_comparison, mpl_comparison, method1, method2, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze network differences and their impact on error and MPL metrics",
        epilog=CONFIG_HELP,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    parser.add_argument(
        "-c",
        "--config",
        required=True,
        help="path to the TOML configuration file",
    )
    
    
    args = parser.parse_args()
    config = load_config(args.config, REQUIRED_KEYS)
    
    print(f"Configuration loaded: {config}\n")
    
    analyze_networks(
        config["networks_csv"],
        config["errors_csv"],
        config["mpl_csv"],
        config["method1"],
        config["method2"],
        config.get("output")
    )
