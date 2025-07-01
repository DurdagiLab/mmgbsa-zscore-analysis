# MM/GBSA Z-Score Analysis Tool
A Python-based utility for statistical post-processing of MM/GBSA binding free energy results in virtual screening (VS) studies. This script calculates Z-scores for MM/GBSA energies, filters ligands based on user-defined statistical thresholds, and generates a structured PDF summary and high-resolution plots suitable for publication or internal reporting.

# Description
The **MM/GBSA Z-Score Analysis Tool** evaluates the statistical significance of binding free energies (in CSV format) derived from MM/GBSA calculations. It is designed to assist in the prioritization of candidate compounds following MM/GBSA (re)scoring stages in VS workflows.

Key features:
- Robust parsing of numerical values, with support for locale-specific decimal formats (e.g., comma separators in European datasets).
- Computation of Z-scores for each compound based on the mean and standard deviation of the entire MM/GBSA energy distribution.
- Filtering of statistically significant hits using a customizable Z-score threshold (default: Z ≤ –1.960).
- Automated generation of result files, including:
  - A CSV file containing compounds passing the statistical threshold, along with their Z-scores;
  - A structured, tabular PDF report summarizing the filtered compounds;
  - A high-resolution plot of the standard normal distribution with the selected threshold indicated;
  - A distribution plot of MM/GBSA energies with the last processed molecule's score annotated.

# Suggested Z-score Thresholds
The script applies statistical filtering using two-tailed Z-score thresholds. The appropriate absolute Z-score (|Z|) can be selected based on the desired confidence level, as outlined below:

- For 50% confidence, use a threshold of |Z| ≥ 0.674
- For 75% confidence, use |Z| ≥ 1.150
- For 90% confidence, use |Z| ≥ 1.645
- For 95% confidence, use |Z| ≥ 1.960
- For 97% confidence, use |Z| ≥ 2.170
- For 99% confidence, use |Z| ≥ 2.576
- For 99.9% confidence, use |Z| ≥ 3.290

**These thresholds are based on the standard normal distribution for two-tailed significance testing and can be modified in the script to accommodate different confidence levels as required by the analysis.

# Requirements & Installation

Python version:
- `Python ≥ 3.6`

Dependencies:
- `numpy`  
- `matplotlib`  
- `scipy`  
- `fpdf`

Install via pip:

> pip install numpy matplotlib scipy fpdf

# Usage
Place the input CSV file (default name: FILE_NAME.csv) in the same directory as the script and run the following command in your terminal:

> python mmgbsa_zscore_analysis.py

**All output files, including the filtered CSV containing statistically significant compounds, the PDF summary table, and the high-resolution distribution plots, are saved automatically in the working directory when the script is run.

Note: The default Z-score threshold is set to –1.960, corresponding to a 95% confidence level. This value can be modified directly in the script to apply different statistical significance criteria depending on the analysis requirements.

# Citation
If you use this tool in your academic work, please cite:

Computational Drug Design Center (HITMER), Faculty of Pharmacy, Bahçeşehir University, Istanbul, Turkey
