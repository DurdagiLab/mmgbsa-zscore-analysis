# MM/GBSA Z-Score Analysis Tool
A Python-based tool for statistical post-processing of MM/GBSA binding free-energy results in virtual screening studies.  
This script computes Z-scores for MM/GBSA energies and filters ligands based on user-defined significance thresholds.  
It also produces a summary PDF table and high-resolution visualizations suitable for publication or internal reporting.

# Description
The **MM/GBSA Z-Score Analysis Tool** processes binding free energy results (CSV format) from MM/GBSA calculations and evaluates the statistical significance of each compound's binding energy. The tool is particularly useful for prioritizing hits following molecular docking and MM/GBSA rescoring workflows in virtual screening campaigns. The workflow includes:
- Robust parsing of numeric fields, including locale-specific decimal separators (e.g., commas in European formats).
- Calculation of Z-scores for MM/GBSA binding energies based on the dataset mean and standard deviation.
- Filtering of statistically significant compounds using a user-defined threshold (default: Z ≤ -1.96).
- Automated output of:
  - A CSV file listing significant compounds and their Z-scores;
  - A structured PDF report summarizing the selected ligands;
  - A high-resolution plot of the standard normal distribution with the threshold marked;
  - An empirical MM/GBSA energy distribution plot highlighting the last molecule’s score.

# Suggested Z-score Thresholds
The script supports statistical filtering based on two-tailed Z-score significance levels. Depending on the desired confidence level, the corresponding absolute Z-score threshold (|Z|) can be selected as follows:

- For 50% confidence, use a threshold of |Z| ≥ 0.674
- For 75% confidence, use |Z| ≥ 1.150
- For 90% confidence, use |Z| ≥ 1.645
- For 95% confidence, use |Z| ≥ 1.960
- For 97% confidence, use |Z| ≥ 2.170
- For 99% confidence, use |Z| ≥ 2.576
- For 99.9% confidence, use |Z| ≥ 3.290

**These values correspond to the standard normal distribution under two-tailed statistical testing and can be adjusted in the script to meet specific confidence requirements.

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

> python docking-zscore-analysis.py

**All generated artefacts - including the filtered CSV file, summary PDF report, and high-resolution PNG figures - will be automatically saved to the current working directory.

Note: The script uses a default Z-score threshold of -1.960 to select statistically significant compounds. You may edit the script to adjust this threshold as needed.

# Citation
If you use this tool in your academic work, please cite:

Computational Drug Design Center (HITMER), Faculty of Pharmacy, Bahçeşehir University, Istanbul, Turkey
