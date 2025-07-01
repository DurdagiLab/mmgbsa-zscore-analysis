"""
###############################################################################
Title: MM/GBSA Z-Score Analysis Tool

Developed by: Mine Isaoglu, Ph.D.
Principal Investigator: Serdar Durdagi, Ph.D.
Affiliation: Computational Drug Design Center (HITMER), Faculty of Pharmacy,
             Bahçeşehir University, Istanbul, Turkey
             
Version: January 2025
################################################################################
"""

import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from fpdf import FPDF


def parse_float(value):
    """Convert numeric strings with commas or dots to float.

    Handles Turkish/European decimal commas gracefully.
    Returns *None* for non-numeric values.
    """
    try:
        # If the value contains a comma and has at least two characters after the comma
        if "," in value and len(value.split(";")[1]) >= 2:
            return float(value.replace(",", "."))
        # If the value contains a comma but has less than two characters after the comma
        elif "," in value:
            return float(value.replace(",", ".") + "00")  # Pad missing decimals
        else:
            return float(value)
    except ValueError:
        return None  # Return None if the value cannot be converted to float


def z_score(data, x):
    """Return the Z-score of *x* with respect to *data*."""
    mean = np.mean(data)
    std_dev = np.std(data)
    return (x - mean) / std_dev


def print_last_molecule_mmgbsa(csv_file_path):
    """Print and return the MMGBSA score of the last molecule in *csv_file_path*."""
    with open(csv_file_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        data = [(row["Title"], float(row["MMGBSA dG Bind"])) for row in reader]

    molecule_name, mmgbsa_score = data[-1]
    print(f"The MMGBSA score of the last molecule ({molecule_name}) is: {mmgbsa_score}")
    return float(mmgbsa_score)


def main():
    # -----------------------------------------------------------------------
    # Configuration
    # -----------------------------------------------------------------------
    file_name = "your_file.csv"
    updated_file_name = "your_file_with_Z_Scores.csv"
    output_pdf_name = "Z-Score_Table.pdf"
    output_image_name = "Z_Score_Distribution_Curve.png"
    mmgbsa_output_image_name = "MMGBSA_Distribution_Curve.png"

    working_directory = os.getcwd()

    file_path = os.path.join(working_directory, file_name)
    updated_file_path = os.path.join(working_directory, updated_file_name)
    output_pdf_path = os.path.join(working_directory, output_pdf_name)
    output_image_path = os.path.join(working_directory, output_image_name)
    mmgbsa_output_image_path = os.path.join(working_directory, mmgbsa_output_image_name)

    molecule_name_column = "Title"
    dg_bind_column = "MMGBSA dG Bind"
    z_score_threshold_lower = -1.960  # Lower bound for Z-score filtering

    # -----------------------------------------------------------------------
    # Data ingestion and preprocessing
    # -----------------------------------------------------------------------
    with open(file_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        data = [(row[molecule_name_column], parse_float(row[dg_bind_column])) for row in reader]

    # Filter out non-numeric entries
    data = [(molecule, value) for molecule, value in data if value is not None]

    # -----------------------------------------------------------------------
    # Z-score calculation
    # -----------------------------------------------------------------------
    dg_bind_values = np.array([item[1] for item in data])
    z_scores = [(x - np.mean(dg_bind_values)) / np.std(dg_bind_values) for x in dg_bind_values]

    selected_data = [
        (item[0], item[1], f"{z:.6f}")
        for item, z in zip(data, z_scores) if z <= z_score_threshold_lower
    ]

    # Append Z-scores to full dataset (for later use in plotting)
    data_with_z = [
        (item[0], item[1], z)
        for item, z in zip(data, z_scores)
    ]

    # -----------------------------------------------------------------------
    # Output: CSV containing only significant hits
    # -----------------------------------------------------------------------
    with open(updated_file_path, "w", newline="") as csvfile:
        fieldnames = [molecule_name_column, dg_bind_column, "Z Score"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for item in selected_data:
            writer.writerow({
                molecule_name_column: item[0],
                dg_bind_column: item[1],
                "Z Score": item[2],
            })

    # -----------------------------------------------------------------------
    # Output: PDF report
    # -----------------------------------------------------------------------
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", style="B", size=12)

    pdf.cell(200, 10, txt="Compounds with Z-Score <= -1.960", ln=True, align="C")
    pdf.ln(10)
    pdf.cell(200, 10, txt=f"Total Number of Selected Compounds: {len(selected_data)}", ln=True, align="C")
    pdf.ln(10)

    pdf.set_font("Arial", style="B", size=8)
    pdf.cell(60, 10, "Compound ID", border=1)
    pdf.cell(60, 10, "MMGBSA dG Bind (kcal/mol)", border=1)
    pdf.cell(60, 10, "Calculated Z-Score", border=1)
    pdf.ln(10)

    pdf.set_font("Arial", size=10)
    for item in selected_data:
        pdf.cell(60, 10, item[0], border=1)
        pdf.cell(60, 10, str(item[1]), border=1)
        pdf.cell(60, 10, item[2], border=1)
        pdf.ln(10)

    pdf.output(output_pdf_path)

    # -----------------------------------------------------------------------
    # Output: Z-score distribution plot
    # -----------------------------------------------------------------------
    plt.figure(figsize=(8, 6))
    plt.title("Z-Score Normal Distribution", fontweight="bold")
    plt.xlabel("Z-Score")
    plt.ylabel("Probability Density")

    x = np.linspace(-5, 5, 100)
    plt.plot(x, norm.pdf(x, 0, 1), label="Z-Score Normal Distribution", linewidth=2.5)
    plt.axvline(z_score_threshold_lower, linestyle="--")
    plt.axvline(0, linestyle="-")
    plt.legend()

    plt.savefig(output_image_path, bbox_inches="tight", dpi=300)
    plt.close()

    # -----------------------------------------------------------------------
    # Output: MM/GBSA energy distribution plot
    # -----------------------------------------------------------------------
    plt.figure(figsize=(8, 6))
    mmgbsa_values = [item[1] for item in data_with_z]
    plt.title("MM/GBSA Binding Free-Energy Distribution", fontweight="bold")
    plt.xlabel("MM/GBSA dG Bind Energy (kcal/mol)")
    plt.ylabel("Probability Density")

    mmgbsa_mean = np.mean(mmgbsa_values)
    mmgbsa_std_dev = np.std(mmgbsa_values)
    mmgbsa_x = np.linspace(min(mmgbsa_values), max(mmgbsa_values), 100)
    plt.plot(mmgbsa_x, norm.pdf(mmgbsa_x, mmgbsa_mean, mmgbsa_std_dev), label="MM/GBSA Energy Distribution", linewidth=2.5)

    last_mmgbsa_value = print_last_molecule_mmgbsa(updated_file_path)
    plt.axvline(last_mmgbsa_value, linestyle="--")
    plt.text(last_mmgbsa_value - 1, 0.01, f"Fitness = {last_mmgbsa_value:.2f}", fontsize=12, fontweight="bold", ha="right")

    plt.legend()
    plt.savefig(mmgbsa_output_image_path, bbox_inches="tight", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
