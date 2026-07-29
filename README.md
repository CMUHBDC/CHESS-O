# CHESS-O Analysis Code

This repository contains the analysis code supporting the manuscript:

**Artificial Intelligence–Based Chest Radiography Screening for Osteoporosis: Model Development, External Validation, Fracture Risk Prediction, and Lead Time Evaluation**

## Overview

Osteoporosis is a global health challenge that contributes substantially to fractures, disability, healthcare costs, and mortality, yet remains underdiagnosed despite effective therapies. We developed Chest X-ray to Early Systemic Disease Screening for Osteoporosis (CHESS-O), a deep learning model for osteoporosis detection from routine chest radiographs. CHESS-O was trained using 16,172 radiographs from Taiwan and validated in 10,794 radiographs from Taiwan, the United Arab Emirates, and the United States, demonstrating robust discrimination across cohorts (area under the receiver operating characteristic curve: 0.862–0.952). In 76,051 individuals aged 50–55 years, CHESS-O–predicted osteoporosis was associated with elevated nontraumatic fracture risk (adjusted hazard ratio: 4.38 at 1 year; 2.46 overall). Among patients with DEXA-confirmed osteoporosis, CHESS-O detected osteoporosis a median of 3.6 years before diagnosis. These findings support AI–enabled opportunistic chest radiography screening for osteoporosis risk stratification and early intervention. Prospective validation and workflow integration are needed to assess real-world impact.

## Repository contents

| File | Description |
|---|---|
| `fig2a_calibration_plot.R` | Code used to generate the Figure 2A. |
| `supp_table4_code.ipynb` | Code used to generate Supplementary Table 4. |
| `table2_code.ipynb` | Code used to generate Table 2. |
| `table3_fracture_risk_analysis.R` | Code used to generate Table 3. |

## Data availability

Individual-level study data are not included because of institutional data-governance and participant-privacy restrictions. The code documents the analytical procedures used in the study. Users must supply appropriately structured analytic data to reproduce the analyses.

## Reproducibility

The public code does not include individual-level data. Where a script refers to an object such as `analytic_data`, that object is a placeholder for the corresponding restricted analytic dataset.
