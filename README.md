# CHESS-O Analysis Code

This repository contains the analysis code supporting the manuscript:

**Fracture Risk and Lead Time With AI-Based Osteoporosis Screening on Routine Chest Radiographs**

## Overview

**Background.** Artificial intelligence (AI) applied to routine chest radiographs (CXR) could identify osteoporosis, but the detection window before diagnosis or fracture remains under-characterized. We evaluated CHESS-O, a CXR-based AI, for early osteoporosis identification and fracture-risk association. 

**Methods.** CHESS-O was developed using 16,172 posteroanterior CXRs paired with dual-energy x-ray absorptiometry (DXA; osteoporosis: T-score ≤−2.5), evaluated in 3,794 held-out radiographs, externally validated in 7,000 radiographs from Taiwan, the United Arab Emirates, and the United States, and compared with Osteoporosis Self-Assessment Tool for Asians index (OSTAi). Lead time from the earliest observed CHESS-O–positive CXR to DXA-confirmed osteoporosis or fracture was evaluated in separate cohorts with prior CXRs. Associations between CHESS-O estimated osteoporosis and incident nontraumatic major osteoporotic fractures were assessed in 76,051 adults aged 50–55 years using cause-specific Cox models accounting for competing death, presented as adjusted hazard ratios (aHRs) and absolute risk increases (aARIs). 

**Findings.** Across international cohorts, areas under the receiver operating characteristic curve ranged from 0.862 (95% CI 0.819–0.901) to 0.952 (0.946–0.957), compared with 0.832 (0.817–0.848) for OSTAi in internal held-out cohort. Among 286 patients with DXA-confirmed osteoporosis and prior CXRs, 68.5% tested positive beforehand (median lead time 5.6 years \[IQR 2.11–8.37]). Among 1,533 patients with fracture and prior CXRs, 23.4% tested positive beforehand (lead time 4.3 years \[IQR 1.98–6.61]). CHESS-O–predicted osteoporosis was associated with fracture within 1 year (aHR: 2.73 \[95% CI 2.00–3.74]; aARI: 7.42 \[95% CI 4.33–10.9] per 1,000 persons) to 5 years (2.83 \[2.35–3.40]; 31.6 \[23.1–39.6]). 

**Interpretation.** CHESS-O identified osteoporosis years before DXA-confirmed diagnosis or fracture and was associated with subsequent fracture risk, with consistent discrimination across international cohorts. Prospective evaluation should determine whether CHESS-O–guided screening improves care and reduces fractures.

## Repository contents

|File|Description|
|-|-|
|`fig2\_calibration\_plot.R`|Code used to generate the Figure 2.|
|`supp\_table5\_code.ipynb`|Code used to generate Supplementary Table 5.|
|`table2\_code.ipynb`|Code used to generate Table 2.|
|`table3\_fracture\_risk\_analysis.R`|Code used to generate Table 3.|

## Data availability

Individual-level study data are not included because of institutional data-governance and participant-privacy restrictions. The code documents the analytical procedures used in the study. Users must supply appropriately structured analytic data to reproduce the analyses.

## Reproducibility

The public code does not include individual-level data. Where a script refers to an object such as `analytic\_data`, that object is a placeholder for the corresponding restricted analytic dataset.

