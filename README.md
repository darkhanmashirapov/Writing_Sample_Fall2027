# Geopolitical Shock Transmission in Emerging Bond Markets

Data, code, and analysis pipeline for a writing sample examining whether
the Russia-Ukraine war altered the relationship between corporate
credit-risk spreads and macro-financial conditions in Kazakhstan and
Russia.

## Structure
- `notebooks/` — seven Jupyter notebooks, run in numeric order, covering
  data collection, cleaning, exploratory analysis, ARDL bounds testing,
  Markov-switching regime detection, and VAR/IRF analysis.
- `data/master_dataset_final.xlsx` — the final merged dataset (Levels,
  Transformed, and Codebook sheets).
- `thesis/` — the full written thesis (LaTeX source).

## Requirements
Python 3.x with pandas, numpy, statsmodels, matplotlib, seaborn, openpyxl.

## Reproducing results
Run notebooks 01 through 06 in order; each depends on outputs saved by
the previous one.
