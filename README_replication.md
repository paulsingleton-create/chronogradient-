# Chronogradient Replication Package

## Paul Singleton, Nottingham, May 2026

### What this is

This package contains everything needed to independently verify the core findings
of the chronogradient hypothesis paper: "An Isotropic Velocity Multiplier in
925,864 Stars and Its Implications for Dark Matter."

### Contents

- `chronogradient_paper_may2026.docx` — The full paper (8 sections)
- `chronogradient_replication.py` — Self-contained Python script that downloads
  the data and reproduces all three key findings
- `gamma_extraction_for_GPT.txt` — Step-by-step method with intermediate values
- `dark_matter_crib_sheet.docx` — Summary scorecard of dark matter evidence addressed
- `README_replication.md` — This file

### Quick start

```
pip install pandas numpy scipy
python chronogradient_replication.py
```

The script will:
1. Download 651,208 stars from GALAH DR4 via public TAP query (~80 MB, ~5 min)
2. Apply metallicity control ([Fe/H] -0.3 to +0.3)
3. Run three tests:
   - Percentile uniformity (expected: ~3.3 pp spread)
   - Mass-dependent sigmoid (expected: R² 0.75 to 0.98)
   - Isotropic gamma extraction (expected: 0.2% spread at 6.5 Gyr)

### Data sources (all public)

- GALAH DR4: https://datacentral.org.au (TAP query)
- Gaia DR3: https://gea.esac.esa.int
- APOGEE DR17: VizieR catalogue J/A+A/673/A155

### Key findings to verify

1. **Percentile uniformity**: P5 through P90 rise by ~51% with 3.32 pp spread
2. **Sigmoid threshold**: drops from 10.6 Gyr at 0.875 Msun to 7.8 Gyr at 1.025 Msun
3. **Isotropic gamma**: at 6-7 Gyr, gamma_U = 1.095, gamma_V = 1.098, gamma_W = 1.095

### Falsifiable predictions (untested)

1. Stellar velocity dispersions in dwarf spheroidals, corrected for chronogradient,
   should reduce inferred dark matter by 50-80%
2. In the same galaxy, DM from old stellar dispersions should exceed DM from gas
   rotation curves. Gap should correlate with stellar population age.
3. Void galaxies measured with stellar kinematics should show less apparent DM
   than cluster galaxies at equal stellar mass.
4. Wide binary anomaly transition radius scales as ~M², not M^0.5 (MOND).

### Contact

Paul Singleton, Nottingham, England
Reddit: r/EmergentAIPersonas, r/FractalTapestry, r/CoherencePhysics
