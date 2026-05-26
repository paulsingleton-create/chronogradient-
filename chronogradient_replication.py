"""
Chronogradient Replication Script
Paul Singleton, Nottingham, May 2026

This script reproduces the core findings from the chronogradient paper:
1. Percentile uniformity (3.32 pp spread)
2. Mass-dependent sigmoid
3. Isotropic gamma extraction (0.2% at onset)

REQUIREMENTS:
- Python 3.8+
- pandas, numpy, scipy
- Internet connection (to download GALAH DR4 via TAP)

The script downloads the data directly from the public GALAH DR4 archive.
No pre-processed files needed. Total download: ~80 MB.
Runtime: approximately 10 minutes including download.
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import urllib.request
import urllib.parse
import time
import os

# ============================================================
# STEP 1: DOWNLOAD GALAH DR4 DATA
# ============================================================

def download_galah():
    """Download GALAH DR4 stars via TAP query from datacentral.org.au"""
    
    tap_url = "https://datacentral.org.au/vo/tap/sync"
    
    mass_ranges = [
        (0.3, 0.6), (0.6, 0.8), (0.8, 0.9), (0.9, 0.95),
        (0.95, 1.0), (1.0, 1.05), (1.05, 1.1), (1.1, 1.2),
        (1.2, 1.5), (1.5, 2.0)
    ]
    
    all_data = []
    header = None
    
    for mlo, mhi in mass_ranges:
        query = f"""SELECT 
            m.sobject_id, m.mass, m.age, m.teff, m.logg, m.fe_h,
            m.rv_comp_1,
            d.U_UVW, d.V_UVW, d.W_UVW,
            d.R_Rzphi, d.z_Rzphi
        FROM galah_dr4.mainstartable AS m
        JOIN galah_dr4.vacdynamicstable AS d ON m.sobject_id = d.sobject_id
        WHERE m.age > 0 AND m.mass >= {mlo} AND m.mass < {mhi}
        AND m.rv_comp_1 IS NOT NULL
        AND m.flag_sp = 0
        AND d.U_UVW IS NOT NULL"""
        
        params = {'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': query}
        url = tap_url + '?' + urllib.parse.urlencode(params)
        
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                data = r.read().decode()
                lines = data.strip().split('\n')
                if header is None:
                    header = lines[0]
                rows = lines[1:]
                all_data.extend(rows)
                print(f"  Mass {mlo}-{mhi}: {len(rows)} stars")
        except Exception as e:
            print(f"  Mass {mlo}-{mhi}: FAILED - {e}")
        
        time.sleep(1)
    
    outpath = 'galah_dr4_stars.csv'
    with open(outpath, 'w') as f:
        f.write(header + '\n')
        for row in all_data:
            f.write(row + '\n')
    
    print(f"\nTotal: {len(all_data)} stars written to {outpath}")
    return outpath


# ============================================================
# STEP 2: PERCENTILE UNIFORMITY
# ============================================================

def test_percentile_uniformity(df):
    """Test whether all percentiles rise by the same amount"""
    
    print("\n" + "="*70)
    print("TEST 1: PERCENTILE UNIFORMITY")
    print("Solar core (0.95-1.05 Msun), [Fe/H] -0.3 to +0.3")
    print("="*70)
    
    solar = df[(df['mass'] >= 0.95) & (df['mass'] < 1.05)]
    young = solar[solar['age'] < 4]
    old = solar[solar['age'] > 8]
    
    print(f"\nYoung (<4 Gyr): {len(young)} stars")
    print(f"Old (>8 Gyr): {len(old)} stars")
    
    print(f"\n  NOTE: Using abs(radial_velocity) = abs(rv_comp_1)")
    print(f"  This is the line-of-sight component only, not total space velocity.")
    print(f"  The percentile uniformity holds in both RV and total velocity.\n")
    print(f"{'Percentile':<12} {'Young':>10} {'Old':>10} {'Rise':>10}")
    print("-"*45)
    
    rises = []
    for pct in [5, 10, 25, 50, 75, 90]:
        py = np.percentile(young['absrv'], pct)
        po = np.percentile(old['absrv'], pct)
        rise = (po - py) / py * 100
        rises.append(rise)
        print(f"  P{pct:<4}      {py:>10.2f} {po:>10.2f} {rise:>9.1f}%")
    
    spread = max(rises) - min(rises)
    print(f"\n  SPREAD: {spread:.2f} percentage points")
    print(f"  (Standard disc heating would produce 15-30+ pp spread)")
    
    return spread


# ============================================================
# STEP 3: MASS-DEPENDENT SIGMOID
# ============================================================

def test_sigmoid(df):
    """Fit sigmoids to age-velocity curves at different masses"""
    
    print("\n" + "="*70)
    print("TEST 2: MASS-DEPENDENT SIGMOID")
    print("Sigmoid: v = v0 + L / (1 + exp(-k * (t - t0)))")
    print("="*70)
    
    def sigmoid(t, v0, L, t0, k):
        return v0 + L / (1 + np.exp(-k * (t - t0)))
    
    mass_groups = [
        ('0.85-0.90', 0.85, 0.90),
        ('0.90-0.95', 0.90, 0.95),
        ('0.95-1.00', 0.95, 1.00),
        ('1.00-1.05', 1.00, 1.05),
    ]
    
    for mname, mlo, mhi in mass_groups:
        gl = df[(df['mass'] >= mlo) & (df['mass'] < mhi)]
        
        p50_values = []
        age_points = []
        
        for alo in range(0, 13):
            ag = gl[(gl['age'] >= alo) & (gl['age'] < alo + 1)]
            if len(ag) >= 50:
                p50_values.append(np.percentile(ag['absrv'], 50))
                age_points.append(alo + 0.5)
        
        if len(age_points) >= 6:
            try:
                popt, _ = curve_fit(sigmoid, np.array(age_points), np.array(p50_values),
                                    p0=[16.0, 10.0, 7.0, 0.5],
                                    bounds=([5, 0, 0, 0.01], [25, 50, 14, 5]),
                                    maxfev=5000)
                v0, L, t0, k = popt
                
                residuals = np.array(p50_values) - sigmoid(np.array(age_points), *popt)
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((np.array(p50_values) - np.mean(p50_values))**2)
                r2 = 1 - ss_res / ss_tot
                
                print(f"\n  {mname} Msun: t0={t0:.1f} Gyr, L={L:.1f} km/s, R2={r2:.3f}")
            except Exception as e:
                print(f"\n  {mname}: fit failed - {e}")


# ============================================================
# STEP 4: ISOTROPIC GAMMA EXTRACTION
# ============================================================

def test_gamma(df):
    """Extract the isotropic gamma multiplier"""
    
    print("\n" + "="*70)
    print("TEST 3: ISOTROPIC GAMMA EXTRACTION")
    print("Method: linear fit to pre-threshold (2-6 Gyr), divide out")
    print("="*70)
    
    # NOTE: Gamma extraction uses tighter metallicity cut (±0.2)
    # than the percentile test (±0.3). This matches the original analysis.
    solar = df[(df['mass'] >= 0.95) & (df['mass'] < 1.05) & 
               (df['fe_h'] >= -0.2) & (df['fe_h'] <= 0.2)]
    
    # Compute dispersions for each age bin
    results = {}
    for alo in range(1, 13):
        ag = solar[(solar['age'] >= alo) & (solar['age'] < alo + 1)]
        if len(ag) >= 100:
            su = ag['U_UVW'].std()
            sv = ag['V_UVW'].std()
            sw = ag['W_UVW'].std()
            results[alo + 0.5] = (su, sv, sw, len(ag))
    
    # Print raw dispersions
    print(f"\n{'Age mid':>8} {'sigma_U':>10} {'sigma_V':>10} {'sigma_W':>10} {'n':>8}")
    print("-"*50)
    for t in sorted(results.keys()):
        su, sv, sw, n = results[t]
        print(f"{t:>8.1f} {su:>10.1f} {sv:>10.1f} {sw:>10.1f} {n:>8}")
    
    # Linear fit to pre-threshold (2.5, 3.5, 4.5, 5.5)
    pre_ages = np.array([2.5, 3.5, 4.5, 5.5])
    pre_su = np.array([results[t][0] for t in pre_ages])
    pre_sv = np.array([results[t][1] for t in pre_ages])
    pre_sw = np.array([results[t][2] for t in pre_ages])
    
    u_fit = np.polyfit(pre_ages, pre_su, 1)
    v_fit = np.polyfit(pre_ages, pre_sv, 1)
    w_fit = np.polyfit(pre_ages, pre_sw, 1)
    
    print(f"\nLinear fits (slope, intercept):")
    print(f"  sigma_U = {u_fit[0]:.4f} * t + {u_fit[1]:.4f}")
    print(f"  sigma_V = {v_fit[0]:.4f} * t + {v_fit[1]:.4f}")
    print(f"  sigma_W = {w_fit[0]:.4f} * t + {w_fit[1]:.4f}")
    
    # Compute gamma
    print(f"\n{'Age':>6} {'gam_U':>8} {'gam_V':>8} {'gam_W':>8} {'mean':>8} {'spread':>8}")
    print("-"*50)
    for t in sorted(results.keys()):
        su, sv, sw, n = results[t]
        pu = u_fit[0] * t + u_fit[1]
        pv = v_fit[0] * t + v_fit[1]
        pw = w_fit[0] * t + w_fit[1]
        gu = su / pu
        gv = sv / pv
        gw = sw / pw
        gm = np.mean([gu, gv, gw])
        spread = (max([gu, gv, gw]) - min([gu, gv, gw])) / gm * 100
        marker = "  *** ONSET ***" if 6.0 <= t <= 7.0 else ""
        print(f"{t:>6.1f} {gu:>8.3f} {gv:>8.3f} {gw:>8.3f} {gm:>8.3f} {spread:>7.1f}%{marker}")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("CHRONOGRADIENT REPLICATION SCRIPT")
    print("Paul Singleton, Nottingham, May 2026")
    print("="*70)
    
    # Check if data already downloaded
    datafile = 'galah_dr4_stars.csv'
    if not os.path.exists(datafile):
        print("\nDownloading GALAH DR4 data...")
        datafile = download_galah()
    else:
        print(f"\nUsing existing {datafile}")
    
    # Load data
    print("\nLoading data...")
    df = pd.read_csv(datafile)
    print(f"Loaded: {len(df)} stars")
    
    # Apply metallicity control and quality flags
    # NOTE: flag_sp = 0 is applied in the TAP query at download time.
    # If using pre-downloaded data, ensure this filter was applied.
    df['absrv'] = abs(df['rv_comp_1'])
    df_met = df[(df['fe_h'] >= -0.3) & (df['fe_h'] <= 0.3) & (df['absrv'] < 150)]
    print(f"After metallicity control: {len(df_met)} stars")
    
    # Run tests
    spread = test_percentile_uniformity(df_met)
    test_sigmoid(df_met)
    test_gamma(df_met)
    
    print("\n" + "="*70)
    print("REPLICATION COMPLETE")
    print(f"Percentile spread: {spread:.2f} pp (expected: ~3.3)")
    print(f"Check gamma at 6.5 Gyr: should be ~1.096 in all three components")
    print(f"Check spread at 6.5 Gyr: should be ~0.2%")
    print("="*70)
