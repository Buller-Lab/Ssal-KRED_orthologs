import argparse
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, skew, kurtosis

parser = argparse.ArgumentParser()
parser.add_argument('--orthologs', required=True)
parser.add_argument('--swissprot', required=True)
parser.add_argument('--max_seqid', type=float, default=0.3)
args = parser.parse_args()

ortho = pd.read_csv(args.orthologs)
swiss = pd.read_csv(args.swissprot)

for df in [ortho, swiss]:
    if 'sequence_identity' in df.columns:
        df['sequence_identity'] = pd.to_numeric(
            df['sequence_identity']
            .astype(str)
            .str.replace(r'[%\, ]', '', regex=True)
            .replace(['', 'nan', 'NaN', 'NA', 'n.d.', 'N.D.', 'nd', 'ND'], np.nan),
            errors='coerce'
        )

ortho = ortho[ortho['sequence_identity'].le(args.max_seqid)]
swiss = swiss[swiss['sequence_identity'].le(args.max_seqid)]

def is_completely_missing(col_series):
    missing_values = {np.nan, None, '', 'nan', 'NaN', 'NA', 'n.d.', 'N.D.', 'nd', 'ND'}
    return col_series.astype(str).str.strip().isin(missing_values).all() or col_series.isna().all()

for df in [ortho, swiss]:
    cols_to_drop = [col for col in df.columns if is_completely_missing(df[col])]
    df.drop(columns=cols_to_drop, inplace=True)

common_numeric = [
    col for col in ortho.columns
    if col in swiss.columns
    and pd.api.types.is_numeric_dtype(ortho[col])
    and pd.api.types.is_numeric_dtype(swiss[col])
]

t_test_results = []
desc_results = []

for col in common_numeric:
    d1 = ortho[col].dropna()
    d2 = swiss[col].dropna()
    
    n1 = len(d1)
    n2 = len(d2)
    
    if n1 >= 2 and n2 >= 2:
        t, p = ttest_ind(d1, d2, equal_var=False)
        t_test_results.append({
            'column': col,
            'mean_orthologs': round(d1.mean(), 4),
            'mean_swissprot': round(d2.mean(), 4),
            't_stat': round(t, 4),
            'p_value': round(p, 6),
            'significant_0.05': p < 0.05,
            'n_orthologs': n1,
            'n_swissprot': n2
        })
    
    if n1 > 0:
        desc_results.append({
            'column': col,
            'group': 'orthologs',
            'count': n1,
            'mean': round(d1.mean(), 4),
            'median': round(d1.median(), 4) if n1 > 1 else np.nan,
            'std': round(d1.std(), 4) if n1 > 1 else np.nan,
            'var': round(d1.var(), 4) if n1 > 1 else np.nan,
            'min': round(d1.min(), 4),
            'max': round(d1.max(), 4),
            'skewness': round(skew(d1), 4) if n1 > 2 else np.nan,
            'kurtosis': round(kurtosis(d1), 4) if n1 > 3 else np.nan
        })
    
    if n2 > 0:
        desc_results.append({
            'column': col,
            'group': 'swissprot',
            'count': n2,
            'mean': round(d2.mean(), 4),
            'median': round(d2.median(), 4) if n2 > 1 else np.nan,
            'std': round(d2.std(), 4) if n2 > 1 else np.nan,
            'var': round(d2.var(), 4) if n2 > 1 else np.nan,
            'min': round(d2.min(), 4),
            'max': round(d2.max(), 4),
            'skewness': round(skew(d2), 4) if n2 > 2 else np.nan,
            'kurtosis': round(kurtosis(d2), 4) if n2 > 3 else np.nan
        })

df_ttest = pd.DataFrame(t_test_results).sort_values('p_value')
df_desc  = pd.DataFrame(desc_results).sort_values(['column', 'group'])

with pd.ExcelWriter('t_test_results.xlsx', engine='openpyxl') as writer:
    df_ttest.to_excel(writer, sheet_name='t_tests', index=False)
    df_desc.to_excel(writer, sheet_name='descriptives', index=False)
    pd.DataFrame({
        'info': ['max_seqid_cutoff', 'orthologs_rows_after_filter', 'swissprot_rows_after_filter'],
        'value': [args.max_seqid, len(ortho), len(swiss)]
    }).to_excel(writer, sheet_name='info', index=False)