import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import japanize_matplotlib
from collections import defaultdict

# 1. 設定
ANNOTATION_FILE_PATH = ""
HIT_LIST_FILE_PATH = ""
OUTPUT_FILE_PATH = ""

# --- 列名の設定 ---
ANNOTATION_ID_COL = 'Name'
ANNOTATION_CHR_COL = 'CHR'
ANNOTATION_POS_COL = 'MAPINFO'
HIT_LIST_ID_COL = 'cpgID'

# --- 解析パラメータ ---
WINDOW_SIZE = 1000000
STEP_SIZE = 100000
SIGNIFICANCE_LEVEL = 0.01

# 2. データの読み込みと前処理
print("--- データの読み込みと前処理を開始 ---")
try:
    all_cpgs_map = pd.read_csv(ANNOTATION_FILE_PATH, usecols=[ANNOTATION_ID_COL, ANNOTATION_CHR_COL, ANNOTATION_POS_COL], low_memory=False)
    all_cpgs_map.rename(columns={'Name': 'CpG_ID', 'CHR': 'Chromosome', 'MAPINFO': 'Start'}, inplace=True)
    hits_list = pd.read_csv(HIT_LIST_FILE_PATH)
    hits_list.rename(columns={HIT_LIST_ID_COL: 'CpG_ID'}, inplace=True)
except FileNotFoundError as e:
    print(f"エラー: ファイルが見つかりません。パスを確認してください。-> {e.name}")
    exit()

all_cpgs_map.dropna(inplace=True)
all_cpgs_map['Chromosome'] = 'chr' + all_cpgs_map['Chromosome'].astype(str)
valid_chroms = [f'chr{i}' for i in range(1, 23)]
all_cpgs_map = all_cpgs_map[all_cpgs_map['Chromosome'].isin(valid_chroms)]
all_cpgs_map['Start'] = pd.to_numeric(all_cpgs_map['Start'], errors='coerce').dropna().astype(int)
all_cpgs_map['End'] = all_cpgs_map['Start'] + 1

# 3. アノテーション
print("--- アノテーションを実行中 ---")
hits_with_coords = pd.merge(hits_list, all_cpgs_map, on='CpG_ID')
total_hits_count = len(hits_with_coords)
total_background_count = len(all_cpgs_map)
print(f"ヒットしたCpG数: {total_hits_count}")
print(f"全CpG数（背景）: {total_background_count}")

# 4. 濃縮解析（intervaltreeを使用）
print("--- 濃縮解析を実行中（intervaltreeを使用） ---")
# 染色体ごとにIntervalTreeを作成
hit_trees = defaultdict(IntervalTree)
background_trees = defaultdict(IntervalTree)

for _, row in hits_with_coords.iterrows():
    hit_trees[row['Chromosome']].addi(row['Start'], row['End'])

for _, row in all_cpgs_map.iterrows():
    background_trees[row['Chromosome']].addi(row['Start'], row['End'])

# スライディングウィンドウを作成
results_data = []
for chrom in sorted(all_cpgs_map['Chromosome'].unique(), key=lambda x: int(x.replace('chr',''))):
    chrom_max_pos = all_cpgs_map[all_cpgs_map['Chromosome'] == chrom]['Start'].max()
    for start in range(0, chrom_max_pos, STEP_SIZE):
        end = start + WINDOW_SIZE
        n_hits = len(hit_trees[chrom].overlap(start, end))
        n_background = len(background_trees[chrom].overlap(start, end))
        
        if n_background > 0: # 領域にCpGが1つもなければスキップ
            results_data.append([chrom, start, end, n_hits, n_background])

results_df = pd.DataFrame(results_data, columns=["Chromosome", "Start", "End", "n_hits", "n_background"])

# フィッシャーの正確確率検定
p_values = []
for _, row in results_df.iterrows():
    a = row['n_hits']
    b = total_hits_count - a
    c = row['n_background'] - a
    d = (total_background_count - total_hits_count) - c
    _, p_value = fisher_exact([[a,b], [c,d]], alternative='greater')
    p_values.append(p_value)
results_df['p_value'] = p_values

# 多重検定補正
_, q_values, _, _ = multipletests(results_df['p_value'].fillna(1.0), method='fdr_bh')
results_df['q_value'] = q_values

# 5. 結果の可視化と保存
print("--- 結果を可視化・保存中 ---")
results_df['log10_q'] = -np.log10(results_df['q_value'].replace(0, 1e-300))
plt.figure(figsize=(18, 7))
colors = ['#1f77b4', '#ff7f0e']
chrom_offsets = {}
current_offset = 0
unique_chroms = sorted(results_df['Chromosome'].unique(), key=lambda x: int(x.replace('chr','')))

for i, chrom in enumerate(unique_chroms):
    chrom_df = results_df[results_df['Chromosome'] == chrom].copy()
    chrom_df['plot_pos'] = chrom_df['Start'] + current_offset
    plt.scatter(chrom_df['plot_pos'], chrom_df['log10_q'], color=colors[i % len(colors)], s=15)
    
    chrom_max_pos = chrom_df['plot_pos'].max() if not chrom_df.empty else current_offset
    chrom_offsets[chrom] = current_offset + (chrom_max_pos - current_offset) / 2
    current_offset = chrom_max_pos

plt.title("CpGサイトの濃縮解析結果（マンハッタンプロット）")
plt.xlabel("染色体")
plt.ylabel("-log10(q-value)")
plt.xticks(list(chrom_offsets.values()), [c.replace('chr','') for c in chrom_offsets.keys()], rotation=45)
plt.axhline(y=-np.log10(SIGNIFICANCE_LEVEL), color='r', linestyle='--', label=f'q-value = {SIGNIFICANCE_LEVEL}')
plt.legend()
plt.tight_layout()
plt.show()

results_df.sort_values('q_value').to_csv(OUTPUT_FILE_PATH, index=False)
print(f"\n解析結果を '{OUTPUT_FILE_PATH}' に保存しました。")
print("\n--- 最も有意に濃縮していた領域 TOP 10 ---")
print(results_df.sort_values('q_value').head(10))