# 必要なライブラリをまとめてインポート
import pandas as pd
from pyHSICLasso import HSICLasso
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import japanize_matplotlib  # Matplotlibで日本語を使用するためにインポート

# 1. データの読み込み
methylation_file = ""
efficiency_file = ""

try:
    df = pd.read_csv(methylation_file, index_col=0)
    Y = pd.read_csv(efficiency_file, header=None).values.flatten()
except FileNotFoundError as e:
    print(f"エラー: ファイルが見つかりません。パスを確認してください。")
    print(e)
    # スクリプトを続行できないため、ダミーデータを作成 (エラー回避のため)
    # df = pd.DataFrame(np.random.rand(100, 10), columns=[f'Sample{i+1}' for i in range(10)], index=[f'CpG{i+1}' for i in range(100)])
    # Y = np.random.rand(10)
    exit() # 本来はここで終了すべき


# 2. データの前処理
X = df.T.values
cg_names = df.index.to_list()
print(f"データ形状: X={X.shape}, Y={Y.shape}")

# 3. HSICLassoの実行
model = HSICLasso()
model.input(X, Y)
model.regression(num_feat=100, B=0) # 上位n個のCpGを選択 (ここはHSICLassoの選択数)
selected_idx = model.get_index()
selected_scores = model.get_index_score()

# 4. HSIClassoの結果を表示 
print("\n--- 4. HSIClassoによって選択されたCpGサイト（スコア） ---")
for i, (idx, score) in enumerate(zip(selected_idx, selected_scores)):
    print(f"順位{i+1}: {cg_names[idx]} (スコア = {score:.4f})")

# 5 & 6 重複を除去しつつグループ化してCSVに保存 
import pandas as pd
import numpy as np

print("\n--- 5. 重複を除去しつつ、選択されたCpGとその近傍特徴をグループ化 ---")

# Step 1: 準備
processed_cpgs = set()
grouped_dfs = []
blank_row_df = pd.DataFrame([[np.nan] * len(df.columns)], columns=df.columns, index=[' '])
df_T = pd.DataFrame(X.T, index=cg_names)

# num_neighbors = 100 # 固定数の設定
correlation_threshold = 0.7 # 相関係数のしきい値（絶対値）を設定

num_seeds = len(selected_idx)


# Step 2: ループ処理 
for i, seed_idx in enumerate(selected_idx):
    seed_name = cg_names[seed_idx]
    print(f"\rProcessing group {i + 1}/{num_seeds}: {seed_name}...", end="")

    if seed_name in processed_cpgs:
        print(f" -> スキップ ({seed_name} は既出のため)")
        continue # 次のシードへ

    # (A) 相関を計算
    seed_series = df_T.iloc[seed_idx]
    correlations = df_T.T.corrwith(seed_series)


    # 相関係数のしきい値(0.9)以上のものを抽出
    high_corr_series = correlations[correlations >= correlation_threshold]

    # シード自体（相関係数1.0）を除外し、CpG名のリストを取得
    # (dropにerrors='ignore'を指定し、万が一シードが含まれない場合もエラーにしない)
    all_neighbor_names = high_corr_series.drop(seed_name, errors='ignore').index.tolist()
    # --- 変更点 (ここまで) ---


    # 近傍CpGの中から、まだファイルに追加されていないものだけを抽出
    unique_neighbor_names = [name for name in all_neighbor_names if name not in processed_cpgs]

    # (B) 今回のグループを構成するCpG名リストを作成（シード + 未使用の近傍）
    group_names = [seed_name] + unique_neighbor_names

    # (C) データを抽出し、リストに追加
    # (dfにgroup_namesが存在するか確認)
    valid_group_names = [name for name in group_names if name in df.index]
    if not valid_group_names:
        print(f" -> スキップ ({seed_name} の近傍が見つからないか、すべて処理済み)")
        continue
        
    group_df = df.loc[valid_group_names]
    grouped_dfs.append(group_df)
    grouped_dfs.append(blank_row_df)

    # 今回グループに追加したCpGを「使用済み」としてsetに記録
    processed_cpgs.update(valid_group_names)

# Step 3 & 4 は変更なし 
print("\nグループの処理が完了しました。最終的なDataFrameを結合します。")

if grouped_dfs:
    # 最後の空行を削除
    final_grouped_df = pd.concat(grouped_dfs[:-1])
else:
    final_grouped_df = pd.DataFrame(columns=df.columns)
    print("有効なグループが見つかりませんでした。")

# 出力ファイル名（しきい値がわかるように変更）
output_file = "" 
final_grouped_df.to_csv(output_file, index=True, header=True)

print(f"\n重複除去済みのグループ化されたCpGデータを「{output_file}」に保存しました。")
print(f"最終的なユニークCpG数: {len(processed_cpgs)}個")
