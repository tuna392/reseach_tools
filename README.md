[reseach_tools-README.md](https://github.com/user-attachments/files/28135362/reseach_tools-README.md)
# research_tools

バイオインフォマティクス研究向けのデータ解析スクリプト集です。遺伝子発現データの統計解析・エンリッチメント解析を自動化します。

## 概要

分子生物学・バイオインフォマティクス研究で繰り返し行う解析処理をPythonスクリプトとして整備したツール集です。研究データの前処理から統計解析・可視化までをパイプライン化し、再現性の高い解析環境を提供します。

## スクリプト一覧

### `Enrichment_Analysis.py`

取得したプローブののエンリッチメント解析用。

### `HSICLasso_Cor0.9.py`

HSIC Lasso（Hilbert-Schmidt Independence Criterion + Lasso）による非線形特徴選択スクリプトです。高次元の遺伝子発現データから、目的変数と非線形な依存関係を持つ重要な遺伝子特徴量を選択します。相関係数0.9以上の特徴量を事前フィルタリングし、解析の精度と計算効率を両立しています。

## 技術スタック

| 領域 | 技術・ライブラリ |
|------|------|
| 言語 | Python 3.x |
| 統計・機械学習 | pyHSICLasso, scikit-learn |
| データ処理 | pandas, numpy |
| 可視化 | matplotlib, seaborn |
| エンリッチメント解析 | gseapy / biomaRt連携 |

## 背景・動機

バイオインフォマティクス解析において、遺伝子発現データは特徴量（遺伝子）数がサンプル数を大幅に上回る高次元問題（p >> n）です。従来の線形相関ベースの特徴選択では捉えられない非線形な遺伝子間相互作用を考慮するため、カーネル法ベースのHSIC Lassoを採用しています。

## 使い方

```bash
# 依存ライブラリのインストール
pip install -r requirements.txt

# エンリッチメント解析の実行
python Enrichment_Analysis.py --input gene_list.csv --output results/

# HSIC Lassoによる特徴選択
python HSICLasso_Cor0.9.py --input expression_matrix.csv --target phenotype.csv
```

## 入力データ形式

`HSICLasso_Cor0.9.py` は以下の形式のCSVを想定しています：

```
# expression_matrix.csv
gene_id, sample_1, sample_2, ...
GENE_A,  1.23,     4.56, ...

# phenotype.csv
sample_id, label
sample_1,  0
sample_2,  1
```
