# Hamming Distance and PAM Clustering Analysis

本项目使用 R 对多序列比对后的可变位点矩阵进行模式分类，主要包括：

- 计算序列间的比例 Hamming 距离；
- 使用 PAM（Partitioning Around Medoids，围绕中心点划分）聚类；
- 比较不同聚类数的平均轮廓系数；
- 使用重复子抽样和调整 Rand 指数评估聚类稳定性；
- 计算各可变位点的熵降低；
- 使用置换检验和 Benjamini–Hochberg 方法进行 FDR 校正；
- 输出模式成员、中心序列、状态概率、统计表和图形。

## 项目结构

```text
Hcov/
├── Hcov.Rproj
├── hcov.R
├── README.md
├── .gitignore
├── data/
│   └── protein_sequence_analysis.xlsx
└── results/
    └── motif_analysis/
```

`results/` 中的文件由脚本运行后自动生成。

## 输入数据

默认输入文件：

```text
data/protein_sequence_analysis.xlsx
```

默认工作表：

```text
variable_matrix
```

输入表需满足：

1. 第一列或其中一列名为 `Sequence_ID`；
2. 每条序列的 `Sequence_ID` 唯一且非空；
3. 其余列为可变位点，列名格式为 `Pos_数字`，例如 `Pos_417`；
4. 单元格保存氨基酸单字母状态，也允许 `-`、`X`、`?` 和 `*`；
5. 不允许存在空值。

## 依赖

```r
install.packages(c("readxl", "openxlsx", "ggplot2", "cluster"))
```

## 运行方法

用 RStudio 打开 `Hcov.Rproj`，然后运行：

```r
source("hcov.R")
```

也可以在终端中运行：

```bash
Rscript hcov.R
```

## 分析方法

两条序列之间的比例 Hamming 距离定义为：

```text
不同状态的位点数 / 总可变位点数
```

脚本比较 `K = 2` 至 `K = 6` 的 PAM 聚类结果。在最小聚类样本数不少于 3 的候选结果中，选择平均轮廓系数最大的聚类数。

聚类稳定性通过 200 次、每次抽取 80% 样本的重复子抽样进行评估，并计算调整 Rand 指数。各位点对模式区分的贡献使用熵降低衡量，并通过 2000 次置换检验获得经验 P 值，随后进行 FDR 校正。

## 主要输出

运行结果保存在：

```text
results/motif_analysis/
```

主要文件包括：

- `hamming_distance_matrix.csv`：Hamming 距离矩阵；
- `pattern_number_selection.csv`：不同聚类数的评价指标；
- `pattern_membership.csv`：每条序列所属模式；
- `major_pattern_profiles.csv`：各模式的中心序列和多数状态；
- `pattern_defining_variables.csv`：模式区分位点及统计结果；
- `pattern_state_probabilities.csv`：各模式中每个位点状态的概率；
- `major_pattern_analysis.xlsx`：汇总工作簿；
- `pattern_number_selection.pdf`：聚类数选择图；
- `pattern_defining_variables.pdf`：主要区分位点图；
- `major_pattern_analysis.rds`：完整 R 分析对象。

## 可重复性参数

- 随机种子：`20260721`
- 子抽样次数：`200`
- 子抽样比例：`80%`
- 置换检验次数：`2000`

## 注意

本脚本将残基字符作为分类状态进行统计分析，不自动进行生物学功能解释。
