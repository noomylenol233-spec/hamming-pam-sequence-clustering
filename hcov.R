# ============================================================
# 多序列比对后的 Hamming 距离与聚类分析
# 说明：本脚本只把残基字符视为分类状态，不进行任何生物学解释。
# ============================================================

# ---------- 1. 集中设置文件路径 ----------
project_dir <- getwd()
input_excel <- file.path(project_dir, "data", "protein_sequence_analysis.xlsx")
output_dir <- file.path(project_dir, "results", "motif_analysis")
input_sheet <- "variable_matrix"

# ---------- 2. 检查依赖包 ----------
required_packages <- c("readxl", "openxlsx", "ggplot2", "cluster")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "缺少所需 R 包：", paste(missing_packages, collapse = ", "),
    "\n请先运行：install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "), "))"
  )
}

if (!file.exists(input_excel)) stop("找不到输入文件：", input_excel)
available_sheets <- readxl::excel_sheets(input_excel)
if (!input_sheet %in% available_sheets) stop("Excel 中缺少工作表：", input_sheet)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 3. Hamming 距离与 PAM 聚类 ----------
cat("[1/6] 正在读取分类矩阵并计算 Hamming 距离...\n")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
pattern_input <- as.data.frame(
  readxl::read_excel(input_excel, sheet = input_sheet),
  check.names = FALSE, stringsAsFactors = FALSE
)
if (nrow(pattern_input) == 0L) stop("variable_matrix 是空表。")
if (!"Sequence_ID" %in% names(pattern_input)) stop("缺少 Sequence_ID 列。")
if (ncol(pattern_input) < 2L) stop("variable_matrix 中没有可变位点列。")
if (anyDuplicated(names(pattern_input))) stop("variable_matrix 存在重复列名。")
if (anyNA(pattern_input$Sequence_ID) || any(trimws(pattern_input$Sequence_ID) == "")) {
  stop("Sequence_ID 存在 NA 或空值。")
}
if (anyDuplicated(pattern_input$Sequence_ID)) stop("Sequence_ID 存在重复值。")
pattern_ids <- as.character(pattern_input$Sequence_ID)
pattern_matrix <- pattern_input[, setdiff(names(pattern_input), "Sequence_ID"),
                                drop = FALSE]
pattern_matrix[] <- lapply(pattern_matrix, function(x) toupper(trimws(as.character(x))))
if (anyNA(pattern_matrix) || any(pattern_matrix == "")) {
  stop("分类矩阵存在 NA 或空状态。")
}
pattern_variables <- names(pattern_matrix)
if (any(!grepl("^Pos_[0-9]+$", pattern_variables))) {
  stop("位点列名必须符合 Pos_数字 格式。")
}
allowed_states <- c(strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]], "-", "X", "?", "*")
unexpected_states <- setdiff(sort(unique(unlist(pattern_matrix, use.names = FALSE))),
                             allowed_states)
if (length(unexpected_states) > 0L) {
  stop("发现非预期状态：", paste(unexpected_states, collapse = ", "))
}
pattern_n <- nrow(pattern_matrix)
if (pattern_n < 4L) stop("PAM 聚类至少需要 4 条序列。")
pattern_p <- ncol(pattern_matrix)
pattern_array <- as.matrix(pattern_matrix)

hamming_matrix <- matrix(0, nrow = pattern_n, ncol = pattern_n,
                         dimnames = list(pattern_ids, pattern_ids))
for (i in seq_len(pattern_n - 1L)) {
  for (j in (i + 1L):pattern_n) {
    distance_value <- mean(pattern_array[i, ] != pattern_array[j, ])
    hamming_matrix[i, j] <- distance_value
    hamming_matrix[j, i] <- distance_value
  }
}
hamming_distance <- stats::as.dist(hamming_matrix)
saveRDS(hamming_distance, file.path(output_dir, "hamming_distance.rds"))
utils::write.csv(hamming_matrix, file.path(output_dir, "hamming_distance_matrix.csv"),
                 row.names = TRUE)

# 调整Rand指数：比较同一批样本在原始分组和重复子抽样分组中的一致性。
adjusted_rand_index <- function(labels_a, labels_b) {
  contingency <- table(labels_a, labels_b)
  choose2 <- function(x) x * (x - 1) / 2
  sum_cells <- sum(choose2(contingency))
  sum_rows <- sum(choose2(rowSums(contingency)))
  sum_cols <- sum(choose2(colSums(contingency)))
  total_pairs <- choose2(sum(contingency))
  if (total_pairs == 0) return(NA_real_)
  expected <- sum_rows * sum_cols / total_pairs
  maximum <- 0.5 * (sum_rows + sum_cols)
  denominator <- maximum - expected
  if (abs(denominator) < .Machine$double.eps) return(1)
  (sum_cells - expected) / denominator
}

cat("[2/6] 正在比较不同模式数量...\n")
k_values <- 2L:min(6L, pattern_n - 1L)
subsample_repeats <- 200L
subsample_size <- max(2L, floor(0.80 * pattern_n))
set.seed(20260721L)
pam_models <- vector("list", length(k_values))
names(pam_models) <- as.character(k_values)
model_rows <- vector("list", length(k_values))

for (k_index in seq_along(k_values)) {
  k <- k_values[k_index]
  fitted <- cluster::pam(hamming_distance, k = k, diss = TRUE)
  pam_models[[as.character(k)]] <- fitted
  cluster_sizes <- table(fitted$clustering)
  ari_values <- replicate(subsample_repeats, {
    selected_rows <- sort(sample(seq_len(pattern_n), subsample_size,
                                 replace = FALSE))
    sub_distance <- stats::as.dist(hamming_matrix[selected_rows, selected_rows,
                                                   drop = FALSE])
    sub_fit <- cluster::pam(sub_distance, k = k, diss = TRUE)
    adjusted_rand_index(fitted$clustering[selected_rows], sub_fit$clustering)
  })
  model_rows[[k_index]] <- data.frame(
    K = k,
    Average_Silhouette = mean(fitted$silinfo$widths[, "sil_width"]),
    Median_Silhouette = stats::median(fitted$silinfo$widths[, "sil_width"]),
    Minimum_Cluster_Size = min(cluster_sizes),
    Maximum_Cluster_Size = max(cluster_sizes),
    Mean_Subsample_ARI = mean(ari_values, na.rm = TRUE),
    SD_Subsample_ARI = stats::sd(ari_values, na.rm = TRUE),
    Subsample_Repeats = subsample_repeats,
    stringsAsFactors = FALSE
  )
}
pattern_number_metrics <- do.call(rbind, model_rows)
eligible_k <- pattern_number_metrics$Minimum_Cluster_Size >= 3L
if (!any(eligible_k)) eligible_k <- rep(TRUE, nrow(pattern_number_metrics))
selected_row <- which(eligible_k)[which.max(
  pattern_number_metrics$Average_Silhouette[eligible_k]
)]
selected_k <- pattern_number_metrics$K[selected_row]
selected_model <- pam_models[[as.character(selected_k)]]
selected_clusters <- as.integer(selected_model$clustering)

pattern_number_metrics$Selected <- pattern_number_metrics$K == selected_k
utils::write.csv(pattern_number_metrics,
                 file.path(output_dir, "pattern_number_selection.csv"),
                 row.names = FALSE)

cat("[3/6] 正在计算每个变量的熵差和置换检验...\n")
entropy_value <- function(x) {
  probabilities <- as.numeric(prop.table(table(x)))
  -sum(probabilities * log2(probabilities))
}
entropy_reduction <- function(values, groups) {
  overall <- entropy_value(values)
  group_levels <- sort(unique(groups))
  within <- sum(vapply(group_levels, function(group) {
    index <- groups == group
    mean(index) * entropy_value(values[index])
  }, numeric(1)))
  c(Overall_Entropy = overall, Within_Pattern_Entropy = within,
    Entropy_Reduction = overall - within,
    Normalized_Entropy_Reduction = if (overall > 0) (overall - within) / overall else 0)
}
observed_entropy <- t(vapply(pattern_matrix, entropy_reduction,
                             numeric(4), groups = selected_clusters))
permutation_repeats <- 2000L
set.seed(20260721L)
exceedance <- integer(pattern_p)
for (permutation_index in seq_len(permutation_repeats)) {
  permuted_clusters <- sample(selected_clusters, replace = FALSE)
  permuted_reduction <- vapply(pattern_matrix, function(values) {
    entropy_reduction(values, permuted_clusters)[["Entropy_Reduction"]]
  }, numeric(1))
  exceedance <- exceedance +
    as.integer(permuted_reduction >= observed_entropy[, "Entropy_Reduction"] - 1e-12)
}
empirical_p <- (exceedance + 1) / (permutation_repeats + 1)
variable_importance <- data.frame(
  Variable = pattern_variables,
  observed_entropy,
  Permutation_P = empirical_p,
  FDR = stats::p.adjust(empirical_p, method = "BH"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
variable_importance <- variable_importance[
  order(-variable_importance$Entropy_Reduction,
        variable_importance$FDR, variable_importance$Variable),
  , drop = FALSE
]
row.names(variable_importance) <- NULL
utils::write.csv(variable_importance,
                 file.path(output_dir, "pattern_defining_variables.csv"),
                 row.names = FALSE)

cat("[4/6] 正在整理模式成员和状态概率...\n")
medoid_rows <- selected_model$id.med
medoid_ids <- pattern_ids[medoid_rows]
membership <- data.frame(
  Sequence_ID = pattern_ids,
  Pattern = paste0("Pattern_", selected_clusters),
  Distance_to_Medoid = vapply(seq_len(pattern_n), function(i) {
    hamming_matrix[i, medoid_rows[selected_clusters[i]]]
  }, numeric(1)),
  Silhouette_Width = selected_model$silinfo$widths[
    match(pattern_ids, rownames(selected_model$silinfo$widths)),
    "sil_width"
  ],
  stringsAsFactors = FALSE
)
utils::write.csv(membership, file.path(output_dir, "pattern_membership.csv"),
                 row.names = FALSE)

state_probability_rows <- list()
probability_row_index <- 1L
modal_rows <- list()
for (cluster_id in seq_len(selected_k)) {
  cluster_index <- selected_clusters == cluster_id
  modal_states <- character(pattern_p)
  modal_frequencies <- numeric(pattern_p)
  for (position_index in seq_along(pattern_variables)) {
    position <- pattern_variables[position_index]
    overall_table <- prop.table(table(pattern_matrix[[position]]))
    cluster_table <- prop.table(table(pattern_matrix[[position]][cluster_index]))
    states <- sort(unique(pattern_matrix[[position]]))
    cluster_probabilities <- setNames(rep(0, length(states)), states)
    cluster_probabilities[names(cluster_table)] <- as.numeric(cluster_table)
    overall_probabilities <- setNames(rep(0, length(states)), states)
    overall_probabilities[names(overall_table)] <- as.numeric(overall_table)
    modal_states[position_index] <- names(cluster_probabilities)[
      which.max(cluster_probabilities)
    ]
    modal_frequencies[position_index] <- max(cluster_probabilities)
    for (state in states) {
      state_probability_rows[[probability_row_index]] <- data.frame(
        Pattern = paste0("Pattern_", cluster_id),
        Variable = position,
        State = state,
        Count = sum(pattern_matrix[[position]][cluster_index] == state),
        Pattern_Frequency = cluster_probabilities[[state]],
        Overall_Frequency = overall_probabilities[[state]],
        Frequency_Difference = cluster_probabilities[[state]] -
          overall_probabilities[[state]],
        stringsAsFactors = FALSE
      )
      probability_row_index <- probability_row_index + 1L
    }
  }
  modal_rows[[cluster_id]] <- data.frame(
    Pattern = paste0("Pattern_", cluster_id),
    Sample_Count = sum(cluster_index),
    Sample_Frequency = mean(cluster_index),
    Medoid_Sequence_ID = medoid_ids[cluster_id],
    Mean_Distance_to_Medoid = mean(membership$Distance_to_Medoid[cluster_index]),
    Mean_Silhouette_Width = mean(membership$Silhouette_Width[cluster_index]),
    Mean_Modal_State_Frequency = mean(modal_frequencies),
    as.list(setNames(modal_states, pattern_variables)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
pattern_state_probabilities <- do.call(rbind, state_probability_rows)
pattern_profiles <- do.call(rbind, modal_rows)
row.names(pattern_profiles) <- NULL
utils::write.csv(pattern_state_probabilities,
                 file.path(output_dir, "pattern_state_probabilities.csv"),
                 row.names = FALSE)
utils::write.csv(pattern_profiles,
                 file.path(output_dir, "major_pattern_profiles.csv"),
                 row.names = FALSE)

cat("[5/6] 正在输出工作簿和必要图形...\n")
pattern_wb <- openxlsx::createWorkbook()
pattern_tables <- list(
  Number_Selection = pattern_number_metrics,
  Pattern_Profiles = pattern_profiles,
  Membership = membership,
  Defining_Variables = variable_importance,
  State_Probabilities = pattern_state_probabilities
)
for (sheet_name in names(pattern_tables)) {
  openxlsx::addWorksheet(pattern_wb, sheet_name)
  openxlsx::writeData(pattern_wb, sheet_name, pattern_tables[[sheet_name]])
  openxlsx::freezePane(pattern_wb, sheet_name, firstRow = TRUE)
  openxlsx::setColWidths(pattern_wb, sheet_name,
                         cols = seq_len(ncol(pattern_tables[[sheet_name]])),
                         widths = "auto")
}
openxlsx::saveWorkbook(pattern_wb,
                       file.path(output_dir, "major_pattern_analysis.xlsx"),
                       overwrite = TRUE)

selection_plot <- ggplot2::ggplot(
  pattern_number_metrics,
  ggplot2::aes(x = K, y = Average_Silhouette)
) +
  ggplot2::geom_line(color = "#3B78A8") +
  ggplot2::geom_point(ggplot2::aes(color = Selected), size = 3) +
  ggplot2::scale_x_continuous(breaks = k_values) +
  ggplot2::labs(title = "Pattern number selection", x = "Number of patterns",
                y = "Average silhouette width") +
  ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(file.path(output_dir, "pattern_number_selection.pdf"),
                selection_plot, width = 7, height = 5, device = "pdf")

top_variables <- utils::head(variable_importance, 20L)
top_variables$Variable <- factor(top_variables$Variable,
                                 levels = rev(top_variables$Variable))
variable_plot <- ggplot2::ggplot(
  top_variables,
  ggplot2::aes(x = Entropy_Reduction, y = Variable, fill = FDR < 0.05)
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(values = c("TRUE" = "#3B78A8", "FALSE" = "#A9A9A9")) +
  ggplot2::labs(title = "Top pattern-defining variables",
                x = "Entropy reduction", y = "Variable", fill = "FDR < 0.05") +
  ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(file.path(output_dir, "pattern_defining_variables.pdf"),
                variable_plot, width = 8, height = 6, device = "pdf")

saveRDS(
  list(
    selected_k = selected_k,
    selected_model = selected_model,
    pattern_number_metrics = pattern_number_metrics,
    membership = membership,
    variable_importance = variable_importance,
    state_probabilities = pattern_state_probabilities,
    profiles = pattern_profiles,
    parameters = list(
      distance = "proportional Hamming distance",
      k_values = k_values,
      selection_rule = "maximum average silhouette among solutions with minimum cluster size >= 3",
      subsample_fraction = subsample_size / pattern_n,
      subsample_repeats = subsample_repeats,
      permutation_repeats = permutation_repeats,
      random_seed = 20260721L
    )
  ),
  file.path(output_dir, "major_pattern_analysis.rds")
)

cat("[6/6] 分析完成。\n")
cat("选择的主要模式数：", selected_k, "\n", sep = "")
cat("平均轮廓系数：",
    sprintf("%.4f", pattern_number_metrics$Average_Silhouette[selected_row]),
    "\n", sep = "")
cat("平均子抽样 ARI：",
    sprintf("%.4f", pattern_number_metrics$Mean_Subsample_ARI[selected_row]),
    "\n", sep = "")
cat("FDR < 0.05 的区分变量数：", sum(variable_importance$FDR < 0.05),
    "\n", sep = "")
cat("结果目录：", normalizePath(output_dir), "\n", sep = "")
warnings()

