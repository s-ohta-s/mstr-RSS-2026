# plot_frobenius_boxplots_v2.R
# dgp = 1111 / 1011 / 0111 について、行列 R / Lambda / A / Sigma とベクトル beta の
# フロベニウスノルム誤差を rep ごとに計算し、
# 横軸 n × paramset (weak, moderate, strong) のボックスプロットを描く。
# 誤差は 2 系統:
#   pAIC : estimate 列      (pAIC 選択 γ でのフィット)
#   pBIC : estimate_pBIC 列 (pBIC 選択 γ でのフィット; 列が無ければスキップ)
# 縦軸は 3 種:
#   abs      : ||theta_hat - theta_true||_F
#   rel_true : ||theta_hat - theta_true||_F / ||theta_true||_F
#   rel_est  : ||theta_hat - theta_true||_F / ||theta_hat||_F
# Sigma は対称行列だが、CSV に格納された上三角要素 ([1,1],[1,2],[2,2]) を
# そのまま用いる（非対角要素は 1 回だけカウントする定義）。
# 入力:  <base>/n<N>/results_A.csv
# 出力:  <figdir>/frob_<crit>_<metric>_<group>_dgp<dgp>.png      (線形軸)
#        <figdir>/frob_<crit>_<metric>_<group>_dgp<dgp>_log.png  (log10 軸)
#        <figdir>/frob_<crit>_rel_true_<group>_dgp<dgp>_log_shared.png
#                 (log10 軸 + 同一グループの dgp 間で縦軸を統一)
#        <figdir>/gamma_<crit>_dgp<dgp>.png                      (γ, log10 軸)
#   (dgp 1011 に Lambda は、dgp 0111 に R は存在しないためスキップ)
# 使い方:
#   Rscript scripts/plot_frobenius_boxplots_v2.R [--base=output/sim_full_v2]
#                                                [--ns=100,400,900]
#                                                [--figdir=<base>/figures_v2]
#                                                [--crits=pAIC,pBIC]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

.args <- commandArgs(trailingOnly = TRUE)
optval <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

base      <- optval("base", "output/sim_full_v2")
ns        <- as.integer(strsplit(optval("ns", "100,400,900"), ",")[[1]])
outdir    <- optval("figdir", file.path(base, "figures_v2"))
crits_req <- strsplit(optval("crits", "pAIC,pBIC"), ",")[[1]]

dirs <- setNames(file.path(base, paste0("n", ns)), ns)
dgps      <- c("1111", "1011", "0111")
# CSV 上の値 -> 図中の表示ラベル
paramsets <- c("weak", "mid", "strong")
param_labels <- c(weak = "weak", mid = "moderate", strong = "strong")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

res <- bind_rows(lapply(dirs, function(d) {
  read.csv(file.path(d, "results_A.csv"),
           colClasses = c(dgp = "character"), stringsAsFactors = FALSE)
})) %>%
  filter(dgp %in% dgps, paramset %in% paramsets)

has_pbic <- "estimate_pBIC" %in% names(res)
if (!has_pbic) message("estimate_pBIC 列が無いため pAIC 系統のみ描画します")

# 1. フロベニウスノルム誤差のボックスプロット (pAIC / pBIC)

grp <- res %>%
  mutate(group = case_when(
    grepl("^R\\[",      parameter) ~ "R",
    grepl("^Lambda\\[", parameter) ~ "Lambda",
    grepl("^A\\[",      parameter) ~ "A",
    grepl("^Sigma\\[",  parameter) ~ "Sigma",
    grepl("^beta_",     parameter) ~ "beta",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group))

# estimator を縦持ちにして pAIC / pBIC を同じ式で処理する
long <- grp %>%
  select(dgp, group, paramset, n, rep, parameter, true, estimate,
         any_of("estimate_pBIC")) %>%
  pivot_longer(cols = c(estimate, any_of("estimate_pBIC")),
               names_to = "crit", values_to = "est") %>%
  mutate(crit = ifelse(crit == "estimate", "pAIC", "pBIC"))

frob <- long %>%
  group_by(crit, dgp, group, paramset, n, rep) %>%
  summarise(abs      = sqrt(sum((est - true)^2)),
            rel_true = abs / sqrt(sum(true^2)),
            rel_est  = abs / sqrt(sum(est^2)),
            .groups = "drop")

n_na <- sum(!complete.cases(frob[, c("abs", "rel_true", "rel_est")]))
if (n_na > 0) {
  message(sprintf("NA の推定値を含む rep を %d 件除外しました", n_na))
  frob <- frob[complete.cases(frob[, c("abs", "rel_true", "rel_est")]), ]
}

frob <- frob %>%
  mutate(n = factor(n, levels = sort(ns)),
         paramset = factor(paramset, levels = paramsets,
                           labels = param_labels[paramsets]))

# dgp 間で縦軸を揃えるための rel_true の値域 (crit × group ごと、全 dgp 共通)
rng_rel_true <- frob %>%
  group_by(crit, group) %>%
  summarise(lo = min(rel_true), hi = max(rel_true), .groups = "drop")

hats <- c(R = "hat(R)", Lambda = "hat(Lambda)", A = "hat(A)",
          Sigma = "hat(Sigma)", beta = "hat(beta)")

ylab_expr <- function(metric, g) {
  h <- hats[[g]]
  b <- sub("hat\\((.*)\\)", "\\1", h)
  err <- sprintf("group('||', %s - %s, '||')[F]", h, b)
  txt <- switch(metric,
    abs      = err,
    rel_true = sprintf("%s / group('||', %s, '||')[F]", err, b),
    rel_est  = sprintf("%s / group('||', %s, '||')[F]", err, h))
  parse(text = txt)[[1]]
}

metric_label <- c(abs      = "Frobenius norm error",
                  rel_true = "relative Frobenius error (vs true)",
                  rel_est  = "relative Frobenius error (vs estimate)")

theme_box <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank())

crits <- intersect(crits_req, c("pAIC", if (has_pbic) "pBIC"))
n_plot <- 0
for (cr in crits) {
  for (d in dgps) {
    for (g in names(hats)) {
      dat <- filter(frob, crit == cr, dgp == d, group == g)
      if (nrow(dat) == 0) next
      for (m in c("abs", "rel_true", "rel_est")) {
        p <- ggplot(dat, aes(x = n, y = .data[[m]], fill = paramset)) +
          geom_boxplot(position = position_dodge(width = 0.8), width = 0.65,
                       outlier.size = 0.8, outlier.alpha = 0.5) +
          scale_fill_brewer(palette = "Blues", direction = 1, name = "paramset") +
          labs(x = "n", y = ylab_expr(m, g),
               title = parse(text = sprintf(
                 "paste('%s: ', %s, '  (%s, dgp %s)')",
                 metric_label[[m]], hats[[g]], cr, d))[[1]]) +
          theme_box

        f <- file.path(outdir, sprintf("frob_%s_%s_%s_dgp%s.png", cr, m, g, d))
        ggsave(f, p, width = 6.5, height = 4.5, dpi = 300)
        n_plot <- n_plot + 1

        # log10 軸版 (外れ値による箱の潰れ対策)
        f_log <- file.path(outdir,
                           sprintf("frob_%s_%s_%s_dgp%s_log.png", cr, m, g, d))
        ggsave(f_log, p + scale_y_log10(), width = 6.5, height = 4.5, dpi = 300)
        n_plot <- n_plot + 1

        # rel_true のみ: 同一グループの dgp 間で縦軸を統一した log10 軸版
        if (m == "rel_true") {
          rng <- filter(rng_rel_true, crit == cr, group == g)
          f_sh <- file.path(outdir,
                            sprintf("frob_%s_%s_%s_dgp%s_log_shared.png",
                                    cr, m, g, d))
          ggsave(f_sh, p + scale_y_log10(limits = c(rng$lo, rng$hi)),
                 width = 6.5, height = 4.5, dpi = 300)
          n_plot <- n_plot + 1
        }
      }
    }
  }
}
message(sprintf("誤差ボックスプロット: %d 枚", n_plot))

# 2. 選択された γ のボックスプロット (pAIC / pBIC)
#    γ は探索グリッド上の離散値しか取らないが、
#    ボックスプロット (縦軸 log10) で示す

gam <- res %>%
  select(dgp, paramset, n, rep, gamma, any_of("gamma_pBIC")) %>%
  distinct() %>%
  mutate(n = factor(n, levels = sort(ns)),
         paramset = factor(paramset, levels = paramsets,
                           labels = param_labels[paramsets]))

gamma_cols <- c(pAIC = "gamma")
if ("gamma_pBIC" %in% names(gam)) gamma_cols <- c(gamma_cols, pBIC = "gamma_pBIC")
gamma_cols <- gamma_cols[names(gamma_cols) %in% crits]

for (d in dgps) {
  for (cr in names(gamma_cols)) {
    dat <- filter(gam, dgp == d)
    p <- ggplot(dat, aes(x = n, y = .data[[gamma_cols[[cr]]]], fill = paramset)) +
      geom_boxplot(position = position_dodge(width = 0.8), width = 0.65,
                   outlier.size = 0.8, outlier.alpha = 0.5) +
      scale_fill_brewer(palette = "Blues", direction = 1, name = "paramset") +
      scale_y_log10() +
      labs(x = "n", y = bquote("selected " * gamma * " (" * .(cr) * ")"),
           title = bquote("selected " * gamma * " by " * .(cr) *
                            "  (dgp " * .(d) * ")")) +
      theme_box

    f <- file.path(outdir, sprintf("gamma_%s_dgp%s.png", cr, d))
    ggsave(f, p, width = 6.5, height = 4.5, dpi = 300)
    n_plot <- n_plot + 1
  }
}

message(sprintf("合計 %d 枚のプロットを出力しました -> %s", n_plot, outdir))
