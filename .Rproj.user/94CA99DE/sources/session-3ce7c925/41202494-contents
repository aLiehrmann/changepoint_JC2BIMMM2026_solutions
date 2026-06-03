plot_strand_summary <- function(strand, start_w, end_w, lfc_list, cov_list, seg_list) {
  strand_symbol <- ifelse(strand == "forward", "+", "-")
  pos      <- start_w:end_w
  lfc      <- lfc_list[[strand]][pos]
  cov_mat  <- cov_list[[strand]][pos, ]
  wt_cols  <- grep("^wt",  colnames(cov_mat), value = TRUE)
  mut_cols <- grep("^mut", colnames(cov_mat), value = TRUE)

  mean_wt  <- rowMeans(log2(cov_mat[, wt_cols]  + pseudo))
  mean_mut <- rowMeans(log2(cov_mat[, mut_cols] + pseudo))
  smt      <- getSMT(seg_list[[strand]])[pos]

  strand_idx <- as.character(strand(rowRanges(dds))) == strand_symbol
  gr_all     <- rowRanges(dds)[strand_idx]
  res_all    <- as.data.frame(mcols(dds))[strand_idx, ]

  window_gr  <- GRanges("Pt", IRanges(start_w, end_w))
  ol         <- overlapsAny(gr_all, window_gr)
  gr_s       <- gr_all[ol]
  res_s      <- res_all[ol, ]
  start(gr_s) <- pmax(start(gr_s), start_w)
  end(gr_s)   <- pmin(end(gr_s),   end_w)

  genes_s <- Pt_ath[as.character(strand(Pt_ath)) == strand_symbol]
  genes_w <- genes_s[overlapsAny(genes_s, window_gr)]

  par(mfrow = c(4, 1), mar = c(1, 5, 1, 1), oma = c(3, 0, 3, 0))

  # --- panel 1 : log2-coverage ---
  ylim_cov <- range(c(mean_wt, mean_mut))
  plot(pos, mean_wt, type = "l", col = "steelblue", lwd = 0.7,
       xlim = c(start_w, end_w), ylim = ylim_cov,
       xlab = "", ylab = "log2-cov", xaxt = "n")
  lines(pos, mean_mut, col = "tomato", lwd = 0.7)
  legend("topright", legend = c("wt", "mut"),
         col = c("steelblue", "tomato"), lwd = 2, bty = "n", cex = 0.8)

  # --- panel 2 : log2-FC + piecewise-constant fit ---
  plot(pos, lfc, type = "l", col = "grey70", lwd = 0.5,
       xlim = c(start_w, end_w),
       xlab = "", ylab = expression(log[2]-FC), xaxt = "n")
  lines(pos, smt, col = "black", lwd = 1.2)
  abline(h = 0, lty = 2, col = "grey40")

  # --- panel 3 : -log10(padj) per segment ---
  neg_log10_padj <- -log10(pmax(res_s$padj, 1e-30))
  der_col        <- ifelse(res_s$DER, "firebrick", "grey70")
  plot(0, 0, type = "n",
       xlim = c(start_w, end_w),
       ylim = c(0, max(c(neg_log10_padj, -log10(0.05)), na.rm = TRUE) + 1),
       xlab = "", ylab = "-log10(padj)", xaxt = "n")
  abline(h = -log10(0.05), lty = 2, col = "grey50")
  if (length(gr_s) > 0)
    rect(start(gr_s), 0, end(gr_s), neg_log10_padj, col = der_col, border = NA)
  legend("topright", legend = c("DER (padj ≤ 0.01)", "non significant"),
         fill = c("firebrick", "grey70"), border = NA, bty = "n", cex = 0.8)

  # --- panel 4 : gene annotations ---
  plot(0, 0, type = "n",
       xlim = c(start_w, end_w), ylim = c(0, 2),
       xlab = "Genomic position (Pt)", ylab = "", yaxt = "n")
  if (length(genes_w) > 0)
    rect(start(genes_w), 0.3, end(genes_w), 1.7,
         col = "orange", border = "white", lwd = 0.3)
  legend("topright", legend = "annotated gene",
         fill = "orange", border = NA, bty = "n", cex = 0.8)
  axis(1)

  mtext(sprintf("Strand %s  [%d–%d]", strand_symbol, start_w, end_w),
        outer = TRUE, cex = 1.1, font = 2)
}
