#' @title volcanoTriangleInator
#'
#' @description Adds triangular shapes to a volcano plot to highlight sides of a contrast in a differential expression analysis.
#' @param df Data frame containing differential expression results with a 'log2FoldChange' column (DESeq2).
#' @param xmin Minimum x-coordinate for the left triangle (negative log2 fold change).
#' @param xmax Maximum x-coordinate for the right triangle (positive log2 fold change).
#' @param ybase Base y-coordinate for the triangles.
#' @param yheight Height of the triangles.
#' @param left_label Label for the left triangle (default: "Left").
#' @param right_label Label for the right triangle (default: "Right").
#' @param triangle_fill Vector of two colors for the triangles (default: c("skyblue", "orange")).
#' @param triangle_alpha Transparency level for the triangles (default: 0.7).
#' @param label_size Size of the labels (default: 5).
#' @param left_label_hjust Horizontal justification for the left label (default: -Inf).
#' @param right_label_hjust Horizontal justification for the right label (default: Inf).
#' @param left_bump Horizontal adjustment for the left label position (default: 0).
#' @param right_bump Horizontal adjustment for the right label position (default: 0).
#' @return A list of ggplot2 layers to add to a volcano plot.
#' @import ggplot2
#'
#' @export

volcanoTriangleInator <- function(df, xmin, xmax, ybase, yheight,
                                    left_label = "Left", right_label = "Right",
                                    triangle_fill = c("skyblue", "orange"),
                                    triangle_alpha = 0.7, label_size = 5,
                                    left_label_hjust = -Inf, right_label_hjust = Inf,
                                    left_bump = 0, right_bump = 0) {
  #default imputation if a df is supplied
  if (is.null(xmin) | is.null(xmax) | is.null(ybase) | is.null(yheight)) {
    ybase <- -1
    yheight <- 1
    xmin <- -max(abs(df$log2FoldChange, na.rm = TRUE))
    xmax <- max(abs(df$log2FoldChange, na.rm = TRUE))
  }

  #triangle coordinates
  triangle_left <- data.frame(
    x = c(xmin, xmin, 0),
    y = c(ybase, ybase + yheight, ybase)
  )
  triangle_right <- data.frame(
    x = c(xmax, xmax, 0),
    y = c(ybase, ybase + yheight, ybase)
  )
  #label positions (at the midpoint of the shortest side)
  left_label_pos <- data.frame(
    x = (xmin + left_bump),
    y = ybase - 0.1 * yheight
  )
  right_label_pos <- data.frame(
    x = (xmax + right_bump) / 2,
    y = ybase - 0.1 * yheight
  )
  list(
    geom_polygon(data = triangle_left, aes(x, y), fill = triangle_fill[1], alpha = triangle_alpha, inherit.aes = FALSE),
    geom_polygon(data = triangle_right, aes(x, y), fill = triangle_fill[2], alpha = triangle_alpha, inherit.aes = FALSE),
    annotate("text", x = left_label_pos$x, y = left_label_pos$y, label = left_label, hjust = left_label_hjust, size = label_size),
    annotate("text", x = right_label_pos$x, y = right_label_pos$y, label = right_label, hjust = right_label_hjust, size = label_size)
  )
}
