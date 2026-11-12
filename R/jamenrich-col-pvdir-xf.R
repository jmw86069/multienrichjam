
#' Bivariate color interpolation function with optional floors
#' 
#' Bivariate color interpolation function with optional floors
#' 
#' @param p_cutoff `numeric` P-value threshold, default 0.01
#' @param p_floor `numeric` P-value lowest recognized value, default 1e-10
#' @param direction_cutoff `numeric` directional score cutoff, default 0
#'    applies no cutoff.
#' @param direction_max `numeric` directional maximum recognized,
#'    default 2.
#' @param color_non_hits `logical` whether to apply color when the P-value
#'    does not meet the `p_cutoff`, default FALSE uses white for non-hits.
#' @param pretty.n `integer` number of labels used with `pretty()`
#'    for rows and columns, respectively. Default is `10` and `5`.
#' @param ... additional arguments are ignored.
#' 
#' @returns `function`
#' 
#' @family jam utility functions
#' 
#' @examples
#' col_fun <- col_pvdir_xf()
#' attr(col_fun, "row_labels")
#' display_colorRamp2D(col_fun, pretty.n=c(10, 5))
#' 
#' col_fun1 <- col_pvdir_xf(direction_cutoff=1)
#' attr(col_fun1, "row_labels")
#' display_colorRamp2D(col_fun1, pretty.n=c(10, 5))
#' 
#' col_fun <- col_pvdir_xf(
#'    p_cutoff=0.05,
#'    direction_cutoff=1
#' )
#' display_colorRamp2D(col_fun)
#' 
#' @export
col_pvdir_xf <- function(
   p_cutoff=0.01,
   p_floor=1e-10,
   direction_cutoff=0,
   direction_max=2,
   color_non_hits=FALSE,
   pretty.n=c(10, 5),
   ...
) {
   #

   ## Previous default behavior
   # dir_colors <- c("royalblue4", "gold3", "firebrick3");
   # dir_colors2 <- c("skyblue", "gold", "indianred1");
   # dir_colors3 <- c("white", "white", "white");
   # dir_colors2.5 <- sapply(dir_colors2, function(i){
   #    colorjam::blend_colors(c(i, "white"))
   # })
   # mcolor2 <- jamba::rbindList(list(
   #    dir_colors3,
   #    dir_colors2.5,
   #    dir_colors2,
   #    dir_colors))
   # jamba::imageByColors(mcolor2)

   n <- 3;
   trimRamp <- c(3, 2);
   mcolor1 <- cbind(
      c("#FFFFFF", jamba::getColorRamp("Blues", n=n, trimRamp=trimRamp)),
      c("#FFFFFF", jamba::getColorRamp("Greens", n=n, trimRamp=trimRamp+c(0, 1))),
      c("#FFFFFF", jamba::getColorRamp("gold", n=n, trimRamp=trimRamp)),
      c("#FFFFFF", jamba::getColorRamp("Oranges", n=n, trimRamp=trimRamp+c(0, 1))),
      c("#FFFFFF", jamba::getColorRamp("Reds", n=n, trimRamp=trimRamp))
   )
   # jamba::imageByColors(mcolor1)
   pretty.n <- rep_len(pretty.n, 2);

   # column breaks
   column_breaks <- seq(
      from=-direction_max,
      to=direction_max,
      length.out=5
   );

   # row breaks
   p_cutoff <- 1;
   p_cutoff <- 0.05
   p_cut_lvl <- 10^(-1 * (ceiling(-log10(p_cutoff)) + 5));
   row_breaks <- c(
      -log10(p_cutoff),
      seq(from=-log10(p_cutoff),
         to=-log10(p_floor),
         length.out=3
      )
   );
   row_at <- pretty(row_breaks, n=pretty.n[1])
   row_at <- row_at[row_at >= -log10(p_cutoff) & row_at <= -log10(p_floor)]
   if (!-log10(p_cutoff) %in% row_at) {
      row_at <- c(
         # -log10(p_cutoff) - p_cut_lvl,
         -log10(p_cutoff),
         row_at
      )
   }
   row_labels <- signif(digits=3, row_breaks);

   mcolor0 <- mcolor[1, , drop=FALSE];
   mcolor0[1, ] <- "#FFFFFF";
   if (p_cutoff == 1) {
      row_at <- pretty(c(0, -log10(p_floor)), n=pretty.n[1])
      row_at <- c(
         row_at[1],
         row_at[1] + p_cut_lvl*2,
         row_at[-1]
      )
      row_label <- signif(digits=3, row_at);
      row_label[2] <- paste0(">", row_at[1])
      row_breaks[2] <- row_breaks[1] + p_cut_lvl;
      row_labels[2] <- paste0(">", row_labels[1])
   } else {
      if (isTRUE(color_non_hits)) {
         #
         row_at <- c(
            0,
            row_at
         )
         row_label <- signif(digits=3, row_at);
         row_breaks <- c(
            0,
            row_breaks#[-1]
         )
         row_labels <- c(
            0,
            row_labels#[-1]
         )
         mcolor1 <- rbind(mcolor0, mcolor1);
      } else {
         row_at <- c(
            0,
            -log10(p_cutoff) - p_cut_lvl * 2,
            row_at
         )
         row_label <- signif(digits=3, row_at);
         row_label <- c(
            paste0("<", signif(digits=3, row_at[3])),
            row_label[-1:-2]
         )
         row_breaks <- c(
            0,
            row_breaks[2] - p_cut_lvl,
            row_breaks
         )
         row_labels <- c(
            0,
            paste0("<", signif(digits=3, row_labels[2])),
            row_labels[-1]
         )
         mcolor1 <- rbind(
            mcolor0,
            rbind(
               mcolor1[1, , drop=FALSE],
               mcolor1)
         )
      }
      colnames(mcolor1) <- column_breaks;
      rownames(mcolor1) <- row_breaks;
   }
   mcolor <- mcolor1;
   
   direction_label <- NULL;
   if (direction_cutoff > direction_max) {
      direction_cutoff <- direction_max;
   }
   direction_at <- seq(0, direction_max);
   column_labels <- NULL;
   if (!direction_max %in% direction_at) {
      if (direction_max - max(direction_at) >= 0.5) {
         direction_at <- c(direction_at, direction_max)
      } else {
         direction_at[length(direction_at)] <- direction_max;
      }
   }
   if (direction_cutoff > 0) {
      direction_at <- direction_at[direction_at >= direction_cutoff];
      # Legend labels and 'at' positions
      if (!direction_cutoff %in% direction_at) {
         direction_label <- c(
            sort(-1 * direction_at),
            signif(digits=3, -1 * direction_cutoff),
            paste0(">", signif(digits=3, -1 * direction_cutoff)),
            0,
            paste0("<", signif(digits=3, direction_cutoff)),
            signif(digits=3, direction_cutoff),
            direction_at
         )
         direction_at <- c(
            sort(-1 * direction_at),
            -1 * direction_cutoff,
            -1 * direction_cutoff + 1e-5,
            0,
            direction_cutoff - 1e-5,
            direction_cutoff,
            direction_at
         )
      } else {
         direction_label <- c(
            sort(-1 * direction_at),
            paste0(">", signif(digits=3, -1 * direction_cutoff)),
            0,
            paste0("<", signif(digits=3, direction_cutoff)),
            direction_at
         )
         direction_at <- c(
            sort(-1 * direction_at),
            -1 * direction_cutoff + 1e-5,
            0,
            direction_cutoff - 1e-5,
            direction_at
         );
      }
      # Color breaks
      if (direction_cutoff >= column_breaks[4]) {
         #
         mcolor1a <- mcolor1[, c(1, 2, 3, 3, 3, 4, 5), drop=FALSE];
         column_breaks1a <- c(
            column_breaks[1],
            -direction_cutoff,
            -direction_cutoff + p_cut_lvl,
            column_breaks[3],
            direction_cutoff - p_cut_lvl,
            direction_cutoff,
            column_breaks[5]
         )
         colnames(mcolor1a) <- column_breaks1a;
         column_breaks <- column_breaks1a;
         column_labels <- c(
            column_breaks[1],
            -direction_cutoff,
            paste0(">", -direction_cutoff),
            column_breaks[3],
            paste0("<", direction_cutoff),
            direction_cutoff,
            column_breaks[5]
         )
         mcolor <- mcolor1a;
         # jamba::imageByColors(mcolor1a)
      } else {
         #
         mcolor1d <- mcolor1[, c(1, 2, 2, 3, 3, 3, 4, 4, 5), drop=FALSE];
         mcolor1add <- t(apply(mcolor1, 1, function(ic){
            circlize::colorRamp2(
               breaks=column_breaks,
               colors=ic
            )(c(-direction_cutoff, direction_cutoff))
         }))
         mcolor1d[, 3] <- mcolor1add[, 1];
         mcolor1d[, 7] <- mcolor1add[, 2];
         column_breaks1d <- c(
            column_breaks[c(1, 2)],
            -direction_cutoff,
            -direction_cutoff + p_cut_lvl,
            column_breaks[3],
            direction_cutoff - p_cut_lvl,
            direction_cutoff,
            column_breaks[c(4, 5)]
         )
         colnames(mcolor1d) <- column_breaks1d;
         column_breaks <- column_breaks1d;
         column_labels <- c(
            column_breaks[c(1, 2)],
            -direction_cutoff,
            paste0(">", -direction_cutoff),
            column_breaks[3],
            paste0("<", direction_cutoff),
            direction_cutoff,
            column_breaks[c(4, 5)]
         )
         mcolor <- mcolor1d;
         # jamba::imageByColors(mcolor1d)
      }
   } else {
      # no direction_cutoff
      column_labels <- column_breaks;
      direction_label <- c(
         paste0("<=", signif(digits=3, direction_at[1])),
         direction_at[c(-1, -length(direction_at))],
         paste0(">=", signif(digits=3, direction_at[length(direction_at)]))
      )
      direction_at <- unique(sort(c(
         -1 * direction_at,
         direction_at
      )))
   }
   col_bivariate_fun <- colorRamp2D(
      column_breaks=column_breaks,
      row_breaks=row_breaks,
      mcolor=mcolor,
      ...);
   # Add labels to attributes
   attr(col_bivariate_fun, "row_labels") <- row_labels;
   attr(col_bivariate_fun, "column_labels") <- column_labels;
   attr(col_bivariate_fun, "row_at") <- row_at;
   attr(col_bivariate_fun, "row_label") <- row_label;
   attr(col_bivariate_fun, "column_at") <- direction_at;
   attr(col_bivariate_fun, "column_label") <- direction_label;
   return(col_bivariate_fun)
}
