
#' Handle row and column split parameters for gene-pathway data
#' 
#' @returns `list` with
#'    * Mem
#'    * row_split
#'    * row_title
#'    * cluster_rows
#'    * column_split
#'    * column_title
#'    * cluster_columns
#'    
#' @keywords internal
handle_rowcol_splits <- function
(Mem,
 auto_split=TRUE,
 row_split,
 row_title,
 max_row_split=12,
 cluster_rows,
 row_method,
 gene_im_weight=0.5,
 column_split,
 column_title,
 max_column_split=8,
 cluster_columns,
 column_method,
 enrich_im_weight=0.3,
 trim_rows=TRUE,
 trim_columns=TRUE,
 p_cutoff=NULL,
 p_floor=NULL,
 seed=123,
 verbose=FALSE,
 debug=FALSE,
 ...)
{
   #
   use_genes <- genes(Mem);
   use_sets <- sets(Mem);
   
   # debug upfront
   if (debug) {
      jamba::printDebug("1row_split:");print(row_split);# debug
      jamba::printDebug("1row_title:");print(row_title);# debug
      jamba::printDebug("1column_split:");print(column_split);# debug
      jamba::printDebug("1column_title:");print(column_title);# debug
   }   
   
   # split functions
   row_split_fn <- function(n, use_values, min_ct=5, max_x_split=12)
   {
      if (n == 0 || n > length(use_values)) {
         if (length(use_values) <= min_ct) {
            n <- 1;
         } else {
            n <- jamba::noiseFloor(floor(length(use_values)^(1/2.5)),
               ceiling=max_x_split)
         }
      }
      n
   }
   column_split_fn <- function(n, use_values,
      min_ct=5, max_x_split=8, power=1/2, div=1.5)
   {
      if (n == 0 || n > length(use_values)) {
         if (length(use_values) <= min_ct) {
            n <- ceiling(sqrt(length(use_values)));
         } else {
            ncol_x <- (min_ct + length(use_values)/div) ^ power;
            n <- jamba::noiseFloor(floor(ncol_x),
               ceiling=max_x_split)
         }
      }
      n
   }

   # helper function to take a peek at clustering results to see if
   # the requested split is valid for this data
   peek_cluster_split <- function
   (use_matrix, 
    x_split, 
    x_title, 
    x_method, 
    x_cluster=NULL,
    seed=123,
    verbose=FALSE)
   {
      #
      if (length(seed) > 0) {
         set.seed(head(seed, 1));
      }
      # pre-calculate the clustering
      if (is.function(x_cluster)) {
         # custom clustering function
         cluster_x <- x_cluster(use_matrix)
         # Todo: Consider converting dendrogram to hclust for cutree(),
         # otherwise inform user to produce data usable by cutree().
      } else {
         cluster_x <- amap::hcluster(
            link="ward",
            use_matrix,
            method=x_method);
      }
      # 0.0.101.900 - use height instead of cutree k
      if (x_split > 1) {
         cc_h <- c(sort(decreasing=TRUE,
            cluster_x$height), 0);
         # if (verbose){jamba::printDebug("cc_h:");print(cc_h);}# debug
         cc_height <- mean(c(cc_h[x_split - 1], cc_h[x_split]));
         
         # if (verbose){jamba::printDebug("cc_height:");print(cc_height);}# debug
         if (cc_height == 0 && any(cc_h > 0)) {
            cc_height <- mean(c(0, min(cc_h[cc_h > 0])))
         }
         ccnum <- stats::cutree(cluster_x, h=cc_height);
         x_split <- length(unique(ccnum));
         # 0.0.101.900 - shorten or expand as needed?
         if (length(x_title) != x_split) {
            x_title <- jamba::makeNames(
               rep(x_title,
                  length.out=x_split),
               ...);
         }
      } else {
         x_split <- NULL;
         if (length(x_title) > 1) {
            x_title <- NULL;
         }
      }
      return(list(x_split=x_split,
         x_title=x_title,
         cluster_x=cluster_x))
   }
   
   
   # helper function for rows or columns
   calc_split <- function
   (x_split,
    x_title, 
    x_values, 
    max_x_split, 
    x_split_fn, 
    trim_values=TRUE,
    auto_split=TRUE,
    verbose=FALSE)
   {
      # x_split is empty
      if (length(x_split) == 0) {
         # auto-define cluster counts
         if (auto_split) {
            if (verbose) jamba::printDebug("auto_split");# verbose
            # for more than 5 rows, use scaling logic up to 12 gene clusters
            x_split <- x_split_fn(n=0,
               use_values=x_values,
               max_x_split=max_x_split)
            if (verbose) {
               jamba::printDebug("handle_rowcol_splits(): ",
                  "auto_split x_split:", x_split,
                  ", x_values:", x_values,
                  ", max_x_split:", max_x_split); # verbose
            }
         } else {
            x_split <- 1;
         }
      }
      # 0.0.105.900: convert x_split as list into a named vector
      if (inherits(x_split, "list")) {
         # expand into factor vector
         # using levels that match the order of list names
         x_split <- jamba::nameVector(
            factor(
               rep(names(x_split), lengths(x_split)),
               levels=names(x_split)),
            unlist(x_split));
      }
      
      retvals <- list();
      retvals$trimmed <- FALSE;
      #
      # method by type of split
      if (is.numeric(x_split) && length(x_split) == 1) {
         #
         # x_split is numeric
         if (verbose) jamba::printDebug("numeric x_split:", x_split);# verbose
         if (x_split > length(x_values)) {
            x_split <- x_split_fn(n=0,
               use_values=x_values,
               max_x_split=max_x_split)
            if (verbose) jamba::printDebug("x_split_fn() x_split:", x_split,
               ", x_values:", x_values,
               ", max_x_split:", max_x_split); # verbose
            
         }
         if (x_split < 1) {
            x_split <- 1;
         }
         if (x_split %in% 1 && length(x_title) > 1) {
            x_title <- NULL;
         }
      } else {
         #
         # x_split is a vector or data.frame
         cluster_x_slices <- FALSE;
         #
         if (is.atomic(x_split)) {
            if (length(names(x_split)) == 0) {
               stop("x_split supplied as vector must have names.")
            }
            # atomic vector of values
            # align names(x_split) with x_values
            if (!any(x_values %in% names(x_split))) {
               stop("names(x_split) do not match any x_values");
            }
            if (TRUE %in% trim_values) {
               if (!length(x_values) == length(x_split) ||
                  !all(x_values == names(x_split))) {
                  retvals$trimmed <- TRUE;
               }
               x_values <- intersect(x_values, names(x_split));
               x_split <- x_split[x_values];
               # memIM <- memIM[x_values, , drop=FALSE];
            } else {
               gmatch <- match(x_values, names(x_split))
               x_split <- x_split[gmatch];
               names(x_split) <- x_values;
            }
            x_title <- as.character(levels(
               factor(x_split, exclude=NULL)));
         } else if (inherits(x_split, "data.frame")) {
            # data.frame of column values
            # align rownames(x_split) with x_values
            if (!any(x_values %in% rownames(x_split))) {
               stop("rownames(x_split) do not match any x_values");
            }
            x_split_v <- jamba::pasteByRowOrdered(x_split);
            names(x_split_v) <- rownames(x_split);
            if (TRUE %in% trim_values) {
               if (!all(x_values == names(x_split))) {
                  retvals$trimmed <- TRUE;
               }
               x_values <- intersect(x_values, names(x_split_v));
               x_split <- x_split_v[x_values];
               # memIM <- memIM[x_values, , drop=FALSE];
            } else {
               gmatch <- match(x_values, names(x_split))
               x_split <- x_split_v[gmatch];
               names(x_split) <- x_values;
            }
            x_title <- as.character(levels(
               factor(x_split, exclude=NULL)));
         } else {
            x_split <- 1;
            if (length(x_title) > 1) {
               x_title <- NULL;
            }
         }
      }
      retvals$x_split <- x_split;
      retvals$x_title <- x_title;
      return(retvals);
   }

   ## Validate row_split, row_title, column_split, column_title
   ############################################
   ## row_split
   ##
   row_split_list <- calc_split(x_split=row_split,
      x_title=row_title,
      x_values=use_genes,
      max_x_split=max_row_split,
      x_split_fn=row_split_fn,
      trim_values=trim_rows,
      auto_split=auto_split,
      verbose=FALSE)
   row_split <- row_split_list$x_split;
   row_title <- row_split_list$x_title;
   
   # Verify: subset rows
   if (TRUE %in% row_split_list$trimmed) {
      Mem <- Mem[names(row_split), , ];
      use_genes <- genes(Mem);
   }
   
   ############################################
   ## column_split
   ##
   column_split_list <- calc_split(x_split=column_split,
      x_title=column_title,
      x_values=use_sets,
      max_x_split=max_column_split,
      x_split_fn=column_split_fn,
      trim_values=trim_columns,
      auto_split=auto_split,
      verbose=isTRUE(debug))
   column_split <- column_split_list$x_split;
   column_title <- column_split_list$x_title;
   
   # Verify: subset columns
   if (TRUE %in% column_split_list$trimmed) {
      Mem <- Mem[, names(column_split), ];
      use_sets <- sets(Mem)
   }
   
   if (debug) {
      jamba::printDebug("2row_split:");print(row_split);# debug
      jamba::printDebug("2row_title:");print(row_title);# debug
      jamba::printDebug("2column_split:");print(column_split);# debug
      jamba::printDebug("2column_title:");print(column_title);# debug
   }   

   ############################################
   ## Cluster rows
   if (length(cluster_rows) == 0 ||
         isTRUE(cluster_rows) ||
         is.function(cluster_rows)) {
      #
      # generate weighted geneIM,memIM matrix
      if (!length(gene_im_weight) == 1 || any(gene_im_weight > 1)) {
         gene_im_weight <- 0.5;
      }
      gene_weight <- round(gene_im_weight * 100) / 10;
      im_weight <- 10 - gene_weight;
      min_weight <- max(c(1, min(c(gene_weight, im_weight))));
      gene_weight <- gene_weight / min_weight;
      im_weight <- im_weight / min_weight;
      use_gene_im <- geneIM(Mem) * jamba::rmNA(naValue=1,
         geneIMdirection(Mem))
      
      # generate the combined matrix to use during clustering
      if (im_weight == 0) {
         row_matrix <- (use_gene_im[use_genes, , drop=FALSE]) * gene_weight;
      } else if (gene_weight == 0) {
         row_matrix <- ((memIM(Mem)[use_genes, use_sets, drop=FALSE]) *
               im_weight); # not incidence
      } else {
         row_matrix <- cbind(
            (use_gene_im[use_genes, , drop=FALSE]) * gene_weight,
            ((memIM(Mem)[use_genes, use_sets, drop=FALSE])) * im_weight)
         # above: does not convert memIM to incidence, keeps counts
         #
         # below converts to incidence matrix with 1 or 0
         # ((memIM(Mem)[use_genes, use_sets, drop=FALSE] != 0) * 1) * im_weight)
      }

      # for numeric split number, peek at data to confirm the number is valid

      if (is.numeric(row_split) && length(row_split) == 1) {
         row_peek <- peek_cluster_split(use_matrix=row_matrix,
            x_split=row_split,
            x_title=row_title,
            x_method=row_method,
            x_cluster=cluster_rows,
            verbose=FALSE,
            seed=seed)
         row_split <- row_peek$x_split;
         row_title <- row_peek$x_title;
         cluster_rows <- row_peek$cluster_x;
      } else {
         cluster_row_slices <- FALSE;
         use_cluster_rows <- cluster_rows;
         cluster_rows <- function(x, ...){
            if (length(seed) > 0) {
               set.seed(head(seed, 1));
            }
            userows <- rownames(x);
            usematch <- match(userows, rownames(row_matrix));
            usematrix <- jamba::rmNA(naValue=0,
               row_matrix[usematch, , drop=FALSE]);
            if (is.function(use_cluster_rows)) {
               use_cluster_rows(usematrix)
            } else {
               amap::hcluster(
                  usematrix,
                  link="ward",
                  method=row_method);
            }
         }
      }
   }
   
   ############################################
   ## Cluster columns
   if (length(cluster_columns) == 0 ||
         isTRUE(cluster_columns) ||
         is.function(cluster_columns)) {
      # Assemble the P-value matrix with gene incidence matrix
      # and cluster altogether, which has the benefit/goal of
      # accentuating similar enrichment profiles which also have
      # similar gene content.
      if (!length(enrich_im_weight) == 1 || any(enrich_im_weight > 1)) {
         enrich_im_weight <- 0.3;
      }
      enrich_weight <- round(enrich_im_weight * 10);
      im_weight <- 10 - enrich_weight;
      min_weight <- max(c(1, min(c(enrich_weight, im_weight))));
      enrich_weight <- enrich_weight / min_weight;
      im_weight <- im_weight / min_weight;
      ## 0.0.31.900 use column_matrix with enrich_im_weight adjustment
      column_matrix <- cbind(
         jamba::noiseFloor(
            -log10(enrichIM(Mem)[use_sets, , drop=FALSE]),
            minimum=-log10(p_cutoff + 1e-5),
            newValue=0,
            ceiling=-log10(p_floor)) * enrich_weight,
         t((memIM(Mem)[use_genes, use_sets, drop=FALSE]) * im_weight) # non-incidence
         # t((memIM(Mem)[use_genes, use_sets, drop=FALSE] != 0) * 1) * im_weight # im
      );
      
      if (is.numeric(column_split) && length(column_split) == 1) {
         column_peek <- peek_cluster_split(use_matrix=column_matrix,
            x_split=column_split,
            x_title=column_title,
            x_method=column_method,
            x_cluster=cluster_columns,
            verbose=FALSE,
            seed=seed)
         column_split <- column_peek$x_split;
         column_title <- column_peek$x_title;
         cluster_columns <- column_peek$cluster_x;
      } else {
         cluster_column_slices <- FALSE;
         use_cluster_columns <- cluster_columns;
         cluster_columns <- function(x, ...){
            if (length(seed) > 0) {
               set.seed(head(seed, 1));
            }
            userows <- rownames(x);
            usematch <- match(userows, rownames(column_matrix));
            use_x <- jamba::rmNA(naValue=0,
               column_matrix[usematch, , drop=FALSE]);
            if (is.function(use_cluster_columns)) {
               use_cluster_columns(use_x)
            } else {
               amap::hcluster(use_x,
                  link="ward",
                  method=column_method)
            }
         }
      }
   }
   if (debug) {
      jamba::printDebug("3row_split:");print(row_split);# debug
      jamba::printDebug("3row_title:");print(row_title);# debug
      jamba::printDebug("3column_split:");print(column_split);# debug
      jamba::printDebug("3column_title:");print(column_title);# debug
   }   
   
   if (debug) {
      withr::with_par(list(mfrow=c(2, 1)), {
         if (!is.function(cluster_columns)) {
            plot(cluster_columns, main="cluster_columns")
         } else {
            jamba::nullPlot(plotAreaTitle="cluster_columns is a function")
         }
         if (!is.function(cluster_rows)) {
            plot(cluster_rows, main="cluster_rows")
         } else {
            jamba::nullPlot(plotAreaTitle="cluster_rows is a function")
         }
      })
   }
   
   retvals <- list();
   retvals$Mem <- Mem;
   retvals$row_split <- row_split;
   retvals$row_title <- row_title;
   retvals$cluster_rows <- cluster_rows;
   retvals$column_split <- column_split;
   retvals$column_title <- column_title;
   retvals$cluster_columns <- cluster_columns;
   retvals;
}

