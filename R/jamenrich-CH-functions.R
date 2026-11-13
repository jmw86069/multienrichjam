

#' ComplexHeatmap method to edit viewport in legend grob
#' @param gtree grid `gTree` object
#' @param ... additional arguments are ignored.
#' @keywords internal
#' @noRd
CH_edit_vp_in_legend_grob <- function(
   gtree,
   ...
) {
   # borrowed from unexported ComplexHeatmap functions
	if(is.null(gtree$vp)) {
		vp_param = list(...)
		nm = names(vp_param)
		if ("valid.just" %in% nm) {
			valid.just = vp_param$valid.just
		} else {
			valid.just = NULL
		}
		vp_param = vp_param[names(vp_param) != "valid.just"]
		gtree$vp = do.call(grid::viewport, vp_param)
		if(!is.null(valid.just)) {
			gtree$vp$valid.just = valid.just
			gtree$vp$justification = valid.just
		}
	} else {
		vp_param = list(...)
		nm = names(vp_param)
		if("just" %in% nm) {
			if(is.numeric(vp_param$just)) {
				vp_param$valid.just = vp_param$just
			} else {
				vp_param$valid.just = CH_valid_just(vp_param$just)
			}
		}
		if("justification" %in% nm) {
			if(is.numeric(vp_param$justification)) {
				vp_param$valid.just = vp_param$justification
			} else {
				vp_param$valid.just = CH_valid_just(vp_param$justification)
			}
		} else {
			vp_param$justification = vp_param$valid.just
		}
		
		for(nm in names(vp_param)) {
			gtree$vp[[nm]] = vp_param[[nm]]
		}
	}

	# gtree$vp$name = legend_vp_name()

	return(gtree)
}

#' ComplexHeatmap check for valid justification
#' @param just `numeric` or `character` justification
#' @keywords internal
#' @noRd
CH_valid_just <- function(
   just
) {
   # borrowed from unexported ComplexHeatmap functions
	if(is.numeric(just)) return(just)
	if(length(just) == 1) {
		just = switch(just,
			"centre" = c("center", "center"),
			"center" = c("center", "center"),
			"left" = c("left", "center"),
			"right" = c("right", "center"),
			"top" = c("center", "top"),
			"bottom" = c("center", "bottom"),
			c("center", "center")
      )
	}
	if(length(just) != 2) {
		cli::cli_abort(
         "{.var just} should be a single character or a vector of length 2.")
	}
	j = c(
      "center" = 0.5,
      "left" = 0,
      "right" = 1,
      "top" = 1,
      "bottom" = 0
   )
	if(is.character(just)) {
		just = j[just]
	} else if(!is.numeric(just)) {
		stop_wrap("`just` can only be character or numeric.")
	}
	return(unname(just))
}
