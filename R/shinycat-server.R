
#' R-shiny server function for shinycat
#'
#' R-shiny server function for shinycat, internal use by launch_shinycat()
#'
#' @family jam shiny functions
#'
#' @param input,output,session arguments used by `shiny::shinyApp()`
#'
#' @export
shinycat_server <- function
(input,
 output,
 session)
{
   #
   # reactive variable with active selected node
   active_node <- shiny::reactiveVal(NULL);
   # reactive variable with active selected nodeset
   active_nodeset <- shiny::reactiveVal(NULL);

   ###############################################
   # adjust_node: default node adjustments
   adjust_node <- shiny::reactiveVal(data.frame(
      node=igraph::V(g)$name,
      x=0,
      y=0,
      label.dist=0,
      label.degree=0))

   ###############################################
   # adjust_nodeset: default nodeset adjustments
   va <- igraph::vertex_attr_names(g);
   nodesets <- NULL;
   if ("nodeType" %in% va) {
      nodesets <- get_cnet_nodeset(g);
      nsdf <- data.frame(node=unlist(nodesets),
         nodeset=rep(names(nodesets), lengths(nodesets)));
      # get nodeset spacing
      default_spacings <- summarize_node_spacing(g,
         scaled=TRUE,
         each_group=TRUE,
         node_groups=nodesets,
         dist_type="nearest_node")
      # default values
      default_spacing <- jamba::nameVector(
         round(
            jamba::rmInfinite(default_spacings$nearest_within[, 'Median'],
               infiniteValue=0) * 10) / 10,
         rownames(default_spacings$nearest_within))

      adjust_nodeset <- shiny::reactiveVal(data.frame(
         nodeset=names(nodesets),
         x=0,
         y=0,
         percent_spacing=default_spacing[names(nodesets)],
         rotate_degrees=0))
   } else {
      adjust_nodeset <- shiny::reactiveVal(NULL)
   }


   #########################################
   # adjust igraph data dynamically
   adjusted_cnet <- shiny::reactive({
      node_adj <- adjust_node();
      nodeset_adj <- adjust_nodeset();

      # data.table::fwrite(node_adj, file="node_adj1.txt", sep="\t");
      # data.table::fwrite(nodeset_adj, file="nodeset_adj1.txt", sep="\t");

      # nodeset adjustments
      adj_cnet <- bulk_cnet_adjustments(g,
         nodeset_adj=nodeset_adj,
         node_adj=node_adj,
         apply_negative=TRUE,
         nodesets=nodesets)

      if (FALSE) {
         # operate on a copy of the original
         adj_cnet <- g;
         if (inherits(nodeset_adj, "data.frame") && nrow(nodeset_adj) > 0) {
            # test the difference from starting point
            percent_spacing_diff <- (
               nodeset_adj$percent_spacing -
                  default_spacing[rownames(nodeset_adj)]);
            nodeset_which <- which(nodeset_adj$x != 0 |
                  nodeset_adj$y != 0 |
                  percent_spacing_diff != 0 |
                  nodeset_adj$rotate_degrees != 0);
            # adjust each
            for (i in nodeset_which) {
               adj_cnet <- adjust_cnet_nodeset(adj_cnet,
                  set_nodes=strsplit(nodeset_adj$nodeset[i], ",")[[1]],
                  x=nodeset_adj$x[i],
                  y=nodeset_adj$y[i],
                  percent_spacing=nodeset_adj$percent_spacing[i],
                  rotate_degrees=nodeset_adj$rotate_degrees[i])
            }
         }

         # node adjustments
         if (inherits(node_adj, "data.frame") && nrow(node_adj) > 0) {
            # node x,y coordinates
            node_which <- which(node_adj$x != 0 |
                  node_adj$y != 0);
            if (length(node_which) > 0) {
               adj_cnet <- nudge_igraph_node(adj_cnet,
                  nodes=node_adj$node[node_which],
                  x=node_adj$x[node_which],
                  y=node_adj$y[node_which])
            }

            # label.dist
            node_which_l <- which(node_adj$label.degree != 0 |
                  node_adj$label.dist != 0);
            if (length(node_which_l) > 0) {
               matchv <- match(node_adj$node[node_which_l],
                  igraph::V(adj_cnet)$name);
               # distance
               igraph::V(adj_cnet)[matchv]$label.dist <- (
                  igraph::V(adj_cnet)[matchv]$label.dist +
                     node_adj$label.dist[node_which_l])
               # angle
               igraph::V(adj_cnet)[matchv]$label.degree <- (
                  igraph::V(adj_cnet)[matchv]$label.degree +
                     jamba::deg2rad(
                        node_adj$label.degree[node_which_l])) %% (pi*2);
            }
         }
      }
      attr(adj_cnet, "nodeset_adj") <- nodeset_adj;
      attr(adj_cnet, "node_adj") <- node_adj;
      assign("adj_cnet",
         value=adj_cnet,
         envir=environment(server));
      # environment(server)$adj_cnet <- adj_cnet;
      return(adj_cnet)
   })

   # get layout coordinate matrix
   layout_coords <- igraph::graph_attr(g, "layout");
   if (length(layout_coords) == 0) {
      g <- multienrichjam::relayout_with_qfr(g,
         repulse=3.7)
      layout_coords <- igraph::graph_attr(g, "layout");
   }
   rownames(layout_coords) <- igraph::V(g)$name;
   colnames(layout_coords)[1:2] <- c("x", "y")

   x_range <- range(layout_coords[, 1], na.rm=TRUE);
   y_range <- range(layout_coords[, 2], na.rm=TRUE);
   max_xy_range <- max(c(
      diff(x_range),
      diff(y_range)));
   jam.vertex.radius <- igraph::V(g)$size * (max_xy_range) / 1.8 / 200;
   names(jam.vertex.radius) <- igraph::V(g)$name;


   ###############################
   # reorder sortAttributes
   output$reorder_output <- shiny::renderPrint({
      input$reorder_attributes
   })


   ###############################
   # igraph plot
   output$networkPlot <- shiny::renderPlot({

      # default behavior
      mark.groups <- NULL;
      mark.expand <- NULL;
      mark.col <- NULL;
      edge_bundling <- "none";
      render_nodelabels <- FALSE;
      use_shadowText <- FALSE;

      # display options
      display_options <- input$display_settings;
      if ("Display Node Labels" %in% display_options) {
         render_nodelabels <- TRUE;
      } else {
         render_nodelabels <- FALSE;
      }
      if ("Use shadowText" %in% display_options) {
         use_shadowText <- TRUE;
      } else {
         use_shadowText <- FALSE;
      }
      if ("Highlight Nodesets" %in% display_options) {
         # mark.groups <- NULL;
      } else {
         # mark.groups <- NULL;
      }
      if ("Bundle Edges" %in% display_options) {
         edge_bundling <- "connections";
      } else {
         edge_bundling <- "none";
      }

      # Logic for highlighted nodeset
      va <- igraph::vertex_attr_names(g);
      use_active_node <- NULL;
      if (!is.null(active_node()) && any(nchar(active_node()) > 0)) {
         use_active_node <- active_node();
         # only for Cnet igraph do we look for nodesets
         if ("nodeType" %in% va) {
            nodesets <- get_cnet_nodeset(g);
            nsdf <- data.frame(node=unlist(nodesets),
               nodeset=rep(names(nodesets), lengths(nodesets)));
            nsdf1 <- subset(nsdf, node %in% use_active_node);
            if (nrow(nsdf1) == 1) {
               nodeset <- nsdf1$nodeset;
               if (length(nodesets[[nodeset]]) > 1) {
                  # update reactive variable active_nodeset
                  active_nodeset(nodeset);
               } else {
                  active_nodeset("");
               }

               mark.groups <- nodesets[nodeset];
               mark.col <- "#FFD7007F";
               if ("Highlight Nodesets" %in% display_options) {
                  nodeset_order <- unique(c(nodeset, names(nodesets)));
                  mark.groups <- nodesets[nodeset_order];
                  mark.col <- rep(c("#FFD7007F", "#AAAAAA44"),
                     c(1, length(nodesets) - 1));
               }
               mark.expand <- NULL;
            } else {
               # update reactive variable active_nodeset
               active_nodeset("");
            }
         }
      } else {
         # update reactive variable active_nodeset
         active_nodeset("");

         if ("Highlight Nodesets" %in% display_options) {
            if ("nodeType" %in% va) {
               mark.groups <- get_cnet_nodeset(g);
            } else {
               mark.groups <- NULL;
            }
         } else {
            mark.groups <- NULL;
         }
         mark.col <- NULL;
      }

      # shadowText options
      options("jam.shadow.n"=32,
         "jam.alphaOutline"=0.2,
         "jam.shadow.r"=0.15);

      # draw the rest of the owl
      # adjusted_cnet()
      use_cnet <- adjusted_cnet()

      # highlight the active node?
      if (length(use_active_node) == 1) {
         #
         if (any(c("pie", "jampie") %in% igraph::V(use_cnet)$shape) &&
               "pie.lwd" %in% igraph::vertex_attr_names(use_cnet)) {
            if ("frame.lwd" %in% igraph::vertex_attr_names(use_cnet)) {
               igraph::V(use_cnet)[use_active_node]$frame.lwd <- (3 +
                     igraph::V(use_cnet)[use_active_node]$frame.lwd);
               igraph::V(use_cnet)$frame.width <- igraph::V(use_cnet)$frame.lwd;
            }
            k <- match(use_active_node, igraph::V(use_cnet)$name);
            pw <- igraph::V(use_cnet)[use_active_node]$pie.lwd;
            igraph::V(use_cnet)[use_active_node]$pie.lwd <- lapply(pw, function(i){
               i + 3;
            })
         } else if ("frame.lwd" %in% igraph::vertex_attr_names(use_cnet)) {
            igraph::V(use_cnet)[use_active_node]$frame.lwd <- (3 +
                  igraph::V(use_cnet)[use_active_node]$frame.lwd);
            igraph::V(use_cnet)$frame.width <- igraph::V(use_cnet)$frame.lwd;
         } else if ("frame.width" %in% igraph::vertex_attr_names(use_cnet)) {
            igraph::V(use_cnet)[use_active_node]$frame.width <- (3 +
                  igraph::V(use_cnet)[use_active_node]$frame.width);
         } else if ("pie.lwd" %in% igraph::vertex_attr_names(use_cnet)) {
            k <- match(use_active_node, igraph::V(use_cnet)$name);
            pw <- igraph::V(use_cnet)[use_active_node]$pie.lwd;
            igraph::V(use_cnet)[use_active_node]$pie.lwd <- lapply(pw, function(i){
               i + 3;
            })
         }
         if (is.na(igraph::V(use_cnet)[use_active_node]$frame.color)) {
            igraph::V(use_cnet)[use_active_node]$frame.color <- "orange";
         }
      }
      use_node_factor <- input$node_factor;
      use_node_factor <- ifelse(use_node_factor >= 0,
         use_node_factor + 1,
         1 / (2^(use_node_factor * -2) / 2))
      use_label_factor <- input$label_factor;
      use_label_factor <- ifelse(use_label_factor >= 0,
         use_label_factor + 1,
         1 / (2^(use_label_factor * -2) / 2))

      withr::with_par(list(mar=c(0.5, 0.5, 0.5, 0.5)), {
         multienrichjam::jam_igraph(
            x=use_cnet,
            # g,
            label_factor=use_label_factor,
            node_factor=use_node_factor,
            mark.groups=mark.groups,
            mark.col=mark.col,
            mark.expand=mark.expand,
            edge_bundling=edge_bundling,
            render_nodelabels=render_nodelabels,
            use_shadowText=use_shadowText)
      })
   })

   ###############################
   # capture mouse click
   output$active_info <- shiny::renderPrint({
      paste0(
         "Node: ", active_node(),
         ", class(environment()$g):", class(environment()$g),
         ", class(environment()$adj_cnet:", class(environment()$adj_cnet),
         ", ls(environment(server)): ", jamba::cPaste(sep=", ", ls(environment(server)))
      )
   });
   output$click_info <- shiny::renderPrint({
      req(input$plot_click)
      click <- input$plot_click

      # get current layout
      layout_coords <- igraph::graph_attr(name="layout",
         shiny::isolate(adjusted_cnet()))

      # Find closest node
      click_distances <- sqrt((layout_coords[,1] - click$x)^2 +
            (layout_coords[,2] - click$y)^2);
      # determine if click is inside each radius
      dist_in_radius <- (click_distances <= jam.vertex.radius)
      if (!any(dist_in_radius)) {
         active_node("")
         return(list(Message='No node clicked, clearing selection.\n'))
      }
      click_distances_in_radius <- click_distances[dist_in_radius];
      closest_node <- which.min(click_distances_in_radius);
      closest_distance <- click_distances_in_radius[closest_node];

      # update reactiveVar for global use
      active_node(head(names(closest_node), 1))

      names(closest_node)
      # remnant debug output
      # list(
      #    Closest_Node=names(closest_node),
      #    # active_node=jamba::rmNULL(nullValue="<NULL>", active_node()),
      #    Clicked_XY=jamba::cPaste(sep=", ",
      #       format(digits=2, c(click$x, click$y))),
      #    Distance_from_Center=closest_distance
      # )
   })

   ######################################
   # Observe changes to active_nodeset()
   # - populate stored values into each input
   observe({
      use_nodeset <- active_nodeset()
      if (is.null(use_nodeset) || any(nchar(use_nodeset) == 0)) {
         return(NULL)
      }
      df <- subset(adjust_nodeset(), nodeset %in% use_nodeset);
      if (nrow(df) == 1) {
         tryCatch({
            shiny::updateNumericInput(
               session=session,
               inputId="nodeset_x",
               value=df$x)
            shiny::updateNumericInput(
               session=session,
               inputId="nodeset_y",
               value=df$y)
            shiny::updateNumericInput(
               session=session,
               inputId="nodeset_rotation",
               value=df$rotate_degrees)
            shiny::updateNumericInput(
               session=session,
               inputId="nodeset_spacing",
               value=df$percent_spacing)
         }, error=function(e){
            # e
         });
      }
   })

   ######################################
   # Observe changes to active_node()
   # - populate stored values into each input
   observe({
      use_node <- active_node()
      if (is.null(use_node) || any(nchar(use_node) == 0)) {
         return(NULL)
      }
      df <- subset(adjust_node(), node %in% use_node);
      if (nrow(df) == 1) {
         tryCatch({
            shiny::updateNumericInput(
               session=session,
               inputId="node_x",
               value=df$x)
            shiny::updateNumericInput(
               session=session,
               inputId="node_y",
               value=df$y)
            shiny::updateNumericInput(
               session=session,
               inputId="label_degrees",
               value=df$label.degree)
            shiny::updateNumericInput(
               session=session,
               inputId="label_distance",
               value=df$label.dist)
         }, error=function(e){
            # e
         })
      }
   })

   ###################################################
   # when active_node is defined, populate UI elements
   # div form-group shiny-input-contained - has margins
   output$node_ui <- shiny::renderUI({
      if (is.null(active_node()) || "" %in% active_node()) {
         return(NULL)
      }
      htmltools::tagList(
         # htmltools::h4(paste0(" Node: ", active_node())),
         htmltools::tags$b("Node: ", style="margin: 0px 0px 0px 10px; font-size: 120%;"),
         htmltools::tags$em(active_node(), style="margin: 0px 0px 0px 10px; font-size: 120%; color: gold; font-weight: bold;"),
         htmltools::br(),
         # div(class="inline-input",
         #    shiny::actionButton("node_x_minus", "<-"),
         #    shiny::actionButton("node_x_plus", "->")),
         div(class="inline-input",
            htmltools::tags$label("x:",
               `for`="node_x",
               style="margin: 0;"),
            shiny::numericInput("node_x",
               label=NULL,
               min=-1, max=1, step=0.005, value=0)
         ),
         div(class="inline-input",
            htmltools::tags$label("y:",
               `for`="node_y",
               style="margin: 0;"),
            shiny::numericInput("node_y",
               label=NULL,
               min=-1, max=1, step=0.005, value=0)
         ),
         div(class="inline-input",
            htmltools::tags$label("Label distance:",
               `for`="label_distance",
               style="margin: 0;"),
            shiny::numericInput("label_distance",
               label=NULL,
               min=-5, max=5, step=0.2, value=0)
         ),
         div(class="inline-input",
            htmltools::tags$label("Label angle:",
               `for`="label_degrees",
               style="margin: 0;"),
            shiny::numericInput("label_degrees",
               label=NULL,
               min=-180, max=180, step=10, value=0)
         )
         # shiny::numericInput("label_distance",
         #    # updateOn="blur",
         #    label="Label distance:",
         #    min=-5, max=5, step=0.2, value=0),
         # shiny::numericInput("label_degrees",
         #    # updateOn="blur",
         #    label="Label angle:",
         #    min=-180, max=180, step=10, value=0)
      )
   })
   ######################################################
   # when active_nodeset is defined, populate UI elements
   output$nodeset_ui <- shiny::renderUI({
      if (is.null(active_nodeset()) || "" %in% active_nodeset()) {
         return(NULL)
      }
      htmltools::tagList(
         # htmltools::h4(paste0(htmltools::tags$b("Nodeset: "), active_nodeset()),
         htmltools::tags$b("Nodeset: ", style="margin: 0px 0px 0px 10px; font-size: 120%;"),
         htmltools::tags$em(active_nodeset(), style="margin: 0px 0px 0px 10px; font-size: 120%; color: gold; font-weight: bold;"),
         htmltools::br(),
         # shiny::verbatimTextOutput("selected_nodeset_info"),
         div(class="inline-input",
            htmltools::tags$label("x:",
               `for`="nodeset_x",
               style="margin: 0;"),
            shiny::numericInput("nodeset_x",
               label=NULL,
               min=-1, max=1, step=0.005, value=0)
         ),
         div(class="inline-input",
            htmltools::tags$label("y:",
               `for`="nodeset_y",
               style="margin: 0;"),
            shiny::numericInput("nodeset_y",
               label=NULL,
               min=-1, max=1, step=0.005, value=0)
         ),
         div(class="inline-input",
            htmltools::tags$label("Spacing:",
               `for`="nodeset_spacing",
               style="margin: 0;"),
            shiny::numericInput("nodeset_spacing",
               label=NULL,
               min=1, max=10, step=0.2, value=3)
         ),
         div(class="inline-input",
            htmltools::tags$label("Spacing:",
               `for`="nodeset_spacing",
               style="margin: 0;"),
            shiny::numericInput("nodeset_spacing",
               label=NULL,
               min=-180, max=180, step=5, value=0)
         )
         # shiny::numericInput("nodeset_x",
         #    # updateOn="blur",
         #    label="Adjust X:",
         #    min=-1, max=1, step=0.005, value=0),
         # shiny::numericInput("nodeset_y",
         #    # updateOn="blur",
         #    label="Adjust Y:",
         #    min=-1, max=1, step=0.005, value=0),
         # shiny::numericInput("nodeset_spacing",
         #    # updateOn="blur",
         #    label="Nodeset Spacing (%):",
         #    min=1, max=10, step=0.2, value=3),
         # shiny::numericInput("nodeset_rotation",
         #    # updateOn="blur",
         #    label="Nodeset Rotation:",
         #    min=-180, max=180, step=5, value=0)
      )
   })

   ######################################
   # observe node adjustments
   shiny::observeEvent(input$node_x, {
      if (is.null(active_node()) || "" %in% active_node()) return(NULL);
      use_node <- active_node();
      df <- adjust_node();
      k <- match(use_node, df$node);
      if (length(k) == 1) {
         df$x[k] <- input$node_x;
         adjust_node(df);
      }
   })
   shiny::observeEvent(input$node_y, {
      if (is.null(active_node()) || "" %in% active_node()) return(NULL);
      use_node <- active_node();
      df <- adjust_node();
      k <- match(use_node, df$node);
      if (length(k) == 1) {
         df$y[k] <- input$node_y;
         adjust_node(df);
      }
   })
   shiny::observeEvent(input$label_degrees, {
      if (is.null(active_node()) || "" %in% active_node()) return(NULL);
      use_node <- active_node();
      df <- adjust_node();
      k <- match(use_node, df$node);
      if (length(k) == 1) {
         df$label.degree[k] <- input$label_degrees;
         adjust_node(df);
      }
   })
   shiny::observeEvent(input$label_distance, {
      if (is.null(active_node()) || "" %in% active_node()) return(NULL);
      use_node <- active_node();
      df <- adjust_node();
      k <- match(use_node, df$node);
      if (length(k) == 1) {
         df$label.dist[k] <- input$label_distance;
         adjust_node(df);
      }
   })

   ######################################
   # observe nodeset adjustments
   shiny::observeEvent(input$nodeset_x, {
      if (is.null(active_nodeset()) || "" %in% active_nodeset()) return(NULL);
      use_nodeset <- active_nodeset();
      df <- adjust_nodeset();
      k <- match(use_nodeset, df$nodeset);
      if (length(k) == 1) {
         df$x[k] <- input$nodeset_x;
         adjust_nodeset(df);
      }
   })
   shiny::observeEvent(input$nodeset_y, {
      if (is.null(active_nodeset()) || "" %in% active_nodeset()) return(NULL);
      use_nodeset <- active_nodeset();
      df <- adjust_nodeset();
      k <- match(use_nodeset, df$nodeset);
      if (length(k) == 1) {
         df$y[k] <- input$nodeset_y;
         adjust_nodeset(df);
      }
   })
   shiny::observeEvent(input$nodeset_spacing, {
      if (is.null(active_nodeset()) || "" %in% active_nodeset()) return(NULL);
      use_nodeset <- active_nodeset();
      df <- adjust_nodeset();
      k <- match(use_nodeset, df$nodeset);
      if (length(k) == 1) {
         df$percent_spacing[k] <- input$nodeset_spacing;
         adjust_nodeset(df);
      }
   })
   shiny::observeEvent(input$nodeset_rotation, {
      if (is.null(active_nodeset()) || "" %in% active_nodeset()) return(NULL);
      use_nodeset <- active_nodeset();
      df <- adjust_nodeset();
      k <- match(use_nodeset, df$nodeset);
      if (length(k) == 1) {
         df$rotate_degrees[k] <- input$nodeset_rotation;
         adjust_nodeset(df);
      }
   })

}
