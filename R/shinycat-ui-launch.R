

#' R-shiny ui function for shinycat
#'
#' R-shiny ui function for shinycat, internal use by launch_shinycat()
#'
#' Internal use by `launch_shinycat()`, this function provides the
#' 'ui' component for the R-shiny app.
#'
#' @returns A UI definition that can be passed to `shiny::shinyApp()`
#'    as argument 'ui'.
#'
#' @family jam shiny functions
#'
#' @export
shinycat_ui <- function
(...)
{
   #
   header <- shinydashboard::dashboardHeader(
      title="ShinyCat: Cnet Adjustment Tool")

   sidebar <- shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
         shinydashboard::menuItem("Cnet Adjustment Tool",
            tabName="cnet_tab1",
            icon=shiny::icon("chart-line")),
         shiny::downloadButton("download_graph",
            icon=shiny::icon("download"),
            style="padding: 5px; margin-left: 10px;",
            label="Save as .RData"),
         shiny::textOutput("save_message")
         # shinydashboard::menuItem("Adjustment Table",
         #    tabName="adjustment_tab1",
         #    icon=shiny::icon("table"))
      ),
      htmltools::hr(),
      htmltools::h4("Visual Options"),
      shiny::checkboxGroupInput("display_settings",
         label="Display Settings",
         choices=c("Display Node Labels",
            "Use shadowText",
            "Highlight Nodesets",
            "Bundle Edges"),
         selected=NULL),
      htmltools::h5("Node Factor"),
      div(class="inline-input",
         htmltools::tags$label("Genes:",
            `for`="gene_node_factor",
            style="margin: 0;"),
         shiny::numericInput("gene_node_factor",
            label=NULL,
            min=0.1, max=3, step=0.1, value=1)
      ),
      div(class="inline-input",
         htmltools::tags$label("Sets:",
            `for`="set_node_factor",
            style="margin: 0;"),
         shiny::numericInput("set_node_factor",
            label=NULL,
            min=0.1, max=3, step=0.1, value=1)
      ),
      shiny::conditionalPanel(
         condition='input.display_settings.includes("Display Node Labels")',
         htmltools::h5("Label Factor"),
         div(class="inline-input",
            htmltools::tags$label("Genes:",
               `for`="gene_label_factor",
               style="margin: 0;"),
            shiny::numericInput("gene_label_factor",
               label=NULL,
               min=0.1, max=3, step=0.1, value=1)
         ),
         div(class="inline-input",
            htmltools::tags$label("Sets:",
               `for`="set_label_factor",
               style="margin: 0;"),
            shiny::numericInput("set_label_factor",
               label=NULL,
               min=0.1, max=3, step=0.1, value=1)
         )
      ),
      htmltools::hr(),
      shiny::uiOutput("node_ui"),
      htmltools::hr(),
      shiny::uiOutput("nodeset_ui")
   )
   ## Prototype using left/right, up/down buttons
      # fluidRow(
      #    column(5, ""),
      #    column(2, actionButton("up_btn", label = "\u2191")),  # Up arrow
      #    column(5, "")),
      # fluidRow(
      #    column(4, ""),
      #    column(2, actionButton("left_btn", label = "\u2190")),  # Left arrow
      #    column(2, actionButton("right_btn", label = "\u2192")),  # Right arrow
      #    column(4, "")),
      # fluidRow(
      #    column(5, ""),
      #    column(2, actionButton("down_btn", label = "\u2193")),  # Down arrow
      #    column(5, "")
      # )

   # adjustment tab to consider:
   # - display the nodeset_adj, node_adj data.frame adjustments
   # - button to save the current igraph object to RData file
   adjustment_tab <- shiny::fluidPage(
      shiny::fluidRow(
         shiny::column(
            width=12,
            style="padding:0px",
            shinydashboardPlus::box(
               shiny::downloadButton("download_graph",
                  icon=shiny::icon("download"),
                  label="Save as .RData"),
               shiny::textOutput("save_message")
            )
         )
      )
   )

   cnet_tab <- shiny::fluidPage(
      shiny::fluidRow(
         shiny::column(
            width=12,
            style="padding:0px",
            shinydashboardPlus::box(
               title="Cnet Plot",
               status="primary",
               solidHeader=TRUE,
               closable=FALSE,
               collapsible=TRUE,
               width=12,
               height="100%",
               shiny::fluidRow(
                  shiny::column(
                     width=12,
                     style="padding:15px",
                     shiny::plotOutput("networkPlot",
                        width="100%",
                        height="800px",
                        click="plot_click"),
                     shiny::verbatimTextOutput("active_info"),
                     shiny::verbatimTextOutput("click_info")
                  )
               )
            )
         )
      )
   )

   body <- shinydashboard::dashboardBody(
      shinydashboard::tabItems(
         shinydashboard::tabItem(tabName="cnet_tab1", cnet_tab)
         # shinydashboard::tabItem(tabName="adjustment_tab1", adjustment_tab)
      )
   )

   # put it together
   shinydashboard::dashboardPage(
      tags$head(
         tags$style(HTML("
         section.sidebar .shiny-bound-input.action-button {
           margin: 0px;
         }
         .btn {
           padding: 0px;
         }
         hr {
           margin-bottom: 2px;
           margin-top: 2px;
         }
         h4 {
           margin-left: 10px;
           margin-bottom: 2px;
           margin-top: 2px;
         }
         h5 {
           margin-left: 10px;
           margin-bottom: 0px;
           font-weight: bold;
         }
         label {
           margin-bottom: 1px;
         }
         section.sidebar .shiny-input-container {
           padding: 5px 10px;
         }
         .shiny-input-container {
           padding: 5px 10px;
         }
         .form-group {
           margin-bottom: 3px;
         }
         .form-control {
           padding: 3px 6px;
         }
         .inline-input {
           display: flex;
           align-items: center;
           gap: 10px;
           margin-left: 15px;
           margin-bottom: 0px;
         }
         .inline-input label {
           margin: 0;
           white-space: nowrap;
         }
         "))
      ),
      header=header,
      sidebar=sidebar,
      body=body,
      skin="blue")
}


#' Launch ShinyCAT: Cnet Adjustment Tool
#'
#' Launch ShinyCAT: Cnet Adjustment Tool
#'
#' This function launches the R-shiny app 'ShinyCat' to manipulate
#' a Cnet `igraph` object, which is expected to contain node (vertex)
#' attribute 'nodeType' with values 'Gene' and 'Set'.
#'
#' There are two ways to retain the resulting igraph object:
#'
#' 1. Capture the output of this function, an `environment` described below.
#' 2. Click "Adjustments" then "Save RData", which stores the 'igraph`
#' object as 'adj_cnet'.
#'
#' The function returns an `environment` inside which are two objects:
#'
#' * 'g': the `igraph` data input
#' * 'adj_cnet': the adjusted Cnet `igraph` output
#'
#' Additionally, the 'adj_cnet' object contains two attributes:
#'
#' * 'nodeset_adj', 'node_adj': `data.frame` objects used by
#' `bulk_cnet_adjustments()` to convert 'g' into 'adj_cnet'.
#'
#'
#' @param g `igraph` object, or can be NULL if the variable 'g' is
#'    defined in 'envir'.
#' @param envir `environment` assigned to the R-shiny function space,
#'    which may be useful when capturing the `igraph` object after
#'    manipulation via the R-shiny app.
#' @param ... additional arguments are ignored.
#' @param options `list` with additional settings, for example:
#'    * host: `character` string with hostname or IP address
#'    to permit incoming connections.
#'    Use '0.0.0.0' to accept all incoming host names.
#'    * port: `numerical` port to listen for incoming connections.
#'    * quiet: `logical` whether to
#'    * launch.browser: `logical` whether to launch a web browser, default
#'    TRUE for interactive sessions.
#'    * display.mode: `character` string with display mode: 'normal' is
#'    recommended; 'showcase' displays the code beside the app.
#'
#' @returns `environment` invisibly, containing
#'    * 'g' the input `igraph`
#'    * 'adj_cnet' the output `igraph` after manipulation in the R-shiny app.
#'    It also contains two attributes for reference:
#'
#'       1. 'nodeset_adj': `data.frame` used to adjust nodesets.
#'       2. 'node_adj': `data.frame` used to adjust nodes.
#'
#'
#' @family jam shiny functions
#'
#' @examples
#' # create Cnet test data
#' g <- make_cnet_test();
#' cnetenv <- new.env();
#'
#' # you must catch the output to use the resulting igraph object
#' output_envir <- launch_shinycat();
#'
#' @export
launch_shinycat <- function
(g=NULL,
 envir=new.env(),
 ...,
 options=list(width=1200))
{
   # Todo:
   # - mechanism to pass igraph object, and '...'
   # - consider new.env() and assigning variables into it
   # - then functions should pass those variables as if they
   #   were passing from '...' (which they are not)
   #

   ## validate that envir contains 'g'
   if (!inherits(envir, "environment")) {
      envir <- new.env();
   }
   ## validate 'g'
   if (length(g) > 0) {
      assign("g", g, envir=envir)
   } else {
      if (!exists("g", envir=envir)) {
         stop(paste0("'g' must be provided either directly, or ",
            "via environment 'envir'"));
      }
   }

   # set this function as parent.env()
   original_parent <- parent.env(envir);
   on.exit({
      parent.env(envir) <- original_parent;
   }, add=TRUE)
   # temporarily set parent.env to this function
   # to make 'environment(server)' work as expected
   parent.env(envir) <- environment();

   ui <- shinycat_ui;
   server <- shinycat_server;
   environment(ui) <- envir;
   environment(server) <- envir;

   # Create the app
   shinycat_app <- shiny::shinyApp(ui=ui,
      server=server,
      options=options);

   # Run the app
   retval <- shiny::runApp(shinycat_app,
      ...)
   return(envir)
}
