library(shiny)
library(shinyjs)
library(BioIndex)

invisible(options(shiny.maxRequestSize = 500 * 1024^2))

# Safely retrieve aggregate_gsas whether it is exported or internal in BioIndex
aggregate_gsas_safe <- function(...) {
  if ("aggregate_gsas" %in% getNamespaceExports("BioIndex")) {
    BioIndex::aggregate_gsas(...)
  } else {
    get("aggregate_gsas", envir = asNamespace("BioIndex"))(...)
  }
}

# Normalize the output structure returned by aggregate_gsas
parse_aggregate_result <- function(x, gsas) {
  agg_gsa <- paste(sort(unique(as.character(gsas))), collapse = "")

  if (is.null(x) || !is.list(x) || length(x) < 5) {
    stop("aggregate_gsas returned an unexpected object.")
  }

  nms <- names(x)

  if (!is.null(nms) && all(c("ta", "tb", "tc") %in% nms)) {
    ta_out <- x$ta
    tb_out <- x$tb
    tc_out <- x$tc
    stratification_out <- if ("stratification" %in% nms) x$stratification else x[[4]]
    strata_scheme_out <- if ("strata_scheme" %in% nms) x$strata_scheme else x[[5]]
    aggregated_gsa_out <- if ("aggregated_gsa" %in% nms) x$aggregated_gsa else agg_gsa
  } else {
    ta_out <- x[[1]]
    tb_out <- x[[2]]
    tc_out <- x[[3]]
    stratification_out <- x[[4]]
    strata_scheme_out <- x[[5]]
    aggregated_gsa_out <- agg_gsa
  }

  list(
    ta = ta_out,
    tb = tb_out,
    tc = tc_out,
    stratification = stratification_out,
    strata_scheme = strata_scheme_out,
    aggregated_gsa = aggregated_gsa_out
  )
}

ui <- fluidPage(
  useShinyjs(),

  tags$head(
    tags$style(HTML("
      body {
        background-color: #f7f7f7;
      }

      .app-header {
        text-align: center;
        padding: 15px 0 10px 0;
      }

      .app-subtitle {
        color: #444444;
        margin-top: -5px;
        margin-bottom: 20px;
      }

      .section-box {
        background-color: white;
        border: 1px solid #dddddd;
        border-radius: 10px;
        padding: 20px;
        margin-bottom: 20px;
        box-shadow: 0 1px 4px rgba(0, 0, 0, 0.05);
      }

      .section-title {
        font-size: 22px;
        font-weight: 600;
        margin-bottom: 15px;
        color: #1f2d3d;
      }

      .btn-main {
        min-width: 180px;
        font-size: 16px;
        font-weight: 600;
      }

      .log-box pre {
        background-color: #fafafa;
        border: 1px solid #d9d9d9;
        border-radius: 6px;
        padding: 12px;
        white-space: pre-wrap;
        word-wrap: break-word;
        max-height: 500px;
        overflow-y: auto;
      }

      .status-box pre {
        background-color: #fafafa;
        border: 1px solid #d9d9d9;
        border-radius: 6px;
        padding: 12px;
        white-space: pre-wrap;
      }
    "))
  ),

  div(
    class = "app-header",
    tags$img(src = "BioIndex_Logo.png", height = "170px"),
    h1("BioIndex"),
    p("Analysis of trawl survey data using MEDITS file format", class = "app-subtitle")
  ),

  fluidRow(
    column(
      width = 12,
      div(
        class = "section-box",
        div(class = "section-title", "Input files"),
        fluidRow(
          column(
            width = 3,
            fileInput("ta_file", "Upload TA file", accept = c(".csv", ".txt"))
          ),
          column(
            width = 3,
            fileInput("tb_file", "Upload TB file", accept = c(".csv", ".txt"))
          ),
          column(
            width = 3,
            fileInput("tc_file", "Upload TC file", accept = c(".csv", ".txt"))
          ),
          column(
            width = 3,
            selectInput(
              "sep",
              "Main file separator",
              choices = c(";" = ";", "," = ","),
              selected = ";"
            )
          )
        ),
        br(),
        fluidRow(
          column(
            width = 4,
            fileInput("strata_file", "Optional custom strata CSV", accept = c(".csv", ".txt"))
          ),
          column(
            width = 4,
            fileInput("stratification_file", "Optional custom stratification CSV", accept = c(".csv", ".txt"))
          ),
          column(
            width = 4,
            selectInput(
              "aux_sep",
              "Custom table separator",
              choices = c(";" = ";", "," = ","),
              selected = ";"
            )
          )
        )
      )
    )
  ),

  fluidRow(
    column(
      width = 12,
      div(
        class = "section-box",
        div(class = "section-title", "Analysis settings"),

        fluidRow(
          column(
            width = 3,
            selectizeInput(
              "sspp",
              "Species code (sspp)",
              choices = NULL,
              selected = NULL,
              options = list(placeholder = "Upload TA/TB/TC files first")
            )
          ),
          column(
            width = 3,
            checkboxInput(
              "use_gsa_aggregation",
              "Aggregate multiple GSAs",
              value = FALSE
            ),
            uiOutput("gsa_selector_ui")
          ),
          column(
            width = 3,
            selectInput(
              "country",
              "Country",
              choices = NULL,
              selected = NULL
            )
          ),
          column(
            width = 3,
            selectInput(
              "sexes",
              "Sexes",
              choices = c("all", "M", "F"),
              selected = "all"
            )
          )
        ),

        fluidRow(
          column(
            width = 3,
            numericInput("rec_threshold", "Recruit threshold (mm)", value = 200, min = 0)
          ),
          column(
            width = 3,
            numericInput("spaw_threshold", "Spawner threshold (mm)", value = 210, min = 0)
          ),
          column(
            width = 3,
            numericInput("haul_threshold", "Haul threshold", value = 30, min = 1, step = 1)
          ),
          column(
            width = 3,
            numericInput("buffer", "Buffer", value = 0.1, min = 0, step = 0.01)
          )
        ),

        fluidRow(
          column(
            width = 3,
            selectInput("depth_min", "Depth lower", choices = NULL)
          ),
          column(
            width = 3,
            selectInput("depth_max", "Depth upper", choices = NULL)
          ),
          column(
            width = 3,
            textInput("depth_lines", "Depth lines", value = "")
          ),
          column(
            width = 3,
            numericInput("resolution", "Resolution", value = NA)
          )
        ),

        fluidRow(
          column(
            width = 3,
            numericInput("xmin", "Longitude lower limit", value = NA)
          ),
          column(
            width = 3,
            numericInput("xmax", "Longitude upper limit", value = NA)
          ),
          column(
            width = 3,
            numericInput("ymin", "Latitude lower limit", value = NA)
          ),
          column(
            width = 3,
            numericInput("ymax", "Latitude upper limit", value = NA)
          )
        ),

        fluidRow(
          style = "margin-top: 10px;",
          column(
            width = 12,
            div(
              style = "text-align: center;",
              checkboxInput("verbose", "Verbose", value = TRUE)
            )
          )
        ),

        fluidRow(
          column(
            width = 12,
            style = "text-align: center; margin-top: 10px;",
            actionButton("run_btn", "Run BioIndex", class = "btn-primary btn-main"),
            tags$span(style = "display:inline-block; width:15px;"),
            downloadButton("download_results", "Download results zip", class = "btn-default btn-main")
          )
        )
      )
    )
  ),

  fluidRow(
    column(
      width = 12,
      div(
        class = "section-box status-box",
        div(class = "section-title", "Status"),
        verbatimTextOutput("status")
      )
    )
  ),

  fluidRow(
    column(
      width = 12,
      div(
        class = "section-box log-box",
        div(class = "section-title", "Run log"),
        verbatimTextOutput("run_log")
      )
    )
  )
)

server <- function(input, output, session) {

  gsa_defaults <- data.frame(
    GSA  = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
             16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
    xmin = c(-5.600000, -3.333333, -5.600000, -2.216667,  0.500000,
             -1.000000,  2.951167,  7.627755,  7.525000, 11.000000,
             6.000000,  8.500000, 10.467833, 10.015361, 13.500000,
             11.000000, 12.141194, 15.169798, 15.082944, 19.166666,
             11.536992, 22.527028, 23.000000, 29.000000, 32.000000,
             25.139999, 34.349998, 26.324205, 27.306167, 34.465444),
    xmax = c(-1.000000, -2.666667, -2.216667,  8.806623,  6.000000,
             6.000000,  8.000000,  9.750000, 13.004581, 16.221250,
             11.000000, 12.000000, 13.500000, 15.300000, 15.300000,
             15.300000, 18.765417, 20.186167, 19.166666, 23.225389,
             25.139999, 29.000000, 26.500000, 36.209250, 35.000000,
             34.349998, 35.991306, 29.939583, 41.779528, 39.305611),
    ymin = c(36.00000, 35.75000, 35.08453, 35.07286, 38.00000,
             37.00000, 41.33000, 41.30000, 41.30000, 37.03458,
             38.00000, 36.71289, 35.00000, 33.16779, 35.00000,
             35.00000, 41.92224, 39.60821, 35.00000, 35.00000,
             30.26622, 34.00000, 34.00000, 34.00000, 34.00000,
             30.82203, 31.44039, 40.07335, 40.90875, 45.23872),
    ymax = c(37.58622, 36.08333, 36.00000, 38.00000, 40.50000,
             42.42947, 43.78458, 43.25000, 44.42711, 41.30000,
             41.78455, 38.50000, 37.00000, 35.00000, 36.50000,
             38.00000, 45.79456, 42.47188, 40.52042, 39.88215,
             35.00000, 41.01206, 36.00000, 36.91292, 35.78333,
             34.00000, 35.88234, 41.15125, 46.88956, 47.28956)
  )

  get_gsa_map_limits <- function(gsa_values, gsa_defaults) {
    gsa_values <- suppressWarnings(as.numeric(gsa_values))
    gsa_values <- gsa_values[!is.na(gsa_values)]

    if (length(gsa_values) == 0) {
      return(NULL)
    }

    rows <- gsa_defaults[gsa_defaults$GSA %in% gsa_values, , drop = FALSE]

    if (nrow(rows) == 0) {
      return(NULL)
    }

    c(
      xmin = min(rows$xmin, na.rm = TRUE),
      xmax = max(rows$xmax, na.rm = TRUE),
      ymin = min(rows$ymin, na.rm = TRUE),
      ymax = max(rows$ymax, na.rm = TRUE)
    )
  }

  rv <- reactiveValues(
    status = "Waiting for input.",
    log = "",
    zipfile = NULL,
    run_dir = NULL,
    running = FALSE,
    ta = NULL,
    tb = NULL,
    tc = NULL,
    strata = BioIndex::strata_scheme,
    stratification_tab = BioIndex::stratification
  )

  output$status <- renderText({
    rv$status
  })

  output$run_log <- renderText({
    rv$log
  })

  output$gsa_selector_ui <- renderUI({
    req(rv$ta)

    gsa_choices <- sort(unique(suppressWarnings(as.numeric(as.character(rv$ta$AREA)))))
    gsa_choices <- gsa_choices[!is.na(gsa_choices)]

    if (isTRUE(input$use_gsa_aggregation)) {
      selectizeInput(
        "GSA_multi",
        "GSAs to aggregate",
        choices = gsa_choices,
        selected = NULL,
        multiple = TRUE,
        options = list(
          placeholder = "Select two or more GSAs",
          plugins = list("remove_button")
        )
      )
    } else {
      selectInput(
        "GSA",
        "GSA",
        choices = gsa_choices,
        selected = if (length(gsa_choices) > 0) gsa_choices[1] else NULL
      )
    }
  })

  read_medits_file <- function(fileinfo, sep) {
    req(fileinfo)
    read.table(
      fileinfo$datapath,
      sep = sep,
      header = TRUE,
      stringsAsFactors = FALSE
    )
  }

  observe({
    if (isTRUE(rv$running)) {
      disable("run_btn")
    } else {
      enable("run_btn")
    }
  })

  observe({
    req(input$strata_file)

    strata_custom <- try(read_medits_file(input$strata_file, input$aux_sep), silent = TRUE)

    if (inherits(strata_custom, "try-error")) {
      rv$status <- "Error while reading custom strata file. Default strata will be used."
      rv$strata <- BioIndex::strata_scheme
    } else {
      rv$strata <- strata_custom
      rv$status <- "Custom strata file loaded."
    }
  })

  observe({
    if (is.null(input$strata_file)) {
      rv$strata <- BioIndex::strata_scheme
    }
  })

  observe({
    req(input$stratification_file)

    strat_custom <- try(read_medits_file(input$stratification_file, input$aux_sep), silent = TRUE)

    if (inherits(strat_custom, "try-error")) {
      rv$status <- "Error while reading custom stratification file. Default stratification will be used."
      rv$stratification_tab <- BioIndex::stratification
    } else {
      rv$stratification_tab <- strat_custom
      rv$status <- "Custom stratification file loaded."
    }
  })

  observe({
    if (is.null(input$stratification_file)) {
      rv$stratification_tab <- BioIndex::stratification
    }
  })

  observe({
    req(input$ta_file, input$tb_file, input$tc_file)

    ta <- try(read_medits_file(input$ta_file, input$sep), silent = TRUE)
    tb <- try(read_medits_file(input$tb_file, input$sep), silent = TRUE)
    tc <- try(read_medits_file(input$tc_file, input$sep), silent = TRUE)

    if (inherits(ta, "try-error") || inherits(tb, "try-error") || inherits(tc, "try-error")) {
      rv$status <- "Error while reading one or more input files for selector update."
      return(NULL)
    }

    rv$ta <- ta
    rv$tb <- tb
    rv$tc <- tc

    if (!all(c("GENUS", "SPECIES") %in% names(tb))) {
      rv$status <- "TB file must contain GENUS and SPECIES columns."
      return(NULL)
    }

    if (!"AREA" %in% names(ta)) {
      rv$status <- "TA file must contain AREA column."
      return(NULL)
    }

    if (!"COUNTRY" %in% names(ta)) {
      rv$status <- "TA file must contain COUNTRY column."
      return(NULL)
    }

    sspp_choices <- sort(unique(paste0(
      trimws(as.character(tb$GENUS)),
      trimws(as.character(tb$SPECIES))
    )))
    sspp_choices <- sspp_choices[!is.na(sspp_choices) & nzchar(sspp_choices)]

    COUNTRY_choices <- sort(unique(trimws(as.character(ta$COUNTRY))))
    COUNTRY_choices <- COUNTRY_choices[!is.na(COUNTRY_choices) & nzchar(COUNTRY_choices)]

    updateSelectizeInput(
      session,
      "sspp",
      choices = sspp_choices,
      selected = if (length(sspp_choices) > 0) sspp_choices[1] else NULL,
      server = TRUE
    )

    updateSelectInput(
      session,
      "country",
      choices = c("all", COUNTRY_choices),
      selected = "all"
    )

    rv$status <- "Files loaded. Select species, GSA and country."
  })

  observe({
    selected_gsas <- if (isTRUE(input$use_gsa_aggregation)) input$GSA_multi else input$GSA
    selected_country <- input$country

    req(selected_gsas)

    # Update map limits
    map_vals <- get_gsa_map_limits(selected_gsas, gsa_defaults)

    if (!is.null(map_vals)) {
      updateNumericInput(session, "xmin", value = unname(map_vals["xmin"]))
      updateNumericInput(session, "xmax", value = unname(map_vals["xmax"]))
      updateNumericInput(session, "ymin", value = unname(map_vals["ymin"]))
      updateNumericInput(session, "ymax", value = unname(map_vals["ymax"]))
    } else {
      updateNumericInput(session, "xmin", value = NA_real_)
      updateNumericInput(session, "xmax", value = NA_real_)
      updateNumericInput(session, "ymin", value = NA_real_)
      updateNumericInput(session, "ymax", value = NA_real_)
    }

    # Update depth range based on strata
    current_strata <- rv$strata
    if (!is.null(current_strata) && "GSA" %in% names(current_strata)) {
      relevant_strata <- current_strata[as.character(current_strata$GSA) %in% as.character(selected_gsas), ]

      if (!is.null(selected_country) && selected_country != "all" && "COUNTRY" %in% names(relevant_strata)) {
        relevant_strata <- relevant_strata[relevant_strata$COUNTRY == selected_country, ]
      }

      if (nrow(relevant_strata) > 0) {
        all_mins <- sort(unique(relevant_strata$MIN_DEPTH))
        all_maxs <- sort(unique(relevant_strata$MAX_DEPTH))

        updateSelectInput(session, "depth_min", choices = all_mins, selected = min(all_mins))
        updateSelectInput(session, "depth_max", choices = all_maxs, selected = max(all_maxs))

        # Update depth lines with sensible defaults from stratification cutoffs
        all_depths <- sort(unique(c(all_mins, all_maxs)))
        d_min <- min(all_depths)
        d_max <- max(all_depths)
        
        if (length(all_depths) >= 3) {
          inner_depths <- all_depths[-c(1, length(all_depths))]
          target <- (d_min + d_max) / 2
          mid_depth <- inner_depths[which.min(abs(inner_depths - target))]
        } else {
          mid_depth <- round((d_min + d_max) / 2)
        }
        
        updateTextInput(session, "depth_lines", value = paste(d_min, mid_depth, d_max, sep = ","))
      }
    }
  })

  observeEvent(input$run_btn, {

    req(rv$ta, rv$tb, rv$tc)
    req(input$sspp, input$country)

    sspp_val <- as.character(input$sspp)
    rec_threshold_val <- suppressWarnings(as.numeric(input$rec_threshold))
    spaw_threshold_val <- suppressWarnings(as.numeric(input$spaw_threshold))
    haul_threshold_val <- suppressWarnings(as.integer(input$haul_threshold))

    resolution_val <- input$resolution
    if (length(resolution_val) == 0 || is.null(resolution_val) || identical(resolution_val, "")) {
      resolution_val <- NA_real_
    } else {
      resolution_val <- suppressWarnings(as.numeric(resolution_val))
    }

    if (is.na(resolution_val) && !is.na(input$resolution)) {
      rv$status <- "Invalid resolution."
      rv$log <- "Resolution must be numeric or empty."
      rv$running <- FALSE
      return(NULL)
    }

    buffer_val <- suppressWarnings(as.numeric(input$buffer))
    sexes_val <- as.character(input$sexes)
    single_gsa_val <- input$GSA
    multi_gsa_val <- input$GSA_multi
    country_val <- as.character(input$country)
    strata_val <- rv$strata
    stratification_tab_val <- rv$stratification_tab
    zip_val <- TRUE
    save_val <- TRUE
    verbose_val <- isTRUE(input$verbose)

    depth_val <- c(
      suppressWarnings(as.numeric(input$depth_min)),
      suppressWarnings(as.numeric(input$depth_max))
    )

    map_lim_val <- c(
      suppressWarnings(as.numeric(input$xmin)),
      suppressWarnings(as.numeric(input$xmax)),
      suppressWarnings(as.numeric(input$ymin)),
      suppressWarnings(as.numeric(input$ymax))
    )

    if (any(is.na(depth_val))) {
      rv$status <- "Invalid depth range."
      rv$log <- "Depth lower and upper must be numeric."
      rv$running <- FALSE
      return(NULL)
    }

    if (any(is.na(map_lim_val))) {
      rv$status <- "Invalid map limits."
      rv$log <- "Map limits are missing. Select a GSA with defaults or insert xmin, xmax, ymin, ymax manually."
      rv$running <- FALSE
      return(NULL)
    }

    depth_lines_val <- suppressWarnings(as.numeric(trimws(strsplit(input$depth_lines, ",")[[1]])))

    if (any(is.na(depth_lines_val))) {
      rv$status <- "Invalid depth lines."
      rv$log <- "Depth lines must be provided as comma-separated numeric values."
      rv$running <- FALSE
      return(NULL)
    }

    if (is.na(rec_threshold_val)) {
      rv$status <- "Invalid recruit threshold."
      rv$log <- "Recruit threshold must be numeric."
      rv$running <- FALSE
      return(NULL)
    }

    if (is.na(spaw_threshold_val)) {
      rv$status <- "Invalid spawner threshold."
      rv$log <- "Spawner threshold must be numeric."
      rv$running <- FALSE
      return(NULL)
    }

    if (!isTRUE(all.equal(as.numeric(input$haul_threshold), as.numeric(haul_threshold_val)))) {
      rv$status <- "Invalid haul threshold."
      rv$log <- "Haul threshold must be an integer value."
      rv$running <- FALSE
      return(NULL)
    }

    if (is.na(buffer_val)) {
      rv$status <- "Invalid buffer."
      rv$log <- "Buffer must be numeric."
      rv$running <- FALSE
      return(NULL)
    }

    if (depth_val[1] >= depth_val[2]) {
      rv$status <- "Invalid depth range."
      rv$log <- "Depth lower must be smaller than depth upper."
      rv$running <- FALSE
      return(NULL)
    }

    if (map_lim_val[1] >= map_lim_val[2] || map_lim_val[3] >= map_lim_val[4]) {
      rv$status <- "Invalid map limits."
      rv$log <- "Longitude and latitude lower limits must be smaller than upper limits."
      rv$running <- FALSE
      return(NULL)
    }

    rv$status <- "Preparing analysis..."
    rv$log <- ""
    rv$zipfile <- NULL
    rv$run_dir <- NULL
    rv$running <- TRUE

    ta <- rv$ta
    tb <- rv$tb
    tc <- rv$tc

    if (isTRUE(input$use_gsa_aggregation)) {

      gsas_run <- sort(unique(as.character(multi_gsa_val)))

      if (length(gsas_run) < 2) {
        rv$status <- "Invalid aggregated GSA selection."
        rv$log <- "Please select at least two GSAs to aggregate."
        rv$running <- FALSE
        return(NULL)
      }

      agg_raw <- tryCatch(
        aggregate_gsas_safe(
          ta = ta,
          tb = tb,
          tc = tc,
          gsas = gsas_run,
          strata_scheme = strata_val,
          stratification = stratification_tab_val
        ),
        error = function(e) e
      )

      if (inherits(agg_raw, "error")) {
        rv$status <- "Aggregation failed."
        rv$log <- paste("ERROR:", conditionMessage(agg_raw))
        rv$running <- FALSE
        return(NULL)
      }

      agg_res <- tryCatch(
        parse_aggregate_result(agg_raw, gsas_run),
        error = function(e) e
      )

      if (inherits(agg_res, "error")) {
        rv$status <- "Aggregation failed."
        rv$log <- paste("ERROR:", conditionMessage(agg_res))
        rv$running <- FALSE
        return(NULL)
      }

      ta <- agg_res$ta
      tb <- agg_res$tb
      tc <- agg_res$tc
      stratification_tab_val <- agg_res$stratification
      strata_val <- agg_res$strata_scheme
      GSA_val <- agg_res$aggregated_gsa

    } else {

      req(single_gsa_val)

      GSA_val <- suppressWarnings(as.numeric(single_gsa_val))

      if (is.na(GSA_val)) {
        rv$status <- "Invalid GSA."
        rv$log <- "GSA must be numeric."
        rv$running <- FALSE
        return(NULL)
      }
    }

    wd_run <- file.path(
      tempdir(),
      paste0("BioIndex_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    )
    dir.create(wd_run, recursive = TRUE, showWarnings = FALSE)
    rv$run_dir <- wd_run

    rv$status <- "BioIndex analysis started..."

    result <- tryCatch({
      log_output <- capture.output(
        BioIndex(
          ta = ta,
          tb = tb,
          tc = tc,
          sspp = sspp_val,
          rec_threshold = rec_threshold_val,
          spaw_threshold = spaw_threshold_val,
          haul_threshold = haul_threshold_val,
          sexes = sexes_val,
          depth = depth_val,
          GSA = GSA_val,
          country = country_val,
          map_lim = map_lim_val,
          depth_lines = depth_lines_val,
          strata = strata_val,
          stratification_tab = stratification_tab_val,
          resolution = resolution_val,
          buffer = buffer_val,
          wd = wd_run,
          zip = zip_val,
          save = save_val,
          verbose = verbose_val
        ),
        type = "output"
      )

      list(
        success = TRUE,
        log_output = log_output
      )
    }, error = function(e) {
      list(
        success = FALSE,
        error_message = conditionMessage(e)
      )
    })

    if (isTRUE(result$success)) {

      zip_candidates <- list.files(
        path = wd_run,
        pattern = "\\.zip$",
        recursive = TRUE,
        full.names = TRUE
      )

      final_log <- result$log_output

      if (length(zip_candidates) > 0) {
        zipfile_found <- zip_candidates[1]
        final_log <- c(
          final_log,
          "",
          "Run completed successfully.",
          paste("ZIP file created:", zipfile_found)
        )
      } else {
        zipfile_found <- NULL
        final_log <- c(
          final_log,
          "",
          "Run completed successfully, but no ZIP file was found."
        )
      }

      rv$zipfile <- zipfile_found
      rv$log <- paste(final_log, collapse = "\n")
      rv$status <- "Analysis completed."

    } else {
      rv$zipfile <- NULL
      rv$log <- paste("ERROR:", result$error_message)
      rv$status <- "Analysis failed."
    }

    rv$running <- FALSE
  })

  output$download_results <- downloadHandler(
    filename = function() {
      zipfile <- rv$zipfile
      if (!is.null(zipfile) && file.exists(zipfile)) {
        basename(zipfile)
      } else {
        "BioIndex_results.zip"
      }
    },
    content = function(file) {
      zipfile <- rv$zipfile
      req(!is.null(zipfile))
      req(file.exists(zipfile))
      ok <- file.copy(zipfile, file, overwrite = TRUE)
      if (!isTRUE(ok)) {
        stop("Unable to copy the ZIP file for download.")
      }
    }
  )
}

shinyApp(ui, server)
