# Jensen–Shannon Divergence Playground (Multinomial)
# --------------------------------------------------
# Requirements: shiny, plotly, dplyr
# Run with:
#   source("jsd_multinomial_shiny_app.R")
# or shiny::runApp("jsd_multinomial_shiny_app.R")

suppressPackageStartupMessages({
  ok <- requireNamespace("shiny", quietly = TRUE); if (!ok) install.packages("shiny")
  ok <- requireNamespace("plotly", quietly = TRUE); if (!ok) install.packages("plotly")
  ok <- requireNamespace("dplyr", quietly = TRUE); if (!ok) install.packages("dplyr")
})

library(shiny)
library(plotly)
library(dplyr)

# -------- Helpers --------

# Dirichlet sampler for a random PMF of length k
rdirichlet1 <- function(k, alpha = NULL) {
  if (is.null(alpha)) alpha <- runif(k, 0.5, 2)
  x <- rgamma(k, shape = alpha, rate = 1)
  x / sum(x)
}

# Adjust guessed vector q by setting q[i] to new_i and distributing the delta
# equally among the other components, with epsilon floors and renormalization
adjust_q <- function(q, i, new_i, eps = 1e-9) {
  k <- length(q)
  if (i < 1 || i > k) return(q)
  new_i <- max(min(new_i, 1 - eps), eps)
  delta <- new_i - q[i]
  if (abs(delta) < 1e-15) return(q)
  others <- setdiff(seq_len(k), i)

  if (delta > 0) {
    available <- sum(pmax(q[others] - eps, 0))
    delta_eff <- min(delta, available)
    if (delta_eff <= 0) return(q)
    q[others] <- pmax(q[others] - delta_eff/length(others), eps)
    q[i] <- q[i] + delta_eff
  } else {
    inc_total <- -delta
    capacity <- sum(pmax((1 - eps) - q[others], 0))
    inc_eff <- min(inc_total, capacity)
    if (inc_eff <= 0) return(q)
    q[others] <- pmin(q[others] + inc_eff/length(others), 1 - eps)
    q[i] <- q[i] - inc_eff
  }
  q <- pmax(q, eps)
  q / sum(q)
}

# Entropy with base option (nats/bits)
entropy <- function(p, base = "nats", eps = 1e-12) {
  p <- pmax(p, eps)
  if (base == "bits") return(-sum(p * log2(p)))
  -sum(p * log(p))
}

# KL(P||Q) with base option
kl_div <- function(p, q, base = "nats", eps = 1e-12) {
  p <- pmax(p, eps); q <- pmax(q, eps)
  if (base == "bits") return(sum(p * (log2(p) - log2(q))))
  sum(p * (log(p) - log(q)))
}

# Jensen–Shannon divergence and distance
# JSD(P||Q) = 0.5*KL(P||M) + 0.5*KL(Q||M), M=(P+Q)/2
# sqrt(JSD) is a metric (the Jensen–Shannon distance)
js_divergence <- function(p, q, base = "nats", eps = 1e-12) {
  m <- 0.5*(p + q)
  0.5*kl_div(p, m, base = base, eps = eps) + 0.5*kl_div(q, m, base = base, eps = eps)
}

fmt_pct <- function(x) sprintf('%.2f%%', 100 * x)
cat_names <- function(k) paste0("C", seq_len(k))

# -------- UI --------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .metric-card {padding: 10px 12px; border-radius: 10px; background: #f6f7fb; margin-top: 8px;}
      .metric-title {font-weight: 600; color: #444; font-size: 13px;}
      .metric-value {font-weight: 700; font-size: 18px; color: #111;}
      .subtle {color: #777; font-size: 12px;}
      .selected-note {font-size: 13px; color: #555;}
    "))
  ),
  titlePanel('Jensen–Shannon Divergence (Multinomial)'),
  sidebarLayout(
    sidebarPanel(
      sliderInput("k", "Number of categories (k)", min = 2, max = 10, value = 5, step = 1),
      radioButtons("base", "Units", choices = c("nats", "bits"), selected = "nats", inline = TRUE),
      actionButton("newP", "NEW (sample TRUE P)", class = "btn btn-primary btn-sm"),
      tags$hr(),
      uiOutput("selectedInfo"),
      sliderInput("adjustVal", "Adjust selected guessed bar", min = 0, max = 1, value = 0, step = 0.001),
      fluidRow(
        column(6, actionButton("decBtn", "– 0.01", class = "btn btn-outline-secondary btn-sm", width = "100%")),
        column(6, actionButton("incBtn", "+ 0.01", class = "btn btn-outline-secondary btn-sm", width = "100%"))
      ),
      helpText("Click a RED bar to select. Use the slider or buttons to adjust. (Plotly bars are not directly draggable.)")
    ),
    mainPanel(
      plotlyOutput("plotPQ", height = "360px"),
      div(class = "metric-card",
          div(class = "metric-title", "Information Metrics"),
          uiOutput("metricsUI"),
          div(class = "subtle", textOutput("boundsText"))
      ),
      div(class = "metric-card",
          div(class = "metric-title", "Values"),
          tableOutput("pqTable")
      )
    )
  )
)

# -------- Server --------
server <- function(input, output, session) {
  eps <- 1e-9
  rv <- reactiveValues(P=NULL, Q=NULL, sel=NULL)

  # Initialize and NEW
  init_reset <- function(k) {
    rv$P <- rdirichlet1(k)
    rv$Q <- rep(1/k, k)
    rv$sel <- NULL
    updateSliderInput(session, "adjustVal", value = 0)
  }
  observeEvent(input$k, { init_reset(input$k) }, priority = 10)
  observeEvent(input$newP, { init_reset(input$k) })

  # Select a guessed bar (only trace #2: red bars)
  observeEvent(event_data("plotly_click", source = "PQ"), {
    d <- event_data("plotly_click", source = "PQ")
    if (is.null(d)) return()
    if (!is.null(d$curveNumber) && d$curveNumber == 1) {
      idx <- d$pointNumber + 1L
      if (idx >= 1 && idx <= input$k) {
        rv$sel <- idx
        updateSliderInput(session, "adjustVal", value = rv$Q[idx])
      }
    } else {
      showNotification("Select a RED (Q) bar to edit.", type = "message", duration = 1.5)
    }
  })

  # Adjust via slider or buttons
  observeEvent(input$adjustVal, {
    if (is.null(rv$sel)) return()
    rv$Q <- adjust_q(rv$Q, rv$sel, input$adjustVal, eps = eps)
    updateSliderInput(session, "adjustVal", value = rv$Q[rv$sel])
  })
  observeEvent(input$incBtn, {
    if (is.null(rv$sel)) return()
    rv$Q <- adjust_q(rv$Q, rv$sel, rv$Q[rv$sel] + 0.01, eps = eps)
    updateSliderInput(session, "adjustVal", value = rv$Q[rv$sel])
  })
  observeEvent(input$decBtn, {
    if (is.null(rv$sel)) return()
    rv$Q <- adjust_q(rv$Q, rv$sel, rv$Q[rv$sel] - 0.01, eps = eps)
    updateSliderInput(session, "adjustVal", value = rv$Q[rv$sel])
  })

  # Metrics (JSD-focused)
  metrics <- reactive({
    req(rv$P, rv$Q)
    base <- input$base
    list(
      HP  = entropy(rv$P, base = base, eps = eps),
      HQ  = entropy(rv$Q, base = base, eps = eps),
      JSD = js_divergence(rv$P, rv$Q, base = base, eps = eps)
    )
  })

  output$metricsUI <- renderUI({
    m <- metrics()
    jsd <- m$JSD
    tagList(
      div(class = "metric-value",
          sprintf("H(P) = %.4f,   H(Q) = %.4f,   JSD(P || Q) = %.4f,   √JSD = %.4f",
                  m$HP, m$HQ, jsd, sqrt(jsd)))
    )
  })

  output$boundsText <- renderText({
    if (input$base == "bits") {
      "JSD is bounded in [0, 1] (bits). √JSD is the Jensen–Shannon distance (a metric in [0, 1])."
    } else {
      "JSD is bounded in [0, ln 2] (nats). √JSD is the Jensen–Shannon distance; its max is √(ln 2)."
    }
  })

  output$selectedInfo <- renderUI({
    k <- input$k
    if (is.null(rv$sel)) {
      div(class = "selected-note", "No selection. Click a RED bar (Q) to select a category.")
    } else {
      nm <- cat_names(k); idx <- rv$sel
      div(class = "selected-note",
          sprintf("Selected: %s (index %d). Current Q = %s",
                  nm[idx], idx, fmt_pct(rv$Q[idx])))
    }
  })

  # Combined plot (grouped bars)
  output$plotPQ <- renderPlotly({
    req(rv$P, rv$Q)
    k <- input$k
    nm <- cat_names(k)
    ymax <- max(0.3, 1.2 * max(rv$P, rv$Q))
    fontSmall <- list(size = 11)
    plot_ly(source = "PQ") |>
      add_bars(x = nm, y = rv$P, name = "TRUE P",
               marker = list(color = "rgba(31,119,180,0.85)"),
               hovertemplate = "%{x}<br>P = %{y:.4f}<extra></extra>") |>
      add_bars(x = nm, y = rv$Q, name = "GUESSED Q",
               marker = list(color = "rgba(214,39,40,0.85)"),
               hovertemplate = "%{x}<br>Q = %{y:.4f}<extra></extra>") |>
      layout(barmode = "group",
             title = list(text = "TRUE (blue) vs GUESSED (red)", font = list(size = 14)),
             legend = list(orientation = "h", y = -0.2, font = fontSmall),
             yaxis = list(range = c(0, ymax), title = "probability", tickfont = fontSmall, titlefont = fontSmall),
             xaxis = list(title = "", tickfont = fontSmall, titlefont = fontSmall),
             margin = list(l = 40, r = 10, t = 40, b = 40))
  })

  output$pqTable <- renderTable({
    req(rv$P, rv$Q)
    tibble::tibble(
      Category = cat_names(input$k),
      `P (TRUE)` = round(rv$P, 6),
      `Q (GUESSED)` = round(rv$Q, 6),
      `P - Q` = round(rv$P - rv$Q, 6)
    )
  })
}

if (interactive()) shinyApp(ui, server)
