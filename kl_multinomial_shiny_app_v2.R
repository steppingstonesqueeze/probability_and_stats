# KL Divergence Playground for Multinomial Distributions — v2
# - TRUE (P) and GUESSED (Q) shown on the SAME grouped bar plot
# - Click a RED bar to select; adjust via slider or +/- buttons (plotly bars aren't directly draggable)
# - Smaller plot fonts

suppressPackageStartupMessages({
  ok <- requireNamespace("shiny", quietly = TRUE); if (!ok) install.packages("shiny")
  ok <- requireNamespace("plotly", quietly = TRUE); if (!ok) install.packages("plotly")
  ok <- requireNamespace("dplyr", quietly = TRUE); if (!ok) install.packages("dplyr")
  ok <- requireNamespace("shinyWidgets", quietly = TRUE); if (!ok) install.packages("shinyWidgets")
})

library(shiny)
library(plotly)
library(dplyr)
library(shinyWidgets)

# ---- Helpers ----

rdirichlet1 <- function(k, alpha = NULL) {
  if (is.null(alpha)) alpha <- runif(k, 0.5, 2)
  x <- rgamma(k, shape = alpha, rate = 1)
  x / sum(x)
}

adjust_q <- function(q, i, new_i, eps = 1e-9) {
  k <- length(q)
  if (i < 1 || i > k) return(q)
  new_i <- max(min(new_i, 1 - eps), eps)
  old_i <- q[i]
  delta <- new_i - old_i
  others <- setdiff(seq_len(k), i)
  if (abs(delta) < 1e-15) return(q)

  if (delta > 0) {
    available <- sum(pmax(q[others] - eps, 0))
    delta_eff <- min(delta, available)
    if (delta_eff <= 0) return(q)
    decrement <- delta_eff / length(others)
    q[others] <- pmax(q[others] - decrement, eps)
    q[i] <- q[i] + delta_eff
  } else {
    inc_total <- -delta
    capacity <- sum(pmax((1 - 1e-9) - q[others], 0))
    inc_eff <- min(inc_total, capacity)
    if (inc_eff <= 0) return(q)
    increment <- inc_eff / length(others)
    q[others] <- pmin(q[others] + increment, 1 - 1e-9)
    q[i] <- q[i] - inc_eff
  }
  q <- pmax(q, eps)
  q / sum(q)
}

entropy <- function(p, base = "nats", eps = 1e-12) {
  p <- pmax(p, eps)
  if (base == "bits") return(-sum(p * log2(p)))
  -sum(p * log(p))
}
cross_entropy <- function(p, q, base = "nats", eps = 1e-12) {
  p <- pmax(p, eps); q <- pmax(q, eps)
  if (base == "bits") return(-sum(p * log2(q)))
  -sum(p * log(q))
}
kl_div <- function(p, q, base = "nats", eps = 1e-12) {
  p <- pmax(p, eps); q <- pmax(q, eps)
  if (base == "bits") return(sum(p * (log2(p) - log2(q))))
  sum(p * (log(p) - log(q)))
}

fmt_pct <- function(x) sprintf('%.2f%%', 100 * x)
cat_names <- function(k) paste0("C", seq_len(k))

# ---- UI ----
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .metric-card {padding: 10px 12px; border-radius: 10px; background: #f8f9fb; margin-bottom: 10px;}
      .metric-title {font-weight: 600; color: #444; font-size: 13px;}
      .metric-value {font-weight: 700; font-size: 18px; color: #111;}
      .subtle {color: #777; font-size: 12px;}
      .selected-note {font-size: 13px; color: #555;}
    "))
  ),
  titlePanel('KL Divergence Playground (Multinomial) — grouped bars'),
  sidebarLayout(
    sidebarPanel(
      sliderInput("k", "Number of categories (k)", min = 2, max = 10, value = 5, step = 1),
      radioButtons("base", "Log base", choices = c("nats", "bits"), selected = "nats", inline = TRUE),
      actionButton("newP", "NEW (sample a new TRUE distribution)", class = "btn btn-primary"),
      tags$hr(),
      uiOutput("selectedInfo"),
      sliderInput("adjustVal", "Adjust selected guessed bar", min = 0, max = 1, value = 0.0, step = 0.001),
      fluidRow(
        column(6, actionButton("decBtn", "– 0.01", class = "btn btn-outline-secondary btn-sm")),
        column(6, actionButton("incBtn", "+ 0.01", class = "btn btn-outline-secondary btn-sm"))
      ),
      helpText("Click a RED bar in the combined plot to select. Plotly doesn't support direct dragging of bar heights, so use the slider or +/- buttons to adjust. ")
    ),
    mainPanel(
      plotlyOutput("plotPQ", height = "370px"),
      br(),
      div(class = "metric-card",
          div(class = "metric-title", "Information Measures"),
          uiOutput("metricsUI"),
          div(class = "subtle", "Identity: H(P,Q) = H(P) + KL(P||Q)")
      ),
      div(class = "metric-card",
          div(class = "metric-title", "Details"),
          tableOutput("pqTable")
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  eps <- 1e-9
  rv <- reactiveValues(P=NULL, Q=NULL, sel=NULL)

  observeEvent(input$k, {
    k <- input$k
    rv$P <- rdirichlet1(k)
    rv$Q <- rep(1/k, k)
    rv$sel <- NULL
    updateSliderInput(session, "adjustVal", value = 0.0)
  }, priority = 10)

  observeEvent(input$newP, {
    k <- input$k
    rv$P <- rdirichlet1(k)
    rv$Q <- rep(1/k, k)
    rv$sel <- NULL
    updateSliderInput(session, "adjustVal", value = 0.0)
  })

  # Click selection on combined plot: only allow selecting Q (trace 2)
  observeEvent(event_data("plotly_click", source = "PQ"), {
    d <- event_data("plotly_click", source = "PQ")
    if (is.null(d)) return()
    # traceNumber is 0-based: 0 = P (blue), 1 = Q (red)
    if (!is.null(d$curveNumber) && d$curveNumber == 1) {
      idx <- d$pointNumber + 1L
      if (idx >= 1 && idx <= input$k) {
        rv$sel <- idx
        updateSliderInput(session, "adjustVal", value = rv$Q[idx])
      }
    } else {
      showNotification("Select the RED (GUESSED Q) bar to edit.", type = "message", duration = 2)
    }
  })

  # Slider adjust
  observeEvent(input$adjustVal, {
    if (is.null(rv$sel)) return()
    q_new <- adjust_q(rv$Q, rv$sel, input$adjustVal, eps = eps)
    rv$Q <- q_new
    updateSliderInput(session, "adjustVal", value = rv$Q[rv$sel])
  })

  # Nudge buttons +/- 0.01
  observeEvent(input$incBtn, {
    if (is.null(rv$sel)) return()
    target <- rv$Q[rv$sel] + 0.01
    rv$Q <- adjust_q(rv$Q, rv$sel, target, eps = eps)
    updateSliderInput(session, "adjustVal", value = rv$Q[rv$sel])
  })
  observeEvent(input$decBtn, {
    if (is.null(rv$sel)) return()
    target <- rv$Q[rv$sel] - 0.01
    rv$Q <- adjust_q(rv$Q, rv$sel, target, eps = eps)
    updateSliderInput(session, "adjustVal", value = rv$Q[rv$sel])
  })

  metrics <- reactive({
    req(rv$P, rv$Q)
    base <- input$base
    list(
      H  = entropy(rv$P, base = base, eps = eps),
      CE = cross_entropy(rv$P, rv$Q, base = base, eps = eps),
      KL = kl_div(rv$P, rv$Q, base = base, eps = eps)
    )
  })

  output$metricsUI <- renderUI({
    m <- metrics()
    div(class = "metric-value",
        sprintf("H(P) = %.4f,   H(P,Q) = %.4f,   KL(P || Q) = %.4f", m$H, m$CE, m$KL))
  })

  output$selectedInfo <- renderUI({
    k <- input$k
    if (is.null(rv$sel)) {
      div(class = "selected-note", "No selection. Click a RED bar (Q) to select a category.")
    } else {
      nm <- cat_names(k)
      idx <- rv$sel
      div(class = "selected-note",
          sprintf("Selected: %s (index %d). Current Q = %s", nm[idx], idx, fmt_pct(rv$Q[idx])))
    }
  })

  # Combined grouped bars
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
