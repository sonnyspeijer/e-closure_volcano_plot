# Libraries
suppressPackageStartupMessages(library(eClosure))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(shiny))

# Example data
load("data/example.RData")

# ---- Helper functions for closed e-BH ----
harmonic_number <- function(m) sum(1 / seq_len(m))

calibrator_BY <- function(p, alpha) {
  k <- length(p)
  h_k <- harmonic_number(k)
  
  indicator <- (h_k * p <= alpha)
  ceil <- pmax(ceiling((k * h_k * p) / alpha), 1)
  
  ifelse(indicator, k / (alpha * ceil), 0)
}

# ---- Frontend ----
ui <- fluidPage(
  titlePanel("e-Closure Volcano Plot"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload your CSV file here", accept = c(".csv")),
      helpText("Ensure that the CSV file contains the columns 'log2FoldChange' and 'pvalue'."),
      radioButtons("mode", "Filter mode:",
                   choices = c("Thresholds" = "threshold",
                               "Top values" = "top"),
                   selected = "threshold"),
      conditionalPanel(
        condition = "input.mode == 'threshold'",
        uiOutput("threshold_lf"),
        uiOutput("threshold_p")
      ),
      conditionalPanel(
        condition = "input.mode == 'top'",
        uiOutput("top_lf"),
        uiOutput("top_p")
      ),
      hr(),
      h4("Summary of selected features"),
      textOutput("n"),
      hr(),
      sliderInput("alpha", "Alpha:", min = 0, max = 1, value = 0.05, step = 0.01),
      uiOutput("results"),
      hr(),
      downloadButton("download", "Download selected features")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Output",
          plotOutput("volcano", height = "600px", width = "800px")
        ),
        tabPanel(
          "Help",
          h4("About"),
          p(HTML("This application implements an interactive volcano plot with <em>simultaneous</em> false discovery rate
                 (FDR) control. Simultaneous control means that the user can adjust the filters and the significance level
                 alpha after seeing the data without compromising FDR control.<br><br>This guarantee is achieved by using
                 methods that are constructed using the e-closure principle. Specifically, this application includes the
                 e-closure variants of three FDR control methods: the e-Benjamini Hochberg (e-BH) procedure, the
                 Benjamini-Yekutieli (BY) procedure, and the Su procedure based on the FDR-linking theorem.")),
          hr(),
          h4("How to use"),
          p(HTML("Features can be selected according to two filter modes. The <em>thresholds</em> mode selects features
                 if they exceed the threshold value for both filters. The <em>top values</em> mode selects features if
                 their absolute fold changes rank among the largest <em>m</em> and the p-values rank among the smallest
                 <em>n</em>.</br></br> Based on the filters, the summary panel returns how many features were selected and if
                 the selected features control the FDR at significance level alpha according to each e-closure method."))
        )
      )
    )
  )
)

# ---- Backend ----
server <- function(input, output, session) {
  # ---- Upload file ----
  react_df <- reactiveVal(df)
  required_cols <- c("log2FoldChange", "pvalue")
  
  observeEvent(input$file, {
    input_df <- read.csv(input$file$datapath, stringsAsFactors = FALSE)
    
    if (!all(required_cols %in% colnames(input_df))) {
      showNotification(
        "Uploaded file must contain columns 'log2FoldChange' and 'pvalue'.",
        type = "error", duration = 10
      )
      return(NULL)
    }
    
    react_df(input_df)
  })
  
  # ---- Interactive filter selection ----
  max_lf <- reactive(ceiling(max(abs(react_df()$log2FoldChange), na.rm = TRUE)))
  max_p  <- reactive(ceiling(max(-log10(react_df()$pvalue), na.rm = TRUE)))
  
  output$threshold_lf <- renderUI({
    sliderInput("thr_lf", "Threshold for absolute fold change (log2 scale):",
                min = 0, max = max_lf(), value = max_lf() / 10, step = 0.05)
  })
  output$threshold_p <- renderUI({
    sliderInput("thr_p", "Threshold for p-value (-log10 scale):",
                min = 0, max = max_p(), value = max_p() / 10, step = 0.5)
  })
  output$top_lf <- renderUI({
    numericInput("top_lf", "Top m largest absolute fold changes:",
                 min = 0, max = nrow(react_df()), value = round(nrow(react_df()) / 10), step = 1)
  })
  output$top_p <- renderUI({
    numericInput("top_p", "Top n smallest p-values:",
                 min = 0, max = nrow(react_df()), value = round(nrow(react_df()) / 10), step = 1)
  })
  
  # ---- Mark selected features ----
  selected <- reactive({
    r_df <- react_df()
    
    if (input$mode == "threshold") {
      req(input$thr_lf, input$thr_p)
      abs(r_df$log2FoldChange) >= input$thr_lf & -log10(r_df$pvalue) >= input$thr_p
    } else {
      req(input$top_lf, input$top_p)
      rank_lf <- rank(-abs(r_df$log2FoldChange), ties.method = "min")
      rank_p  <- rank(r_df$pvalue, ties.method = "min")
      rank_lf <= input$top_lf & rank_p <= input$top_p
    }
  })
  
  # ---- Volcano plot ----
  output$volcano <- renderPlot({
    r_df <- react_df()
    r_df$Selected <- ifelse(selected(), "Selected", "Not selected")
    
    p <- ggplot(r_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Selected)) +
      geom_point(size = 1) +
      scale_color_manual(values = c("Not selected" = "#969696",
                                    "Selected" = "#785EF0"),
                         name = "Status") +
      labs(x = expression(log[2] * "(fold change)"),
           y = expression(-log[10] * "(p-value)")) +
      theme_test(base_size = 14) +
      theme(legend.position = "right")
    
    if (input$mode == "threshold") {
      p <- p +
        geom_vline(xintercept = c(-input$thr_lf, input$thr_lf), linetype = "dashed", color = "black") +
        geom_hline(yintercept = input$thr_p, linetype = "dashed", color = "black")
    } else {
      lf_cutoff <- sort(abs(r_df$log2FoldChange), decreasing = TRUE)[input$top_lf]
      p_cutoff <- sort(r_df$pvalue)[input$top_p]
      
      p <- p +
        geom_vline(xintercept = c(-lf_cutoff, lf_cutoff), linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black")
    }
    
    p
  })
  
  # ---- Summary panel ----
  output$n <- renderText({
    paste("Features selected:", sum(selected()))
  })
  
  output$results <- renderUI({
    r_df  <- react_df()
    alpha <- input$alpha
    set   <- selected()
    
    # Note: the calibrator uses its own alpha, not the global alpha!
    e <- calibrator_BY(r_df$pvalue, alpha = 0.05)
    
    r_eBH <- closedeBH(e = e, set = set, alpha = alpha)
    r_BY  <- closedBY(p = r_df$pvalue, set = set, alpha = alpha)
    r_Su  <- closedSu(p = r_df$pvalue, set = set, alpha = alpha)
    
    validity <- function(method, valid) {
      color <- if (isTRUE(valid)) "#009465" else "#CB2F43"
      text  <- if (isTRUE(valid)) {
        paste("Valid under closed", method)
      } else {
        paste("Not valid under closed", method)
      }
      div(
        style = "display: flex; align-items: center; gap: 8px; margin-bottom: 4px;",
        div(style = paste0(
          "width: 16px; height: 16px; border-radius: 4px; background-color:", color, ";"
        )),
        span(text)
      )
    }
    
    tagList(
      validity("e-BH", r_eBH),
      validity("BY", r_BY),
      validity("Su", r_Su)
    )
  })
  
  output$download <- downloadHandler(
    filename = function() paste0("selected_features_", Sys.Date(), ".csv"),
    content = function(file) {
      r_df <- react_df()
      out  <- r_df[selected(), , drop = FALSE]
      write.csv(out, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
