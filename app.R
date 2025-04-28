# ----------------------------------------------------
# Full Shiny App: Dose–Response Meta-Analysis with Advanced Features
# (Meta-regression, Custom Spline, and Prediction Intervals)
# ----------------------------------------------------

library(shiny)
library(bs4Dash)
library(dosresmeta)
library(ggplot2)
library(DT)
library(rms)       # for rcs()
library(shinyjs)   # for download functionality
library(metafor)   # for meta-regression and funnel/Egger tests
library(dplyr)

# ---------------------------
# Helper Functions
# ---------------------------
filter_data <- function(data) {
  # Remove studies with fewer than two non-reference doses.
  dat_list <- split(data, data$id)
  dat_list <- lapply(dat_list, function(x) {
    ref <- min(x$dose)
    if(sum(x$dose != ref) < 2) return(NULL) else return(x)
  })
  do.call(rbind, dat_list)
}

# Standardize predicted objects to ensure matching lengths
fix_pred <- function(pred) {
  if (is.null(pred)) return(NULL)
  n <- min(length(pred$dose), length(pred$pred))
  pred$dose <- pred$dose[1:n]
  pred$pred <- pred$pred[1:n]
  if(!is.null(pred$ci.lb)) pred$ci.lb <- pred$ci.lb[1:n]
  if(!is.null(pred$ci.ub)) pred$ci.ub <- pred$ci.ub[1:n]
  if(!is.null(pred$pi.lb)) pred$pi.lb <- pred$pi.lb[1:n]
  if(!is.null(pred$pi.ub)) pred$pi.ub <- pred$pi.ub[1:n]
  pred
}

clean_pred <- function(pred) {
  if (is.null(pred)) return(NULL)
  pred <- fix_pred(pred)
  good <- complete.cases(pred$dose, pred$pred)
  if (!is.null(pred$ci.lb)) good <- good & complete.cases(pred$ci.lb)
  if (!is.null(pred$ci.ub)) good <- good & complete.cases(pred$ci.ub)
  if (!is.null(pred$pi.lb)) good <- good & complete.cases(pred$pi.lb)
  if (!is.null(pred$pi.ub)) good <- good & complete.cases(pred$pi.ub)
  pred$dose <- pred$dose[good]
  pred$pred <- pred$pred[good]
  if (!is.null(pred$ci.lb)) pred$ci.lb <- pred$ci.lb[good]
  if (!is.null(pred$ci.ub)) pred$ci.ub <- pred$ci.ub[good]
  if (!is.null(pred$pi.lb)) pred$pi.lb <- pred$pi.lb[good]
  if (!is.null(pred$pi.ub)) pred$pi.ub <- pred$pi.ub[good]
  if(length(pred$dose) < 2) return(NULL)
  pred
}

make_forest_data <- function(dat) {
  dat$RR    <- exp(dat$logrr)
  dat$RR_lb <- exp(dat$logrr - 1.96 * dat$se)
  dat$RR_ub <- exp(dat$logrr + 1.96 * dat$se)
  dat
}

make_data_summary <- function(dat) {
  n_obs <- nrow(dat)
  n_studies <- length(unique(dat$id))
  dose_range <- range(dat$dose, na.rm = TRUE)
  avg_dose <- round(mean(dat$dose, na.rm = TRUE), 2)
  paste("Number of observations:", n_obs, "\n",
        "Number of studies:", n_studies, "\n",
        "Dose range:", paste(dose_range, collapse = " - "), "\n",
        "Average dose:", avg_dose)
}

# ---------------------------
# UI
# ---------------------------
ui <- bs4DashPage(
  header = bs4DashNavbar(title = "Dose–Response Meta-Analysis App"),
  sidebar = bs4DashSidebar(
    skin = "light",
    bs4SidebarMenu(
      bs4SidebarMenuItem("Data", tabName = "data", icon = icon("table")),
      bs4SidebarMenuItem("Data Summary", tabName = "datasum", icon = icon("info")),
      bs4SidebarMenuItem("Linear Model", tabName = "linear", icon = icon("chart-line")),
      bs4SidebarMenuItem("Quadratic Model", tabName = "quadratic", icon = icon("project-diagram")),
      bs4SidebarMenuItem("Spline Model", tabName = "spline", icon = icon("wave-square")),
      bs4SidebarMenuItem("Residuals", tabName = "residual", icon = icon("chart-area")),
      bs4SidebarMenuItem("Fitted vs. Observed", tabName = "fitObs", icon = icon("line-chart")),
      bs4SidebarMenuItem("Q-Q Plot", tabName = "qq", icon = icon("sliders-h")),
      bs4SidebarMenuItem("Residual Histogram", tabName = "hist", icon = icon("bar-chart")),
      bs4SidebarMenuItem("Bubble Plot", tabName = "bubble", icon = icon("circle")),
      bs4SidebarMenuItem("Study Trends", tabName = "studytrend", icon = icon("th")),
      bs4SidebarMenuItem("Forest Plot", tabName = "forest", icon = icon("list")),
      bs4SidebarMenuItem("Model Comparison", tabName = "compare", icon = icon("balance-scale")),
      bs4SidebarMenuItem("Model Coefficients", tabName = "coef", icon = icon("table")),
      bs4SidebarMenuItem("CSV Template", tabName = "csv", icon = icon("file")),
      bs4SidebarMenuItem("Options", tabName = "options", icon = icon("sliders-h")),
      bs4SidebarMenuItem("Download Plots", tabName = "download", icon = icon("download")),
      bs4SidebarMenuItem("About", tabName = "about", icon = icon("info-circle")),
      bs4SidebarMenuItem("Advanced Analysis", tabName = "advanced", icon = icon("cogs"))
    )
  ),
  body = bs4DashBody(
    useShinyjs(),
    bs4TabItems(
      # --- Data Tab ---
      bs4TabItem(
        tabName = "data",
        fluidRow(
          box(title = "Upload Data (CSV)", width = 6, status = "primary", solidHeader = TRUE,
              fileInput("file", "Choose CSV File", 
                        accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
              helpText("Required columns: id, dose, logrr, se, type, cases, n. If none is uploaded, the built-in alcohol_cvd dataset is used.")
          ),
          box(title = "Data Preview", width = 6, status = "info", solidHeader = TRUE,
              DTOutput("dataTable"))
        )
      ),
      # --- Data Summary Tab ---
      bs4TabItem(
        tabName = "datasum",
        fluidRow(
          box(title = "Data Summary", width = 12, status = "info", solidHeader = TRUE,
              verbatimTextOutput("dataSummary"))
        )
      ),
      # --- Linear Model Tab ---
      bs4TabItem(
        tabName = "linear",
        fluidRow(
          box(title = "Linear Model Summary", width = 12, status = "primary", solidHeader = TRUE,
              verbatimTextOutput("linSummary"))
        ),
        fluidRow(
          box(title = "Linear Model Prediction Plot", width = 12, status = "primary", solidHeader = TRUE,
              plotOutput("linPlot", height = "400px"),
              downloadButton("downloadLinPlot", "Download Linear Plot")
          )
        )
      ),
      # --- Quadratic Model Tab ---
      bs4TabItem(
        tabName = "quadratic",
        fluidRow(
          box(title = "Quadratic Model Summary", width = 12, status = "primary", solidHeader = TRUE,
              verbatimTextOutput("quadSummary"))
        ),
        fluidRow(
          box(title = "Quadratic Model Prediction Plot", width = 12, status = "primary", solidHeader = TRUE,
              plotOutput("quadPlot", height = "400px"),
              downloadButton("downloadQuadPlot", "Download Quadratic Plot")
          )
        )
      ),
      # --- Spline Model Tab ---
      bs4TabItem(
        tabName = "spline",
        fluidRow(
          box(title = "Spline Model Summary", width = 12, status = "warning", solidHeader = TRUE,
              verbatimTextOutput("splineSummary"))
        ),
        fluidRow(
          box(title = "Spline Model Prediction Plot", width = 12, status = "warning", solidHeader = TRUE,
              plotOutput("splinePlot", height = "400px"),
              downloadButton("downloadSplinePlot", "Download Spline Plot")
          )
        )
      ),
      # --- Residuals Tab ---
      bs4TabItem(
        tabName = "residual",
        fluidRow(
          box(title = "Residual Plot (Linear Model)", width = 12, status = "danger", solidHeader = TRUE,
              plotOutput("residPlot", height = "400px"),
              downloadButton("downloadResidPlot", "Download Residual Plot")
          )
        )
      ),
      # --- Fitted vs. Observed Tab ---
      bs4TabItem(
        tabName = "fitObs",
        fluidRow(
          box(title = "Fitted vs. Observed (Linear Model)", width = 12, status = "info", solidHeader = TRUE,
              plotOutput("fitObsPlot", height = "400px"),
              downloadButton("downloadFitObsPlot", "Download Fitted vs. Observed Plot")
          )
        )
      ),
      # --- Q-Q Plot Tab ---
      bs4TabItem(
        tabName = "qq",
        fluidRow(
          box(title = "Q-Q Plot of Residuals (Linear Model)", width = 12, status = "info", solidHeader = TRUE,
              plotOutput("qqPlot", height = "400px"),
              downloadButton("downloadQQPlot", "Download Q-Q Plot")
          )
        )
      ),
      # --- Residual Histogram Tab ---
      bs4TabItem(
        tabName = "hist",
        fluidRow(
          box(title = "Histogram of Residuals (Linear Model)", width = 12, status = "info", solidHeader = TRUE,
              plotOutput("histPlot", height = "400px"),
              downloadButton("downloadHistPlot", "Download Residual Histogram")
          )
        )
      ),
      # --- Bubble Plot Tab ---
      bs4TabItem(
        tabName = "bubble",
        fluidRow(
          box(title = "Bubble Plot of Raw Data", width = 12, status = "info", solidHeader = TRUE,
              plotOutput("bubblePlot", height = "400px"),
              downloadButton("downloadBubblePlot", "Download Bubble Plot")
          )
        )
      ),
      # --- Study Trends Tab ---
      bs4TabItem(
        tabName = "studytrend",
        fluidRow(
          box(title = "Study-Specific Trends", width = 12, status = "primary", solidHeader = TRUE,
              plotOutput("studyTrendPlot", height = "400px"),
              downloadButton("downloadStudyTrendPlot", "Download Study Trends Plot")
          )
        )
      ),
      # --- Forest Plot Tab ---
      bs4TabItem(
        tabName = "forest",
        fluidRow(
          box(title = "Forest Plot (Raw Data)", width = 12, status = "primary", solidHeader = TRUE,
              plotOutput("forestPlot", height = "400px"),
              downloadButton("downloadForestPlot", "Download Forest Plot")
          )
        )
      ),
      # --- Model Comparison Tab ---
      bs4TabItem(
        tabName = "compare",
        fluidRow(
          box(title = "Model Fit Comparison", width = 12, status = "success", solidHeader = TRUE,
              DTOutput("modelTable"))
        )
      ),
      # --- Model Coefficients Tab ---
      bs4TabItem(
        tabName = "coef",
        fluidRow(
          box(title = "Model Coefficients", width = 12, status = "primary", solidHeader = TRUE,
              DTOutput("coefTable"))
        )
      ),
      # --- CSV Template Tab ---
      bs4TabItem(
        tabName = "csv",
        fluidRow(
          box(title = "CSV Template and Instructions", width = 12, status = "info", solidHeader = TRUE,
              HTML("<h3>CSV Template</h3>
                    <p>Your CSV file must include the following columns (with exactly these names):</p>
                    <ul>
                      <li><b>id</b>: Study identifier (e.g., 1, 2, 3...)</li>
                      <li><b>dose</b>: The dose value (numeric)</li>
                      <li><b>logrr</b>: The natural logarithm of the relative risk (numeric)</li>
                      <li><b>se</b>: Standard error for logrr (numeric)</li>
                      <li><b>type</b>: Study design type (e.g., 'ci')</li>
                      <li><b>cases</b>: Number of cases (numeric)</li>
                      <li><b>n</b>: Total number of subjects (numeric)</li>
                    </ul>
                    <p>A sample CSV content is shown below:</p>
                    <pre>
id,dose,logrr,se,type,cases,n
1,0,0,0.05,ci,20,200
1,10,-0.1,0.07,ci,25,200
1,20,-0.2,0.10,ci,28,200
2,0,0,0.06,ci,15,220
2,15,-0.05,0.08,ci,18,220
2,30,-0.1,0.10,ci,22,220
                    </pre>")
          )
        )
      ),
      # --- Options Tab ---
      bs4TabItem(
        tabName = "options",
        fluidRow(
          box(title = "Plot Options", width = 12, status = "info", solidHeader = TRUE,
              selectInput("method", "Estimation Method", choices = c("fixed", "reml"), selected = "fixed"),
              numericInput("xref", "Spline Reference Dose (xref)", value = 0, min = 0),
              sliderInput("nPoints", "Number of Prediction Points", min = 50, max = 200, value = 100),
              sliderInput("ciLevel", "Confidence Level", min = 0.80, max = 0.99, value = 0.95, step = 0.01),
              checkboxInput("showCI", "Show Confidence Interval Lines", value = TRUE),
              checkboxInput("showRaw", "Overlay Raw Data on Prediction Plots", value = FALSE),
              selectInput("themeChoice", "Plot Theme", choices = c("Minimal", "Black & White"), selected = "Minimal"),
              sliderInput("ptSize", "Raw Data Point Size", min = 1, max = 6, value = 3, step = 0.5),
              sliderInput("alpha", "Raw Data Transparency", min = 0.1, max = 1, value = 0.7),
              sliderInput("lineWidth", "Prediction Line Width", min = 1, max = 5, value = 2, step = 0.5),
              checkboxInput("includeSpline", "Include Spline in Combined Trends", value = TRUE),
              checkboxInput("showCombinedLegend", "Show Legend in Combined Trends", value = TRUE),
              sliderInput("textSize", "Plot Text Size (for ggplot)", min = 8, max = 20, value = 12),
              checkboxInput("showGrid", "Show Grid Lines (for ggplot)", value = TRUE),
              textInput("bgColor", "Plot Background Color", value = "white")
          )
        )
      ),
      # --- Download Plots Tab ---
      bs4TabItem(
        tabName = "download",
        fluidRow(
          box(title = "Download All Plots", width = 12, status = "info", solidHeader = TRUE,
              downloadButton("downloadLinPlot", "Download Linear Plot"),
              downloadButton("downloadQuadPlot", "Download Quadratic Plot"),
              downloadButton("downloadSplinePlot", "Download Spline Plot"),
              downloadButton("downloadResidPlot", "Download Residual Plot"),
              downloadButton("downloadFitObsPlot", "Download Fitted vs. Observed Plot"),
              downloadButton("downloadQQPlot", "Download Q-Q Plot"),
              downloadButton("downloadHistPlot", "Download Residual Histogram"),
              downloadButton("downloadBubblePlot", "Download Bubble Plot"),
              downloadButton("downloadStudyTrendPlot", "Download Study Trends Plot"),
              downloadButton("downloadForestPlot", "Download Forest Plot")
          )
        )
      ),
      # --- About Tab ---
      bs4TabItem(
        tabName = "about",
        fluidRow(
          box(title = "About This App", width = 12, status = "info", solidHeader = TRUE,
              HTML("<h3>Dose–Response Meta-Analysis App</h3>
                    <p>This app implements dose–response meta-analysis using the dosresmeta package.
                    It fits linear, quadratic, and restricted cubic spline models and includes advanced analyses such as meta‐regression and custom spline modeling.
                    Developed by [Your Name].</p>")
          )
        )
      ),
      # --- Advanced Analysis Tab ---
      bs4TabItem(
        tabName = "advanced",
        fluidRow(
          box(title = "Advanced Analysis", width = 12, status = "primary", solidHeader = TRUE,
              h4("Heterogeneity & Model Diagnostics"),
              verbatimTextOutput("hetStats"),
              
              h4("Subgroup Analysis"),
              helpText("Choose a column in your dataset for subgrouping (e.g., sex or region)."),
              uiOutput("subgroupUI"),
              actionButton("runSubgroup", "Run Subgroup Models"),
              verbatimTextOutput("subgroupResults"),
              
              h4("Quartile-Based Analysis"),
              helpText("Classify 'dose' into quartiles and meta-analyze each quartile vs. Q1 reference."),
              actionButton("runQuartile", "Run Quartile Analysis"),
              plotOutput("quartileForest", height = "400px"),
              
              h4("Publication Bias Checks (Funnel, Egger)"),
              helpText("Use each study's highest vs. lowest dose to assess publication bias."),
              actionButton("runPubBias", "Check Publication Bias"),
              plotOutput("funnelPlot", height = "400px"),
              verbatimTextOutput("eggerTest"),
              
              tags$hr(),
              h4("Leave-One-Out Sensitivity Analysis (Linear Model)"),
              helpText("Removes each study in turn, refits the linear model, and records the 'dose' coefficient."),
              actionButton("runLoo", "Run Leave-One-Out"),
              verbatimTextOutput("looResults"),
              div(style = "display:none;", textInput("dummyLOOTxt", "", value = "")),
              
              tags$hr(),
              h4("Meta-Regression Analysis"),
              helpText("Select one or more numeric moderator variables (e.g., mean_age, follow_up_years) to examine their influence on the dose–response relationship."),
              selectizeInput("metaRegVars", "Select Moderator(s) for Meta-Regression:", choices = NULL, multiple = TRUE),
              actionButton("runMetaReg", "Run Meta-Regression"),
              verbatimTextOutput("metaRegSummary"),
              
              tags$hr(),
              h4("Custom Spline Model Analysis"),
              helpText("Specify the number of knots and their percentiles (comma separated, e.g., 10,50,90) for a custom spline model."),
              numericInput("customNumKnots", "Number of Knots:", value = 3, min = 2, max = 5),
              textInput("customKnotsPercentiles", "Knots Percentiles (comma separated):", "10,50,90"),
              actionButton("runCustomSpline", "Run Custom Spline Model"),
              verbatimTextOutput("customSplineSummary"),
              plotOutput("customSplinePlot", height = "400px")
          )
        )
      )
    )
  ),
  controlbar = bs4DashControlbar(),
  footer = bs4DashFooter()
)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session) {
  
  # Reactive Data: Load uploaded CSV or built-in dataset.
  dataInput <- reactive({
    if (is.null(input$file)) {
      data("alcohol_cvd", package = "dosresmeta")
      dat <- alcohol_cvd
    } else {
      req(input$file)
      dat <- read.csv(input$file$datapath, stringsAsFactors = FALSE)
      dat$id    <- as.factor(dat$id)
      dat$dose  <- as.numeric(dat$dose)
      dat$logrr <- as.numeric(dat$logrr)
      dat$se    <- as.numeric(dat$se)
      dat$type  <- as.character(dat$type)
      dat$cases <- as.numeric(dat$cases)
      dat$n     <- as.numeric(dat$n)
    }
    dat <- dat[order(dat$id, dat$dose), ]
    dat <- filter_data(dat)
    if(is.null(dat) || nrow(dat) == 0)
      stop("No studies with at least two non-reference dose groups remain after filtering.")
    dat
  })
  
  output$dataTable <- renderDT({ dataInput() })
  output$dataSummary <- renderPrint({ make_data_summary(dataInput()) })
  
  ## Reactive Models
  linearModel <- reactive({
    dat <- dataInput()
    dosresmeta(
      formula = logrr ~ dose,
      id = id, type = type, se = se, cases = cases, n = n,
      data = dat, method = input$method
    )
  })
  output$linSummary <- renderPrint({ summary(linearModel()) })
  
  output$linPlot <- renderPlot({
    mod <- linearModel(); dat <- dataInput()
    dose_seq <- seq(min(dat$dose), max(dat$dose), length.out = input$nPoints)
    newdata <- data.frame(dose = dose_seq)
    pred <- tryCatch({
      predict(mod, newdata, order = TRUE, exp = TRUE, ci.level = input$ciLevel, pi = TRUE)
    }, error = function(e) { showNotification("Linear model prediction failed.", type = "error"); return(NULL) })
    pred <- clean_pred(pred)
    if(is.null(pred)) return(NULL)
    
    plot(pred$dose, pred$pred, type = "l", col = "blue", lwd = input$lineWidth,
         cex.axis = input$textSize/12, cex.lab = input$textSize/12, cex.main = input$textSize/12,
         ylim = range(c(pred$ci.lb, pred$ci.ub, pred$pi.lb, pred$pi.ub), na.rm = TRUE),
         xlab = "Dose (grams/day)", ylab = "Relative Risk",
         main = "Linear Model Predictions", bg = input$bgColor)
    if(input$showCI && !is.null(pred$ci.lb)){
      lines(pred$dose, pred$ci.lb, lty = 2, col = "blue", lwd = input$lineWidth)
      lines(pred$dose, pred$ci.ub, lty = 2, col = "blue", lwd = input$lineWidth)
    }
    if(!is.null(pred$pi.lb)){
      lines(pred$dose, pred$pi.lb, lty = 3, col = "darkblue", lwd = input$lineWidth)
      lines(pred$dose, pred$pi.ub, lty = 3, col = "darkblue", lwd = input$lineWidth)
    }
    if(input$showRaw){
      points(dat$dose, exp(dat$logrr), pch = 16,
             col = if(input$themeChoice=="Black & White") "black" else "gray",
             cex = input$ptSize, alpha = input$alpha)
    }
    if(input$showGrid) grid()
  })
  
  quadraticModel <- reactive({
    dat <- dataInput()
    dosresmeta(
      formula = logrr ~ dose + I(dose^2),
      id = id, type = type, se = se, cases = cases, n = n,
      data = dat, method = input$method
    )
  })
  output$quadSummary <- renderPrint({ summary(quadraticModel()) })
  
  output$quadPlot <- renderPlot({
    mod <- quadraticModel(); dat <- dataInput()
    dose_seq <- seq(min(dat$dose), max(dat$dose), length.out = input$nPoints)
    newdata <- data.frame(dose = dose_seq)
    pred <- tryCatch({
      predict(mod, newdata, order = TRUE, exp = TRUE, ci.level = input$ciLevel, pi = TRUE)
    }, error = function(e) { showNotification("Quadratic model prediction failed.", type = "error"); return(NULL) })
    pred <- clean_pred(pred)
    if(is.null(pred)) return(NULL)
    
    plot(pred$dose, pred$pred, type = "l", col = "darkgreen", lwd = input$lineWidth,
         cex.axis = input$textSize/12, cex.lab = input$textSize/12, cex.main = input$textSize/12,
         ylim = range(c(pred$ci.lb, pred$ci.ub, pred$pi.lb, pred$pi.ub), na.rm = TRUE),
         xlab = "Dose (grams/day)", ylab = "Relative Risk",
         main = "Quadratic Model Predictions", bg = input$bgColor)
    if(input$showCI && !is.null(pred$ci.lb)){
      lines(pred$dose, pred$ci.lb, lty = 2, col = "darkgreen", lwd = input$lineWidth)
      lines(pred$dose, pred$ci.ub, lty = 2, col = "darkgreen", lwd = input$lineWidth)
    }
    if(!is.null(pred$pi.lb)){
      lines(pred$dose, pred$pi.lb, lty = 3, col = "darkgreen", lwd = input$lineWidth)
      lines(pred$dose, pred$pi.ub, lty = 3, col = "darkgreen", lwd = input$lineWidth)
    }
    if(input$showRaw){
      points(dat$dose, exp(dat$logrr), pch = 16,
             col = if(input$themeChoice=="Black & White") "black" else "gray",
             cex = input$ptSize, alpha = input$alpha)
    }
    if(input$showGrid) grid()
  })
  
  splineModel <- reactive({
    dat <- dataInput()
    knots <- quantile(dat$dose, probs = c(0.05, 0.35, 0.65, 0.95), na.rm = TRUE)
    dosresmeta(
      formula = logrr ~ rcs(dose, knots),
      id = id, type = type, se = se, cases = cases, n = n,
      data = dat, method = input$method
    )
  })
  output$splineSummary <- renderPrint({ summary(splineModel()) })
  
  output$splinePlot <- renderPlot({
    mod <- splineModel(); dat <- dataInput()
    dose_seq <- seq(min(dat$dose), max(dat$dose), length.out = input$nPoints)
    newdata <- data.frame(dose = dose_seq)
    pred <- tryCatch({
      predict(mod, newdata, xref = input$xref, order = TRUE, exp = TRUE, ci.level = input$ciLevel, pi = TRUE)
    }, error = function(e) { showNotification("Spline model prediction failed.", type = "error"); return(NULL) })
    pred <- clean_pred(pred)
    if(is.null(pred)) return(NULL)
    
    plot(pred$dose, pred$pred, type = "l", col = "purple", lwd = input$lineWidth,
         cex.axis = input$textSize/12, cex.lab = input$textSize/12, cex.main = input$textSize/12,
         ylim = range(c(pred$ci.lb, pred$ci.ub, pred$pi.lb, pred$pi.ub), na.rm = TRUE),
         xlab = "Dose (grams/day)", ylab = "Relative Risk",
         main = "Spline Model Predictions", bg = input$bgColor)
    if(input$showCI && !is.null(pred$ci.lb)){
      lines(pred$dose, pred$ci.lb, lty = 2, col = "purple", lwd = input$lineWidth)
      lines(pred$dose, pred$ci.ub, lty = 2, col = "purple", lwd = input$lineWidth)
    }
    if(!is.null(pred$pi.lb)){
      lines(pred$dose, pred$pi.lb, lty = 3, col = "purple", lwd = input$lineWidth)
      lines(pred$dose, pred$pi.ub, lty = 3, col = "purple", lwd = input$lineWidth)
    }
    if(input$showRaw){
      points(dat$dose, exp(dat$logrr), pch = 16,
             col = if(input$themeChoice=="Black & White") "black" else "gray",
             cex = input$ptSize, alpha = input$alpha)
    }
    if(input$showGrid) grid()
  })
  
  output$residPlot <- renderPlot({
    mod <- linearModel(); dat <- dataInput()
    res <- residuals(mod); fitted_vals <- mod$fitted.values
    n_obs <- min(length(res), length(fitted_vals), nrow(dat))
    res <- res[1:n_obs]; xvals <- dat$dose[1:n_obs]
    plot(xvals, res, xlab = "Dose (grams/day)", ylab = "Residuals",
         main = "Residuals vs Dose (Linear Model)", pch = 16, col = "red",
         cex.axis = input$textSize/12, cex.lab = input$textSize/12, cex.main = input$textSize/12)
    abline(h = 0, lty = 2)
    if(input$showGrid) grid()
  })
  
  output$fitObsPlot <- renderPlot({
    mod <- linearModel(); dat <- dataInput()
    fitted_vals <- mod$fitted.values; observed <- exp(dat$logrr)
    n_obs <- min(length(fitted_vals), length(observed))
    df <- data.frame(Fitted = fitted_vals[1:n_obs], Observed = observed[1:n_obs])
    p <- ggplot(df, aes(x = Fitted, y = Observed)) +
      geom_point(size = 3, color = "blue") +
      geom_smooth(method = "loess", se = FALSE, color = "black") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
      labs(x = "Fitted Relative Risk", y = "Observed Relative Risk",
           title = "Fitted vs. Observed (Linear Model)") +
      theme_minimal(base_size = input$textSize) +
      theme(panel.grid.major = if(input$showGrid) element_line() else element_blank(),
            panel.grid.minor = if(input$showGrid) element_line() else element_blank(),
            plot.background = element_rect(fill = input$bgColor))
    print(p)
  })
  
  output$qqPlot <- renderPlot({
    mod <- linearModel()
    res <- residuals(mod)
    qqnorm(res, main = "Q-Q Plot of Residuals (Linear Model)",
           cex.axis = input$textSize/12, cex.lab = input$textSize/12, cex.main = input$textSize/12)
    qqline(res, col = "red", lty = 2)
  })
  
  output$histPlot <- renderPlot({
    mod <- linearModel()
    res <- residuals(mod)
    hist(res, breaks = 20, main = "Histogram of Residuals (Linear Model)",
         xlab = "Residuals", col = "lightblue", border = "white",
         cex.axis = input$textSize/12, cex.lab = input$textSize/12, cex.main = input$textSize/12)
  })
  
  output$bubblePlot <- renderPlot({
    dat <- dataInput()
    ggplot(dat, aes(x = dose, y = logrr, size = 1/se, color = id)) +
      geom_point(alpha = input$alpha) +
      scale_size_area(max_size = 10) +
      labs(x = "Dose (grams/day)", y = "Log Relative Risk", size = "1/SE", color = "Study ID") +
      theme_minimal(base_size = input$textSize) +
      theme(panel.grid.major = if(input$showGrid) element_line() else element_blank(),
            panel.grid.minor = if(input$showGrid) element_line() else element_blank(),
            plot.background = element_rect(fill = input$bgColor)) +
      ggtitle("Bubble Plot of Raw Data")
  })
  
  output$studyTrendPlot <- renderPlot({
    dat <- dataInput()
    ggplot(dat, aes(x = dose, y = logrr)) +
      geom_point(size = input$ptSize, color = if(input$themeChoice=="Black & White") "black" else "darkblue") +
      geom_line(color = if(input$themeChoice=="Black & White") "black" else "darkblue") +
      facet_wrap(~ id) +
      labs(x = "Dose (grams/day)", y = "Log Relative Risk") +
      theme_minimal(base_size = input$textSize) +
      theme(panel.grid.major = if(input$showGrid) element_line() else element_blank(),
            panel.grid.minor = if(input$showGrid) element_line() else element_blank(),
            plot.background = element_rect(fill = input$bgColor)) +
      ggtitle("Study-Specific Trends")
  })
  
  output$forestPlot <- renderPlot({
    dat <- dataInput()
    fdat <- make_forest_data(dat)
    ggplot(fdat, aes(x = RR, y = reorder(as.factor(id), as.numeric(id)))) +
      geom_point(size = 3, color = "darkred") +
      geom_errorbarh(aes(xmin = RR_lb, xmax = RR_ub), height = 0.2, color = "darkred") +
      labs(x = "Relative Risk", y = "Study ID", title = "Forest Plot (Raw Data)") +
      theme_minimal(base_size = input$textSize) +
      theme(panel.grid.major = if(input$showGrid) element_line() else element_blank(),
            panel.grid.minor = if(input$showGrid) element_line() else element_blank(),
            plot.background = element_rect(fill = input$bgColor))
  })
  
  output$modelTable <- renderDT({
    models <- list(
      Linear = tryCatch(linearModel(), error = function(e) NULL),
      Quadratic = tryCatch(quadraticModel(), error = function(e) NULL),
      Spline = tryCatch(splineModel(), error = function(e) NULL)
    )
    comp <- lapply(names(models), function(name) {
      mod <- models[[name]]
      if(!is.null(mod)) {
        aic_val <- tryCatch(AIC(mod), error = function(e) NA)
        bic_val <- tryCatch(BIC(mod), error = function(e) NA)
      } else {
        aic_val <- NA; bic_val <- NA
      }
      data.frame(Model = name, AIC = aic_val, BIC = bic_val)
    })
    do.call(rbind, comp)
  })
  
  output$coefTable <- renderDT({
    lin_df <- tryCatch({
      coefs <- coef(linearModel())
      ses   <- sqrt(diag(vcov(linearModel())))
      data.frame(Model = "Linear", Term = names(coefs), Estimate = round(coefs, 3),
                 SE = round(ses, 3), stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    quad_df <- tryCatch({
      coefs <- coef(quadraticModel())
      ses   <- sqrt(diag(vcov(quadraticModel())))
      data.frame(Model = "Quadratic", Term = names(coefs), Estimate = round(coefs, 3),
                 SE = round(ses, 3), stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    spl_df <- tryCatch({
      coefs <- coef(splineModel())
      ses   <- sqrt(diag(vcov(splineModel())))
      data.frame(Model = "Spline", Term = names(coefs), Estimate = round(coefs, 3),
                 SE = round(ses, 3), stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    do.call(rbind, list(lin_df, quad_df, spl_df))
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste("filtered_data-", Sys.Date(), ".csv", sep = "") },
    content = function(file) { write.csv(dataInput(), file, row.names = FALSE) }
  )
  
  # ----------------------
  # DOWNLOAD HANDLERS for Plots
  # ----------------------
  output$downloadLinPlot <- downloadHandler(
    filename = function() { paste("LinearPlot-", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 800, height = 600); on.exit(dev.off())
      mod <- linearModel(); dat <- dataInput()
      dose_seq <- seq(min(dat$dose), max(dat$dose), length.out = input$nPoints)
      newdata <- data.frame(dose = dose_seq)
      pred <- tryCatch({
        predict(mod, newdata, order = TRUE, exp = TRUE, ci.level = input$ciLevel, pi = TRUE)
      }, error = function(e) NULL)
      pred <- clean_pred(pred)
      if(is.null(pred)) return()
      plot(pred$dose, pred$pred, type = "l", col = "blue", lwd = input$lineWidth,
           ylim = range(c(pred$ci.lb, pred$ci.ub, pred$pi.lb, pred$pi.ub), na.rm = TRUE),
           xlab = "Dose (grams/day)", ylab = "Relative Risk",
           main = "Linear Model Predictions", bg = input$bgColor,
           cex.axis = input$textSize/12, cex.lab = input$textSize/12, cex.main = input$textSize/12)
      if(input$showCI && !is.null(pred$ci.lb)){
        lines(pred$dose, pred$ci.lb, lty = 2, col = "blue", lwd = input$lineWidth)
        lines(pred$dose, pred$ci.ub, lty = 2, col = "blue", lwd = input$lineWidth)
      }
      if(!is.null(pred$pi.lb)){
        lines(pred$dose, pred$pi.lb, lty = 3, col = "darkblue", lwd = input$lineWidth)
        lines(pred$dose, pred$pi.ub, lty = 3, col = "darkblue", lwd = input$lineWidth)
      }
      if(input$showRaw){
        points(dat$dose, exp(dat$logrr), pch = 16,
               col = if(input$themeChoice=="Black & White") "black" else "gray",
               cex = input$ptSize, alpha = input$alpha)
      }
      if(input$showGrid) grid()
      dev.off()
    }
  )
  
  output$downloadQuadPlot <- downloadHandler(
    filename = function() { paste("QuadraticPlot-", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 800, height = 600); on.exit(dev.off())
      mod <- quadraticModel(); dat <- dataInput()
      dose_seq <- seq(min(dat$dose), max(dat$dose), length.out = input$nPoints)
      newdata <- data.frame(dose = dose_seq)
      pred <- tryCatch({
        predict(mod, newdata, order = TRUE, exp = TRUE, ci.level = input$ciLevel, pi = TRUE)
      }, error = function(e) NULL)
      pred <- clean_pred(pred)
      if(is.null(pred)) return()
      plot(pred$dose, pred$pred, type = "l", col = "darkgreen", lwd = input$lineWidth,
           ylim = range(c(pred$ci.lb, pred$ci.ub, pred$pi.lb, pred$pi.ub), na.rm = TRUE),
           xlab = "Dose (grams/day)", ylab = "Relative Risk",
           main = "Quadratic Model Predictions", bg = input$bgColor,
           cex.axis = input$textSize/12, cex.lab = input$textSize/12, cex.main = input$textSize/12)
      if(input$showCI && !is.null(pred$ci.lb)){
        lines(pred$dose, pred$ci.lb, lty = 2, col = "darkgreen", lwd = input$lineWidth)
        lines(pred$dose, pred$ci.ub, lty = 2, col = "darkgreen", lwd = input$lineWidth)
      }
      if(!is.null(pred$pi.lb)){
        lines(pred$dose, pred$pi.lb, lty = 3, col = "darkgreen", lwd = input$lineWidth)
        lines(pred$dose, pred$pi.ub, lty = 3, col = "darkgreen", lwd = input$lineWidth)
      }
      if(input$showRaw){
        points(dat$dose, exp(dat$logrr), pch = 16,
               col = if(input$themeChoice=="Black & White") "black" else "gray",
               cex = input$ptSize, alpha = input$alpha)
      }
      if(input$showGrid) grid()
      dev.off()
    }
  )
  
  output$downloadSplinePlot <- downloadHandler(
    filename = function() { paste("SplinePlot-", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 800, height = 600); on.exit(dev.off())
      mod <- splineModel(); dat <- dataInput()
      dose_seq <- seq(min(dat$dose), max(dat$dose), length.out = input$nPoints)
      newdata <- data.frame(dose = dose_seq)
      pred <- tryCatch({
        predict(mod, newdata, xref = input$xref, order = TRUE, exp = TRUE, ci.level = input$ciLevel, pi = TRUE)
      }, error = function(e) NULL)
      pred <- clean_pred(pred)
      if(is.null(pred)) return()
      plot(pred$dose, pred$pred, type = "l", col = "purple", lwd = input$lineWidth,
           ylim = range(c(pred$ci.lb, pred$ci.ub, pred$pi.lb, pred$pi.ub), na.rm = TRUE),
           xlab = "Dose (grams/day)", ylab = "Relative Risk",
           main = "Spline Model Predictions", bg = input$bgColor,
           cex.axis = input$textSize/12, cex.lab = input$textSize/12, cex.main = input$textSize/12)
      if(input$showCI && !is.null(pred$ci.lb)){
        lines(pred$dose, pred$ci.lb, lty = 2, col = "purple", lwd = input$lineWidth)
        lines(pred$dose, pred$ci.ub, lty = 2, col = "purple", lwd = input$lineWidth)
      }
      if(!is.null(pred$pi.lb)){
        lines(pred$dose, pred$pi.lb, lty = 3, col = "purple", lwd = input$lineWidth)
        lines(pred$dose, pred$pi.ub, lty = 3, col = "purple", lwd = input$lineWidth)
      }
      if(input$showRaw){
        points(dat$dose, exp(dat$logrr), pch = 16,
               col = if(input$themeChoice=="Black & White") "black" else "gray",
               cex = input$ptSize, alpha = input$alpha)
      }
      if(input$showGrid) grid()
      dev.off()
    }
  )
  
  output$downloadResidPlot <- downloadHandler(
    filename = function() { paste("ResidualPlot-", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 800, height = 600); on.exit(dev.off())
      mod <- linearModel(); dat <- dataInput()
      res <- residuals(mod); fitted_vals <- mod$fitted.values
      n_obs <- min(length(res), length(fitted_vals), nrow(dat))
      res <- res[1:n_obs]; xvals <- dat$dose[1:n_obs]
      plot(xvals, res, xlab = "Dose (grams/day)", ylab = "Residuals",
           main = "Residuals vs Dose (Linear Model)", pch = 16, col = "red",
           cex.axis = input$textSize/12, cex.lab = input$textSize/12, cex.main = input$textSize/12)
      abline(h = 0, lty = 2)
      if(input$showGrid) grid()
      dev.off()
    }
  )
  
  output$downloadFitObsPlot <- downloadHandler(
    filename = function() { paste("FitObsPlot-", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 800, height = 600); on.exit(dev.off())
      mod <- linearModel(); dat <- dataInput()
      fitted_vals <- mod$fitted.values; observed <- exp(dat$logrr)
      n_obs <- min(length(fitted_vals), length(observed))
      df <- data.frame(Fitted = fitted_vals[1:n_obs],
                       Observed = observed[1:n_obs])
      p <- ggplot(df, aes(x = Fitted, y = Observed)) +
        geom_point(color = "blue", size = 3) +
        geom_smooth(method = "loess", se = FALSE, color = "black") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
        labs(x = "Fitted Relative Risk", y = "Observed Relative Risk",
             title = "Fitted vs. Observed (Linear Model)") +
        theme_minimal(base_size = input$textSize) +
        theme(panel.grid.major = if(input$showGrid) element_line() else element_blank(),
              panel.grid.minor = if(input$showGrid) element_line() else element_blank(),
              plot.background = element_rect(fill = input$bgColor))
      print(p)
      dev.off()
    }
  )
  
  output$downloadQQPlot <- downloadHandler(
    filename = function() { paste("QQPlot-", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 800, height = 600); on.exit(dev.off())
      mod <- linearModel(); res <- residuals(mod)
      qqnorm(res, main = "Q-Q Plot of Residuals (Linear Model)",
             cex.axis = input$textSize/12, cex.lab = input$textSize/12, cex.main = input$textSize/12)
      qqline(res, col = "red", lty = 2)
      dev.off()
    }
  )
  
  output$downloadHistPlot <- downloadHandler(
    filename = function() { paste("HistPlot-", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 800, height = 600); on.exit(dev.off())
      mod <- linearModel(); res <- residuals(mod)
      hist(res, breaks = 20, main = "Histogram of Residuals (Linear Model)",
           xlab = "Residuals", col = "lightblue", border = "white",
           cex.axis = input$textSize/12, cex.lab = input$textSize/12, cex.main = input$textSize/12)
      dev.off()
    }
  )
  
  output$downloadBubblePlot <- downloadHandler(
    filename = function() { paste("BubblePlot-", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 800, height = 600); on.exit(dev.off())
      dat <- dataInput()
      p <- ggplot(dat, aes(x = dose, y = logrr, size = 1/se, color = id)) +
        geom_point(alpha = input$alpha) +
        scale_size_area(max_size = 10) +
        labs(x = "Dose (grams/day)", y = "Log Relative Risk", size = "1/SE", color = "Study ID") +
        theme_minimal(base_size = input$textSize) +
        theme(panel.grid.major = if(input$showGrid) element_line() else element_blank(),
              panel.grid.minor = if(input$showGrid) element_line() else element_blank(),
              plot.background = element_rect(fill = input$bgColor)) +
        ggtitle("Bubble Plot of Raw Data")
      print(p)
      dev.off()
    }
  )
  
  output$downloadStudyTrendPlot <- downloadHandler(
    filename = function() { paste("StudyTrendPlot-", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 800, height = 600); on.exit(dev.off())
      dat <- dataInput()
      p <- ggplot(dat, aes(x = dose, y = logrr)) +
        geom_point(size = input$ptSize, color = if(input$themeChoice=="Black & White") "black" else "darkblue") +
        geom_line(color = if(input$themeChoice=="Black & White") "black" else "darkblue") +
        facet_wrap(~ id) +
        labs(x = "Dose (grams/day)", y = "Log Relative Risk") +
        theme_minimal(base_size = input$textSize) +
        theme(panel.grid.major = if(input$showGrid) element_line() else element_blank(),
              panel.grid.minor = if(input$showGrid) element_line() else element_blank(),
              plot.background = element_rect(fill = input$bgColor)) +
        ggtitle("Study-Specific Trends")
      print(p)
      dev.off()
    }
  )
  
  output$downloadForestPlot <- downloadHandler(
    filename = function() { paste("ForestPlot-", Sys.Date(), ".png", sep="") },
    content = function(file) {
      png(file, width = 800, height = 600); on.exit(dev.off())
      dat <- dataInput(); fdat <- make_forest_data(dat)
      p <- ggplot(fdat, aes(x = RR, y = reorder(as.factor(id), as.numeric(id)))) +
        geom_point(size = 3, color = "darkred") +
        geom_errorbarh(aes(xmin = RR_lb, xmax = RR_ub), height = 0.2, color = "darkred") +
        labs(x = "Relative Risk", y = "Study ID", title = "Forest Plot (Raw Data)") +
        theme_minimal(base_size = input$textSize) +
        theme(panel.grid.major = if(input$showGrid) element_line() else element_blank(),
              panel.grid.minor = if(input$showGrid) element_line() else element_blank(),
              plot.background = element_rect(fill = input$bgColor))
      print(p)
      dev.off()
    }
  )
  
  # ---------------------------------------------------------
  # ADVANCED ANALYSIS SECTION
  # ---------------------------------------------------------
  
  # 1. Heterogeneity & Model Diagnostics
  output$hetStats <- renderPrint({
    lin <- linearModel(); quad <- quadraticModel(); spl <- splineModel()
    list(
      Linear    = c(Q = lin$QE, df = lin$df.Q, pval = lin$QEp, I2 = lin$I2),
      Quadratic = c(Q = quad$QE, df = quad$df.Q, pval = quad$QEp, I2 = quad$I2),
      Spline    = c(Q = spl$QE, df = spl$df.Q, pval = spl$QEp, I2 = spl$I2)
    )
  })
  
  # 2. Subgroup Analysis
  output$subgroupUI <- renderUI({
    dat <- dataInput()
    vars <- setdiff(names(dat), c("id", "dose", "logrr", "se", "cases", "n", "type"))
    if(length(vars) == 0){
      helpText("No extra columns for subgroup analysis.")
    } else {
      selectInput("subgroupVar", "Subgroup Variable:", choices = vars)
    }
  })
  
  observeEvent(input$runSubgroup, {
    req(input$subgroupVar)
    dat <- dataInput(); gvar <- input$subgroupVar
    subs <- split(dat, dat[[gvar]])
    results <- lapply(names(subs), function(sg) {
      subdat <- subs[[sg]]
      if(nrow(subdat) < 5) {
        return(list(group = sg, linear = NULL, quad = NULL, spline = NULL))
      }
      lm_model <- tryCatch(
        dosresmeta(logrr ~ dose, id = id, se = se, type = type,
                   cases = cases, n = n, method = input$method, data = subdat),
        error = function(e) NULL
      )
      qm_model <- tryCatch(
        dosresmeta(logrr ~ dose + I(dose^2), id = id, se = se, type = type,
                   cases = cases, n = n, method = input$method, data = subdat),
        error = function(e) NULL
      )
      knots <- tryCatch(quantile(subdat$dose, probs = c(0.05, 0.35, 0.65, 0.95), na.rm = TRUE),
                        error = function(e) NULL)
      sm_model <- NULL
      if(!is.null(knots)) {
        sm_model <- tryCatch(
          dosresmeta(logrr ~ rcs(dose, knots), id = id, se = se, type = type,
                     cases = cases, n = n, method = input$method, data = subdat),
          error = function(e) NULL
        )
      }
      list(group = sg, linear = lm_model, quad = qm_model, spline = sm_model)
    })
    txt <- ""
    for(r in results){
      txt <- paste0(txt, "\n---- Subgroup: ", r$group, " ----\n")
      if(!is.null(r$linear))
        txt <- paste0(txt, "[Linear]\n", capture.output(summary(r$linear)), "\n")
      else
        txt <- paste0(txt, "[Linear] Not enough data or error.\n")
      if(!is.null(r$quad))
        txt <- paste0(txt, "[Quadratic]\n", capture.output(summary(r$quad)), "\n")
      else
        txt <- paste0(txt, "[Quadratic] Not enough data or error.\n")
      if(!is.null(r$spline))
        txt <- paste0(txt, "[Spline]\n", capture.output(summary(r$spline)), "\n")
      else
        txt <- paste0(txt, "[Spline] Not enough data or error.\n")
      txt <- paste0(txt, "------------------------------------\n")
    }
    updateTextInput(session, "dummySubgroupTxt", value = txt)
  })
  
  output$subgroupResults <- renderText({ input$dummySubgroupTxt })
  observe({ if(is.null(input$dummySubgroupTxt)) updateTextInput(session, "dummySubgroupTxt", value = "") })
  outputOptions(output, "subgroupResults", suspendWhenHidden = FALSE)
  
  # 3. Quartile-Based Analysis
  quartileData <- eventReactive(input$runQuartile, {
    dat <- dataInput()
    dat$quartile <- cut(dat$dose, 
                        breaks = quantile(dat$dose, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
                        include.lowest = TRUE,
                        labels = c("Q1", "Q2", "Q3", "Q4"))
    dfq <- dat %>% group_by(id, quartile) %>% summarize(logrr = mean(logrr), se = mean(se), .groups = "drop")
    dfq
  })
  
  output$quartileForest <- renderPlot({
    dfq <- quartileData()
    if(is.null(dfq) || nrow(dfq) < 1) return(NULL)
    meta_df <- dfq %>% group_by(quartile) %>% summarize(m_logrr = mean(logrr),
                                                        se_pooled = sd(logrr) / sqrt(n()),
                                                        n_studies = n())
    meta_df$RR <- exp(meta_df$m_logrr)
    meta_df$lb <- exp(meta_df$m_logrr - 1.96 * meta_df$se_pooled)
    meta_df$ub <- exp(meta_df$m_logrr + 1.96 * meta_df$se_pooled)
    ggplot(meta_df, aes(x = RR, y = quartile)) +
      geom_point(size = 3) +
      geom_errorbarh(aes(xmin = lb, xmax = ub), height = 0.2) +
      labs(title = "Quartile-Based (Naive) Summary", x = "Relative Risk", y = "Quartile")
  })
  
  # 4. Publication Bias Checks
  pubBiasData <- eventReactive(input$runPubBias, {
    dat <- dataInput()
    sums <- dat %>% group_by(id) %>% summarize(
      minLog = logrr[dose == min(dose)][1],
      minSe  = se[dose == min(dose)][1],
      maxLog = logrr[dose == max(dose)][1],
      maxSe  = se[dose == max(dose)][1],
      .groups = "drop")
    sums$diffLog <- sums$maxLog - sums$minLog
    sums$diffVar <- sums$maxSe^2 + sums$minSe^2
    sums
  })
  
  output$funnelPlot <- renderPlot({
    pb <- pubBiasData(); pb <- na.omit(pb)
    if(nrow(pb) < 2) {
      plot.new(); text(0.5, 0.5, "Not enough data for funnel plot", cex = 1.2)
      return(NULL)
    }
    res <- rma(yi = pb$diffLog, vi = pb$diffVar, method = "REML")
    funnel(res, main = "Funnel Plot (Highest vs. Lowest Dose)")
  })
  
  output$eggerTest <- renderPrint({
    pb <- pubBiasData(); pb <- na.omit(pb)
    if(nrow(pb) < 2) return("Not enough data for Egger's test.")
    res <- rma(yi = pb$diffLog, vi = pb$diffVar, method = "REML")
    regtest(res, model = "rma")
  })
  
  # 5. Leave-One-Out Sensitivity (Linear Model)
  observeEvent(input$runLoo, {
    dat <- dataInput(); all_studies <- unique(dat$id)
    results <- lapply(all_studies, function(study) {
      subdat <- subset(dat, id != study)
      if(nrow(subdat) < 5) {
        return(data.frame(StudyOut = as.character(study),
                          Coef = NA, SE = NA,
                          Error = "Too few data after removal",
                          row.names = NULL))
      }
      mod <- tryCatch({
        dosresmeta(logrr ~ dose, id = id, se = se, type = type,
                   cases = cases, n = n, data = subdat, method = input$method)
      }, error = function(e) e)
      if(inherits(mod, "error")){
        return(data.frame(StudyOut = as.character(study),
                          Coef = NA, SE = NA,
                          Error = mod$message,
                          row.names = NULL))
      } else {
        coefs <- coef(mod); ses <- sqrt(diag(vcov(mod)))
        data.frame(StudyOut = as.character(study),
                   Coef = coefs["dose"],
                   SE   = ses["dose"],
                   Error = NA,
                   row.names = NULL)
      }
    })
    results_df <- do.call(rbind, results)
    row.names(results_df) <- NULL
    txt <- capture.output(print(results_df))
    updateTextInput(session, "dummyLOOTxt", value = paste(txt, collapse = "\n"))
  })
  
  output$looResults <- renderText({ input$dummyLOOTxt })
  observe({ if(is.null(input$dummyLOOTxt)) updateTextInput(session, "dummyLOOTxt", value = "") })
  outputOptions(output, "looResults", suspendWhenHidden = FALSE)
  
  # 6. Meta-Regression Analysis (Extra Moderator(s))
  observeEvent(input$runMetaReg, {
    req(length(input$metaRegVars) > 0 && input$metaRegVars[1] != "None")
    dat <- dataInput()
    # Build moderator formula: using vi = se^2
    mods_formula <- as.formula(paste("~", paste(input$metaRegVars, collapse = " + ")))
    mod_reg <- tryCatch({
      rma(yi = logrr, vi = se^2, mods = mods_formula, data = dat, method = "REML")
    }, error = function(e) e)
    if(inherits(mod_reg, "error")){
      output$metaRegSummary <- renderPrint({ paste("Error:", mod_reg$message) })
    } else {
      output$metaRegSummary <- renderPrint({ summary(mod_reg) })
    }
  })
  
  # 7. Custom Spline Model Analysis with User-Defined Knots
  observeEvent(input$runCustomSpline, {
    dat <- dataInput()
    percentiles <- as.numeric(unlist(strsplit(input$customKnotsPercentiles, ","))) / 100
    custom_knots <- quantile(dat$dose, probs = percentiles, na.rm = TRUE)
    custom_spline <- tryCatch({
      dosresmeta(logrr ~ rcs(dose, custom_knots), id = id, se = se, type = type,
                 cases = cases, n = n, data = dat, method = input$method)
    }, error = function(e) e)
    if(inherits(custom_spline, "error")){
      output$customSplineSummary <- renderPrint({ paste("Error:", custom_spline$message) })
      output$customSplinePlot <- renderPlot({ plot.new(); text(0.5, 0.5, "Custom Spline Plot Error", cex = 1.2) })
    } else {
      output$customSplineSummary <- renderPrint({ summary(custom_spline) })
      output$customSplinePlot <- renderPlot({
        dose_seq <- seq(min(dat$dose), max(dat$dose), length.out = input$nPoints)
        newdata <- data.frame(dose = dose_seq)
        pred <- tryCatch({
          predict(custom_spline, newdata, xref = input$xref, order = TRUE,
                  exp = TRUE, ci.level = input$ciLevel, pi = TRUE)
        }, error = function(e) NULL)
        pred <- clean_pred(pred)
        if(is.null(pred)){
          plot.new(); text(0.5, 0.5, "No valid predictions", cex = 1.2)
          return()
        }
        plot(pred$dose, pred$pred, type = "l", col = "orange", lwd = input$lineWidth,
             ylim = range(c(pred$ci.lb, pred$ci.ub, pred$pi.lb, pred$pi.ub), na.rm = TRUE),
             xlab = "Dose (grams/day)", ylab = "Relative Risk",
             main = "Custom Spline Model Predictions")
        if(input$showCI && !is.null(pred$ci.lb)){
          lines(pred$dose, pred$ci.lb, lty = 2, col = "orange", lwd = input$lineWidth)
          lines(pred$dose, pred$ci.ub, lty = 2, col = "orange", lwd = input$lineWidth)
        }
        if(!is.null(pred$pi.lb)){
          lines(pred$dose, pred$pi.lb, lty = 3, col = "orange", lwd = input$lineWidth)
          lines(pred$dose, pred$pi.ub, lty = 3, col = "orange", lwd = input$lineWidth)
        }
        if(input$showRaw){
          points(dat$dose, exp(dat$logrr), pch = 16, col = "gray")
        }
        grid()
      })
    }
  })
  
}

shinyApp(ui, server)
