# Explore any within-cohort comparison interactively.
# Open this file in RStudio and click "Run App" (top right of the editor).
# Needs results/analysis_data.rds, written at the end of notebook 5.
library(shiny)
source("config.R")
source("R/functions.R")

dat <- readRDS(file.path(RESULTS, "analysis_data.rds"))
GROUPING <- c("Diagnostic status" = "dx_status", "Gene group" = "gene_group", "Variant class" = "variant_class",
              "Gene" = "gene", "Epilepsy type" = "epilepsy_type", "Drug resistant" = "drug_resistant", "Sex" = "sex")
pt_labels <- sort(unique(as.numeric(sub("^epilepsy_Pt_(.*)_z$", "\\1", grep("^epilepsy_Pt_.*_z$", names(dat), value = TRUE)))))

ui <- fluidPage(
  titlePanel("PRS: compare any two groups within the cohort"),
  sidebarLayout(
    sidebarPanel(width = 4,
      selectInput("trait", "Polygenic score", c("Epilepsy" = "epilepsy", "Height (negative control)" = "height")),
      selectInput("pt", "GWAS P value threshold", pt_labels, selected = PT_MAIN),
      selectInput("var", "Define the groups by", GROUPING),
      uiOutput("level_pickers"),
      checkboxInput("restrict", "Only people with a monogenic diagnosis", FALSE),
      checkboxGroupInput("covs", "Adjust for", c("Sex" = "sex", "Ancestry PCs" = "PCs"), selected = c("sex", "PCs")),
      helpText("Every extra comparison you try is another test: treat what you find here as exploration,",
               "and write down your hypotheses before the real analysis.")
    ),
    mainPanel(width = 8,
      plotOutput("plot", height = "420px"),
      h4(textOutput("verdict")),
      tableOutput("table")
    )
  )
)

server <- function(input, output, session) {
  output$level_pickers <- renderUI({
    lv <- sort(unique(na.omit(dat[[input$var]])))
    tagList(
      checkboxGroupInput("case_lv", "Case group", lv, selected = lv[1]),
      checkboxGroupInput("ctrl_lv", "Comparison group", lv, selected = lv[min(2, length(lv))])
    )
  })

  sub <- reactive({
    d <- if (input$restrict) filter(dat, dx_status == "Monogenic") else dat
    d$.score <- d[[score_col(d, input$trait, as.numeric(input$pt))]]
    d
  })

  covariates <- reactive({
    cv <- character()
    if ("sex" %in% input$covs) cv <- c(cv, "sex")
    if ("PCs" %in% input$covs) cv <- c(cv, paste0("PC", 1:N_PCS))
    if (length(cv) == 0) "1" else cv
  })

  result <- reactive({
    req(input$case_lv, input$ctrl_lv)
    validate(need(length(intersect(input$case_lv, input$ctrl_lv)) == 0,
                  "A level cannot be in both groups: untick it in one of them."))
    var <- input$var; case_lv <- input$case_lv; ctrl_lv <- input$ctrl_lv
    res <- tryCatch(
      compare_groups(sub(), .data[[var]] %in% case_lv, .data[[var]] %in% ctrl_lv, prs = ".score",
                     covariates = covariates(),
                     name = paste(paste(case_lv, collapse = " + "), "vs", paste(ctrl_lv, collapse = " + "))),
      error = function(e) e)
    validate(need(!inherits(res, "error"), if (inherits(res, "error")) conditionMessage(res) else ""))
    res
  })

  output$plot <- renderPlot({
    req(result())
    var <- input$var; case_lv <- input$case_lv; ctrl_lv <- input$ctrl_lv
    plot_groups(sub(), .data[[var]] %in% case_lv, .data[[var]] %in% ctrl_lv, prs = ".score",
                covariates = covariates(), title = result()$contrast,
                labels = c(paste(ctrl_lv, collapse = " + "), paste(case_lv, collapse = " + "))) +
      labs(y = paste(tools::toTitleCase(input$trait), "PRS at P <", input$pt, "(SD, adjusted)"))
  })

  output$verdict <- renderText({
    r <- result()
    sprintf("The case group's PRS is %.2f SD %s (95%% CI %.2f to %.2f), P = %s.",
            abs(r$diff_SD), ifelse(r$diff_SD > 0, "higher", "lower"), r$low, r$high, signif(r$P, 2))
  })

  output$table <- renderTable(result() %>% select(-prs), digits = 3)
}

shinyApp(ui, server)
