library(shiny)
library(ggplot2)
library(purrr)

# ==========================================
# 1. UNDERLYING LOGIC (From your R script)
# ==========================================
ncp_est = function(N_inp, beta_inp, maf_inp, cc_inp){
  return(2*maf_inp*(1-maf_inp)*N_inp*cc_inp*(1-cc_inp)*beta_inp^2)
}

pwr_est = function(thresh, N, OR, maf, cc_rat){
  beta = log(OR)
  pwr = pchisq(thresh, df = 1, ncp = ncp_est(N, beta, maf, cc_rat), lower = FALSE)
  return(cbind.data.frame("Sample_Size" = N, "Power" = pwr, "OR" = OR, "CC_ratio" = cc_rat))
}

# ==========================================
# 2. USER INTERFACE (UI)
# ==========================================
ui <- fluidPage(
  theme = theme_minimal(),
  titlePanel("Genetic Association Study Power Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      # Slider for MAF parameter
      numericInput("maf", "Minor Allele Frequency (MAF):", 
                  min = 0.001, max = 0.5, value = 0.01, step = 0.005),
      
      # Numeric input for custom significance thresholds
      numericInput("alpha", "Significance Alpha Threshold:", 
                   value = 1e-4, min = 1e-8, max = 0.05, step = 1e-4),
      
      hr(),
      tags$b("Simulation Grid Ranges:"),
      p("Adjust the slider sequences to recalculate the power matrix."),
      
      # Dynamic sequences mapping to your grid variables
      sliderInput("sample_sizes", "Sample Size Range (N):", 
                  min = 50, max = 10000, value = c(100, 3000), step = 50),
      
      sliderInput("odds_ratios", "Odds Ratio Range (OR):", 
                  min = 1.0, max = 10.0, value = c(1.0, 4.0), step = 0.25),
      
      # Dropdown to isolate specific OR views for plotting
      hr(),
      uiOutput("or_selector")
    ),
    
    mainPanel(
      # Visual layout tabs for outputs
      tabsetPanel(
        tabPanel("Power Curves", plotOutput("power_plot")),
        tabPanel("Raw Simulation Data", dataTableOutput("data_table"))
      )
    )
  )
)

# ==========================================
# 3. SERVER LOGIC
# ==========================================
server <- function(input, output, session) {
  
  # Reactive calculations: runs matrix grid simulation only when settings change
  sim_data <- reactive({
    # Calculate threshold based on active UI alpha choice
    sig_thres <- qchisq(input$alpha, df = 1, ncp = 0, lower.tail = FALSE)
    
    # Establish dynamic vectors based on user inputs
    ns <- seq(input$sample_sizes[1], input$sample_sizes[2], 50)
    ors <- seq(input$odds_ratios[1], input$odds_ratios[2], 0.25)
    ccs <- seq(0.1, 0.5, 0.1) # Kept standard case-control ratio range
    
    df_sim <- expand.grid(Var1 = ns, Var2 = ors, Var3 = ccs)
    
    # Grid calculation mapping
    res_cmb <- pmap(df_sim, function(Var1, Var2, Var3) {
      pwr_est(sig_thres, N = Var1, OR = Var2, maf = input$maf, cc_rat = Var3)
    })
    
    do.call("rbind.data.frame", res_cmb)
  })
  
  # Dynamic UI selector to target an individual OR filter value for plotting
  output$or_selector <- renderUI({
    df <- sim_data()
    selectInput("selected_or", "Select OR to View in Graph:", 
                choices = sort(unique(df$OR)), 
                selected = max(df$OR))
  })
  
  # Render the resulting ggplot
  output$power_plot <- renderPlot({
    req(input$selected_or)
    
    # Filter the reactive dataset based on chosen OR slice
    plot_df <- sim_data()
    plot_df <- plot_df[plot_df$OR == as.numeric(input$selected_or), ]
    
    ggplot(plot_df, aes(x = Sample_Size, y = Power, colour = factor(CC_ratio))) + 
      geom_line(linewidth = 1) + 
      labs(x = "Sample Size (N)", y = "Statistical Power", colour = "Case-Control Ratio") +
      theme_minimal() +
      ylim(0, 1)
  })
  
  # Render raw preview data table tab
  output$data_table <- renderDataTable({
    sim_data()
  }, options = list(pageLength = 10))
}

# ==========================================
# 4. APP LAUNCHER
# ==========================================
shinyApp(ui = ui, server = server)