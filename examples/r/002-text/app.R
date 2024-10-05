#library(plotly)
#library(dplyr)
#library(shiny)
#library(combinat)
#library(htmltools)
#library(ggplot2)
#library(gridExtra)
#library(broom)
#library(tidyr)
#library(combinat)
#library(AlgebraicHaploPackage)
#library(huxtable)
#library(kableExtra)
#library(polynom)
ManifoldDestiny::wasmconload()
#Shiny
ui <- fluidPage(
  titlePanel("R2 Simulator"),
  tabsetPanel(
    tabPanel("Graph", 
             sidebarLayout(
               sidebarPanel(
                 selectInput("variable", "Select a form:", c("Normal form" = "normal", "Hybrid form" = "hybrid", "Opposition form" = "opposition")),
                 numericInput("nprec", "Number of precincts:", value = 300),
                 textInput("regv", "Registered voters (mean,std)", value='3.15, 0.25'),
                 textInput("minmax", "Min and max values number of voters", value='400, 4000'),
                 textInput("turn", "Probability of voting:", value = '0.5, 0.10'),
                 textInput("invper", "Election system variability:", value ='0.5, 0.10'),
                 textInput("u", "Probability of x/g/n:", value = '0.6, 0.10'),
                 textInput("dv", "Diff probability of x/g/n:", value = '-0.2, 0.08'),
                 numericInput("draws", "Number of draws", value = 30),
                 actionButton("run", "Run Simulation")
               ),
               mainPanel(
                 plotOutput("plot"),
                 tableOutput("table")
               )
             ))
    ),
    tabPanel("Tables", 
             tabsetPanel(
               tabPanel("Table 1",
                        sidebarLayout(
                          sidebarPanel(),
                          mainPanel(
                            tableOutput("table1")
                          )
                        )
               ),
               tabPanel("Table 2",
                        sidebarLayout(
                          sidebarPanel(),
                          mainPanel(
                            tableOutput("table2")
                          )
                        )
               )
             )
    )
)
server <- function(input, output) {
  # Run simulation
  dfgp <- eventReactive(input$run, {
    # Input DF1
    regv <- as.numeric(strsplit(input$regv, ",")[[1]])
    minmax <- as.numeric(strsplit(input$minmax, ",")[[1]])
    turn <- as.numeric(strsplit(input$turn, ",")[[1]])
    invper <- as.numeric(strsplit(input$invper, ",")[[1]])
    u <- as.numeric(strsplit(input$u, ",")[[1]])
    dv <- as.numeric(strsplit(input$dv, ",")[[1]])
    tf <- replicate(input$draws,r2simn(nprec = input$nprec,
          regs = c(regv[1], regv[2]),
          minmax = c(minmax[1],minmax[2]), turn = c(turn[1], turn[2]), Invper = c(invper[1], invper[2]),
          u = c(u[1], u[2]),
          dv = c(dv[1], dv[2]),
          form = 1)[c(1,2,3)])
    dfgp <- data.frame(r2a=unlist(tf[seq(1,length(tf),3)]),r2b=unlist(tf[seq(2,length(tf),3)])) %>% mutate(perc = ntile(r2a, 100)) 
    # Input DF2
    percentiles <- c(90, 95, 99)
    nstd <- c(1,2,5)
    std <- mean(dfgp$r2a)+nstd*sd(dfgp$r2a)
    perc1 <- quantile(dfgp$r2a,probs = percentiles / 100)
    perc2 <- quantile(dfgp$r2b,probs = percentiles / 100)
    percdf <- data.frame(perc1,perc2,nstd,std) %>% data.table::setnames(c("Perc r2a","Perc r2b","Nstd","Vstd")) 
    list(dfgp,percdf)
  })
  # Create plot
  output$plot <- renderPlot({
	  #browser()
    dfp <- dfgp()[[1]] %>% tidyr::pivot_longer(cols=c("r2a","r2b")) %>% dplyr::arrange(name,perc)
    ggplot(dfp,aes(x=value, fill=name)) + 
      geom_histogram(position = "identity", alpha = 0.5, bins = 30) + 
      labs(title = "Histogram of Values by Category", x = "Value", y = "Count") +
      geom_vline(xintercept = as.numeric(dfgp()[[2]][1,1]), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = as.numeric(dfgp()[[2]][2,1]), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = as.numeric(dfgp()[[2]][3,1]), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = as.numeric(dfgp()[[2]][3,4]), linetype = "solid", color = "red") +
      geom_label(y=0,x=as.numeric(dfgp()[[2]][1,1]),label="*",geom="label") +
      geom_label(y=0,x=as.numeric(dfgp()[[2]][2,1]),label="**",geom="label") +
      geom_label(y=0,x=as.numeric(dfgp()[[2]][3,1]),label="***",geom="label") +
      theme_minimal() +
      scale_fill_manual(values = c("#0072B2", "#E69F00"))  # set fill colors
  })
  
  ## Create table
  output$table1 <- renderUI({
    DT::datatable(round(dfgp()[[2]],digits=4))
  })
  output$table2 <- renderUI({
    DT::datatable(round(dfgp()[[1]], digits=4), options = list(pageLength = 20))
  })
}
# Run app
shinyApp(ui = ui, server = server)

