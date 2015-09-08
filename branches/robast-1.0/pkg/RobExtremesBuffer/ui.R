shinyUI(fluidPage(
  titlePanel("Smoother"),
  
  fluidRow(
    column(1, 
      radioButtons("gridName", 
                  label = "Grid auswählen",
                  choices = c("Sn", "OMSE", "RMXE", "MBRE")
                  ),
      
      radioButtons("familyName", 
                   label = "Familie auswählen",
                   choices = c("Generalized Pareto Family", 
                               "GEV Family",
                               "GEVU Family",
                               "Gamma family",
                               "Weibull Family")
                   ),
      htmlOutput("resetNote")
    ), # column
    
    column(3, 
      # checkboxInput("withSmooth", label="Plot with Smooth"),
      sliderInput("whichLM", label="L-Multiplikator auswählen", min=1, max=10, value=1),
      
      numericInput("df", "DF", 10),
      selectInput("ranges", label="Glättung-Ausschluss-Intervalle", choices=NULL, size=10, selectize=FALSE),
      actionButton("deleteRange", label="Interval löschen"),
      
      fluidRow( 
        column(6, actionButton("saveGrid", label="Speichere Grid (CSV)")),
        column(6, actionButton("addToHistory", label="Zu History hinzufügen"))
      )
    ), # column
  
    column(8,
      plotOutput("out", brush=brushOpts("plot_brush", delay=300, resetOnNew=TRUE), # brush for zoom
                        dblclick="plot_dblclick", # double click to select ranges
                        height="600px"
                 ), 
      fluidRow(
        column(1, actionButton("zoomOut", label="Zoom Out", icon=icon("zoom-out", lib="glyphicon"))),
        column(2, offset=1, checkboxInput("withLegend", label="Legende anzeigen", value=TRUE))
      )
    ) # column
  ) # fluidRow
))