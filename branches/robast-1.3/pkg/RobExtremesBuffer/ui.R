shinyUI(fluidPage(
  shinyjs::useShinyjs(),
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
      sliderInput("whichLM", label="L-Multiplikator auswählen", min=1, max=10, value=1),
      
      checkboxInput("takeUsed", label="Gespeicherten Spline verwenden."),
      numericInput("df", "DF", 10),
      selectInput("ranges", label="Glättung-Ausschluss-Intervalle", choices=NULL, size=10, selectize=FALSE),
      
      fluidRow(
        column(5, actionButton("deleteRange", label="Interval löschen"))
      ),
      
      fluidRow(
        column(5, actionButton("saveGrid", label="Speichere Grid (CSV)"))
      ),
      
      fluidRow(
        column(5, actionButton("addToHistory", label="Zu History hinzufügen")),
        column(5, actionButton("loadFromHistory", label="Aus History laden"))
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
