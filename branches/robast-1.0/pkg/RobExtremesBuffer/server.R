source("config.R")

options(shiny.error=NULL)
options(shiny.trace=F)


source("utilities.R")
checkRequiredPackages()

RESET_NOTE_TEXT <- "<strong>Durch aendern der Grid & Familie, gehen alle nicht-gespeicherte Intervalle verloren.</strong>"
DEFAULT_DEGREE_OF_FREEDOM <- 10

loadedData <- loadRDataFileToEnv("sysdata.rda")
zoomHistory <<- NULL
prev.deleted <<- ""
# We have grid as global, since we want to do testing.
smoothed.totalgrid <<- NULL



shinyServer(function(input, output, session){
######################################################################################
## Set parameters
  output$resetNote <- renderText(RESET_NOTE_TEXT)
  
  ## prepare data
  loadGrids <- reactive({ # Depends on input${gridName, familyName}
    loadGridsIntoEnv(loadedData, gridName=getCurrentGridName(), familyName=getCurrentFamilyName())
  })

######################################################################################
## get current configuration
  getCurrentLM <- reactive({
    lm = input$whichLM
    
    if (isSnGridWithInvalidLM(getCurrentGridName(), lm)){
      updateSliderInput(session, inputId="whichLM", value=1)
      return(1)
    }else{
      return(lm)
    }
  })
  
  getCurrentGridName <- reactive({return(input$gridName)})
  getCurrentFamilyName <- reactive({return(input$familyName)})

######################################################################################
## Construct plot parameters
  observe({ # Depends on getCurrentGridName(), getCurrentFamilyName(). Sets input$whichLM
    maxNumMultiplicators <- getNumMultiplicators(getCurrentGridName(), getCurrentFamilyName())
    updateSliderInput(session, inputId="whichLM", max=maxNumMultiplicators)
    isLoadFromHistoryEnabled <- hasHistory(getCurrentFamilyName(), getCurrentGridName())
    toogleAvailabilityComponents("loadFromHistory", isLoadFromHistoryEnabled)
  }, label="form components setting on grid & family")
  
  getOriginalGrid <- reactive({
    res <- loadGrids()[["orig"]]
    return(res)
    }) # Depends on loadGrids()
  getOriginalSmoothedGrid <- reactive({loadGrids()[["smoothed"]]}) # Depends on loadGrids()
  
  getEditingGrid <- reactive({return(getOriginalGrid())})
  
  getPostSmoothedEditingGrid <- reactive({ # Depends on getEditingGrid(), getCurrentSmoothRestrictions(), getCurrentDf()
    grid <- getEditingGrid()
    result <- applySmoothing(grid=grid, df=getCurrentDf(), gridRestrictions=getCurrentSmoothRestrictions())
    return(result)
  })
  
  plotGridRestriction <- reactive({
    return(list(rep(TRUE, nrow(getPostSmoothedEditingGrid()))))
  }) # Depends on getPostSmoothedEditingGrid()
  
  
######################################################################################
## smooth restrictions
  
  # Represents three states: 
  # (1) no current restriction (i.e. currently: both NULL; first will get some value)
  # (2) restriction started    (i.e. currently: first HAS some value; second will get)
  # (3) restriction finished   (i.e  currently: both HAS some value, will remove all)
  #
  # Removing will be done by clicking on existing first.  (i.e. currently: first HAS some value; it will be removed)
  # Then we will get to state (1)
  # The sw
  #
  #
  # So the automata look like: (-1) -> (1) <-> (2)
  #                                     \--<<--/
  beginRestrictionInterval <<- NULL # We do not want them as responsive events.
  restrictionState <- reactiveValues(state=NULL)
  configuration <- reactiveValues(ranges=NULL, df=NULL, useExisting=NULL)
  
  getCurrentDf <- reactive({ # Depends on configuration$df, getCurrentLM()
    return(configuration$df[[getCurrentLM()]])
  })
  
  getCurrentUseExisting <- reactive({ # Depends on configuration${useExisting, df}
    return(configuration$useExisting[[getCurrentLM()]])
  })
  
  getCurrentSmoothRestrictions <- reactive({  # Depends on configuration$ranges, getEditingGrid(), getCurrentLM()
    result <- get.restrictions.for.smooth(getCurrentLM(), configuration$ranges, getEditingGrid()[,1])
    return(result)
  })
  
  getInputDf <- reactive({
    if(is.na(input$df))
      return(1)
    
    return(input$df)
    })
  
  resetConfiguration <- function(){
    numMultiplicators <- getNumMultiplicators(isolate(getCurrentGridName()), isolate(getCurrentFamilyName()))
    
    beginRestrictionInterval <<- NULL
    
    configuration$ranges <- create.list.of.empty.lists(numMultiplicators)
    configuration$df <- as.list(rep(DEFAULT_DEGREE_OF_FREEDOM, numMultiplicators))
    configuration$useExisting <- as.list(rep(FALSE, numMultiplicators))
  }
  
  
  observe({ # Depends on getInputDf(). Outputs output$df, configuration$df
    if(getInputDf() > 0) {
      whichLM <- isolate(getCurrentLM())
      configuration$df[[whichLM]] <- getInputDf()
    } else {
      updateNumericInput(session, "df", value=1)
    }
  }, label="update df")
  
  
  observe({ # input$takeUsed, getCurrentLM()
    COMPONENTS <- c('df', 'ranges', 'deleteRange')
    whichLM <- isolate(getCurrentLM())
    
    configuration$useExisting[[whichLM]] <- input$takeUsed
    
    toogleAvailabilityComponents(COMPONENTS, !input$takeUsed)
  }, label="Use saved grid checkbox")
  
  
  observe({ # depends on: configuration$ranges. Sets: output$ranges, getCurrentLM()
    lm <- getCurrentLM()
    df <- getCurrentDf()
    useExisting <- getCurrentUseExisting()
    
    ranges <- configuration$ranges[[lm]] # TODO: Create getCurrentRanges() function
    
    if(df != isolate(getInputDf())){
      updateNumericInput(session, "df", value=df)
    }
    
    if(useExisting != isolate(input$takeUsed)) {
      updateCheckboxInput(session, 'takeUsed', value=useExisting)
    }
    
    updateSelectInput(session, "ranges", choices=update.ranges.output(ranges))
  }, label="Set configuration for current LM")
  
  
  getCurrentStateForRestrictions <- reactive({ # depends on input$plot_dblclick
    click <- input$plot_dblclick
    result <- list(click=isolate(click), state=-1)
    
    if(!is.null(click)){
      result$state <- if(is.null(beginRestrictionInterval)) 1 else 2
    }
    
    return(result)
  }, label="get restrictions state")
  
  observe({ # depends on getCurrentStateForRestrictions(), restrictionState${ranges, state}, getCurrentLM()
    whichLM <- isolate(getCurrentLM())
    state <- getCurrentStateForRestrictions()
    click <- state$click
    
    if(state$state == 1){ # (1)=>(2)
      beginRestrictionInterval <<- list(x=click$x, y=click$y)
      restrictionState$state <- 2
    }else if((state$state == 2) && (!is.null(beginRestrictionInterval$x))){ # (2)=> (1)
      # If we click close to the first point, the first point is removed and we are in the state 1 again...
      close.to.frist.Pts <- nearPoints(as.data.frame(beginRestrictionInterval), click, xvar="x", yvar="y")
      if(dim(close.to.frist.Pts)[1] == 0){ # Otherwise: we add the the point to our ranges
        ranges <- isolate(configuration$ranges[[whichLM]])
        ranges <- updateRanges(ranges, c(beginRestrictionInterval$x, click$x))
        configuration$ranges[[whichLM]] <- ranges
      } 
      restrictionState$state <- 1
      beginRestrictionInterval <<- NULL
    }
  }, label="restrictions state changer")
  
  
  prev.deleted <<- "" # Somehow the deleteRange is multiple times not null. So we need to store a value to distinguish them.
  observe({ # depends on input${deleteRange. ranges, whichLM}. Sets configuration$ranges
    if(!is.null(input$deleteRange)){
      whichLM <- isolate(getCurrentLM())
      res <- delete.ranges(whichLM, isolate(configuration$ranges), isolate(input$ranges))
      if (!is.null(res))
        configuration$ranges[[whichLM]] <- res
    }
  }, label="delete ranges")
  
  
  observe({
    if(input$addToHistory){
      addToHistory(isolate(getCurrentFamilyName()), isolate(getCurrentGridName()), 
                   isolate(configuration$df), isolate(configuration$ranges),
                   isolate(configuration$useExisting))
    }
  }, label="save local grid")
  
  observe({
    if(input$loadFromHistory){
      ##1
      values <- loadFromHistory(getCurrentFamilyName(), getCurrentGridName())
      configuration$df <- values$df
      configuration$ranges <- values$ranges
      configuration$useExisting <- values$useExisting
    }
  }, label="load from grid")
######################################################################################
  # zoom
  zoom <- reactiveValues(xlim=NULL, ylim=NULL)
  reset.zoom <- function(){
    zoomHistory <<- list() # HACK: only global variables can be used to pass states between reactives
    zoom$xlim <- NULL
    zoom$ylim <- NULL
  }
    
  observe({ # depends on input${plot_brush}, zoom, 
    res <- zoomIn(input$plot_brush, isolate(reactiveValuesToList(zoom)), zoomHistory)
    if(!is.null(res)){
      zoom$xlim <- res$xlim
      zoom$ylim <- res$ylim
      
      updateNumericInput(session, "zoomYlimMin", value=res$ylim[1])
      updateNumericInput(session, "zoomYlimMax", value=res$ylim[2])
    }
  }, label="Zoom in")
  
  # zoom out
  observe({ # depends on input${zoomOut}, modifies zoom
    if (input$zoomOut){
      res <- zoomOut()
      if(!is.null(res)){
        zoom$xlim <- res$xlim
        zoom$ylim <- res$ylim
        
        updateNumericInput(session, "zoomYlimMin", value=res$ylim[1])
        updateNumericInput(session, "zoomYlimMax", value=res$ylim[2])
        # The event for replotting should be fired now
      }
    }
  }, label="Zoom Out")
  
  # zoom by numeric input fields
  observe({
#     if(input$zoomYlimMin){
#       zoom$ylim[1] <- isolate(input$zoomYlimMin)
#     }
#     if(input$zoomYlimMax){
#       zoom$ylim[2] <- isolate(input$zoomYlimMax)
#     }
  }, label="zoom by numeric input fields")
  
######################################################################################
  ## Reset function
  observe({ # depends on zoom, input${whichLM, familyName, gridName}
    reset <- function(){ 
      reset.zoom(); 
      resetConfiguration()
    }
    
    if(!is.null(getCurrentFamilyName())) reset()
    if(!is.null(getCurrentGridName())) reset()
  }, label="reset")
  
######################################################################################
  ## Save to grid
  observe({ # Depends: input${saveGrid, familyName, gridName, editingGrid}, configuration$df, getEditingGrid()
    
    if(input$saveGrid){
      saveGridToCsv(familyName=isolate(getCurrentFamilyName()),
                       gridName=isolate(getCurrentGridName()),
                       editingGrid=isolate(getEditingGrid()),
                       origSmoothedGrid=isolate(getOriginalSmoothedGrid()),
                       useExisting=isolate(configuration$useExisting),
                       dfs=isolate(configuration$df),
                       ranges=isolate(configuration$ranges))
      
      
      ####################################################
      # TEST of saveGridToCsv
      ####################################################
      if(TEST.save.grid){
        whichLM <- getCurrentLM()
        grid <- isolate(getPostSmoothedEditingGrid())
        if(!is.null(smoothed.totalgrid) && whichLM >= 1 && whichLM <= (ncol(smoothed.totalgrid)+1)){
          
          # Test dims
          if (!all(dim(grid), dim(smoothed.totalgrid))){
            print("not equal sizes: dim(grid)=", paste(dim(grid), collapse=","), 
                                  "dim(st.grid)=", paste(dim(smoothed.totalgrid), collapse=","))
          }else{
            # print("dims OK")
          }
          
          # Test first column
          diff.param <- abs(grid[,1] - smoothed.totalgrid[,1])
          if (! all(diff.param <= 10^-7)){
            print("invalid params:")
            print(paste("grid", paste(head(grid[,1]), collapse=",")))
            print(paste("st.grid", paste(head(smoothed.totalgrid[,1]), collapse=",")))
          }else{
            # print("params OK")
          }
          
          # test cur lm column
          whichLM <- isolate(getCurrentLM())+1
          diff.lm.grid <- abs(grid[,whichLM] - smoothed.totalgrid[,whichLM])
          if (! all(diff.lm.grid <= 10^-7)){
            print(paste("invalid lm grid:", whichLM))
            print(paste("grid", paste(head(grid[,whichLM]), collapse=",")))
            print(paste("st.grid", paste(head(smoothed.totalgrid[,whichLM]), collapse=",")))
          }else{
            # print(paste("grid", whichLM, "OK"))
          }
          
        }
      }
      ####################################################
    }
  }, label="save grids to csv")
  
  ######################################################################################
  ## plot
  output$out <- renderPlot({
    args <- list(grid          = getPostSmoothedEditingGrid(),
                 grid.orig     = getOriginalGrid(),
                 grid.smoothed = getOriginalSmoothedGrid(),
                 idxCol        = getCurrentLM() + 1, # +1 because 1. col are the points of xi
                 xlab          = expression(xi),
                 lty           = c(2, 3, 1),
                 lwd           = c(0.8, 0.8, 1.8),
                 col           = 1:3,
                 main          = getMainTitle(getCurrentGridName(), getCurrentFamilyName()),
                 ylab          = getLMName(getCurrentLM()),
                 restriction   = plotGridRestriction()[[1]],
                 withLegend    = input$withLegend)
    # Zoom
    args[["xlim"]] <- zoom$xlim
    args[["ylim"]] <- zoom$ylim
    
    ## plot
    do.call(matlines.internal, args)
    ## smooth restricton selector
    # (1) => plot first.
    # (2) => plot first.
    # (3) => nothing.
    
    if(!is.null(restrictionState$state)) {
      state <- isolate(restrictionState$state)
      if(state == 2){
        click <- beginRestrictionInterval
        points(click$x, click$y, pch=4, col="red", lwd=3)
      }
    }
  })
})
