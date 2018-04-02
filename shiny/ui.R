library(shiny)

shinyUI(fluidPage(
  
    sidebarLayout(
    position="right",
    sidebarPanel(
      helpText("Optionally, enter a word for highlighting, or two words to highlight shortest paths between them:"),
      textInput("word1", 
                   label = "First word:", 
                   value = ""),
      textInput("word2", 
                label = "Second word:", 
                value = "")
    ),
    
    mainPanel(  
      helpText("The main plot..."),
      checkboxInput("giant", label = "Restrict to giant component?", value = FALSE),
      helpText("Select and drag to zoom in:"),
      plotOutput("mainplot", height=600,
                  brush = brushOpts(id = "plot_brush")
       )
    )
  ),
  sidebarLayout(
    position="right",
    sidebarPanel(
      helpText(""),
      plotOutput("plot_hover")
      #, verbatimTextOutput("test")
    ),
    mainPanel(  
      plotOutput("subplot", height=600,
                 hover = hoverOpts(id="plot_hover", delayType="throttle")
      ),
      wellPanel(selectInput("layout", 
                            label = "Zoom graph layout:",
                            choices = c("As in main plot above",
                                        "Fruchterman-Reingold",
                                        "Kamada-Kawai",
                                        "Reingold-Tilford"),
                            selected = "Fruchterman-Reingold")
      )
    )
  )
  )
  )