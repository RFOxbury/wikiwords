library(shiny)
source("helpers.R")

shinyServer(
  
  function(input, output, session) {

    vset1 <- reactive({
      w <- input$word1
      if(w=="") return()
      findnodes(w)
    })
    vset2 <- reactive({
      w <- input$word2
      if(w=="") return()
      findnodes(w)
    })
    subv <- reactive({
      res <- brushedPoints(lout, input$plot_brush, "X", "Y")
      if (nrow(res) == 0){
        return()
      } else{
        as.numeric(row.names(res))
      }
    })
    currentgraph <- reactive({
      induced.subgraph(mg, subv())
    })
    subsetdf <- reactive({
      g <- currentgraph()
      df <- switch(input$layout,
               "As in main plot above" = lout[subv(),],
               "Fruchterman-Reingold" = data.frame(layout.fruchterman.reingold(g)),
               "Kamada-Kawai" = data.frame(layout.kamada.kawai(g)),
               "Reingold-Tilford" = data.frame(layout.reingold.tilford(g))  )
      row.names(df) <- subv()
      names(df) <- c("X","Y")
      df
    })
    output$mainplot <- renderPlot({
      if(input$giant) handplot(giant, lout[giantv,])
      else handplot(mg, lout)
      # draw paths:
      overlay_paths(vset1(), vset2(), lout, cex=0.7)
    })
    output$subplot <- renderPlot({
      if(is.null(input$plot_brush)) return()
      g <- currentgraph()
      handsubplot(g, subsetdf())
      # draw paths:
      overlay_paths(vset1(), vset2(), subsetdf(), cex=4)
    })
    pts <- reactive({
      if(is.null(input$plot_brush)) return()
      if(is.null(input$plot_hover)) return()
      res <- nearPoints(subsetdf(), input$plot_hover, 
                        "X", "Y", threshold=10, maxpoints=1)
      if (nrow(res) == 0)
        return()
      idx <- as.numeric(row.names(res))
      cluster.set[[idx]]$cluster
    })
    output$plot_hover <- renderPlot({
      if(is.null(pts())) return()
      wordcloud(word[pts()], 
                freq[pts()], 
                random.order=TRUE, 
                scale=c(1.8,1.0), rot.per=0
                )
    })
    #output$test <- renderPrint({
    #   cat(str(input$plot_hover))
    #})
  }
)

