setwd("C:/Users/Elaheh/Downloads/MichaelBrentLab/netprophet_0.1")
library(shiny)

shinyServer(function(input, output) {
  
    inputDataset <- reactive({
      
     
 
      targetExpressionFile <- isolate({input$data$datapath})
     differentialExpressionMatrixFile <- isolate(input$DE$datapath)
     perturbationMatrixFile <- isolate(input$pert$datapath)
     regulatorGeneNamesFileName <- isolate(input$tf_orfs$datapath)
     targetGeneNamesFileName <- isolate(input$orfs$datapath)
   
    
    
    if (input$get==0)
      return(NULL)
   

    source("./CODE/run_netprophet.r", local = TRUE)
    
    cat("I'm done")
#     as.table(interactions)
   read.table("./DATA/YEAST_SUBNETWORK/OUTPUT/GLOBAL_SHRINKAGE/combined_model.adjlst", header=TRUE,sep='\t')
   
   })
  
  output$table <- renderDataTable({
      inputDataset()
  
  }, options = list(aLengthMenu = c(5, 30, 50), iDisplayLength = 10))

  output$downloadData <- downloadHandler(
    filename = function() { paste('combined_model', '.csv', sep='') },
    content = function(file) {
    write.csv(inputDataset(), file)
  }
)
  
})