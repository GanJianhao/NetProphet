# 
# library("shiny")
# library("igraph")
# #library("ENA")
# library("shinyBS")
# library("knitr")

`%then%` <- shiny:::`%OR%`
reactiveAdjacencyMatrix <- function(func){
    reactive({
    
    val <- func()
    
    if (is.null(val)){
      return(list(names=character(), links=list(source=-1, target=-1)))
    }
    
    #TODO: re-arrange columns if necessary
    if (!all(rownames(val) == colnames(val))){
      stop("Colnames and rownames of your matrix must be identical")
    }
    
    diag(val) <- 0
    
    #make the matrix symmetric
    val <- symmetricize(val, method="avg")
    
    #now consider only the upper half of the matrix
    val[lower.tri(val)] <- 0
    
    conns <- cbind(source=row(val)[val>0]-1, target=col(val)[val>0]-1, weight=val[val>0])
    
    if (nrow(conns) == 0){
      conns <- list (source=-1, target=-1, weight=0)
    }
    
    list(names=rownames(val), links=conns)
    
  })
  
}

check_files_uploaded <- function(inp) {
  if (is.null(inp$data) |  is.null(inp$pert) | is.null(inp$DE) | is.null(inp$tf_orfs) | is.null(inp$orfs)) {
    "One or more file(s) not uploaded. Please upload your data files"
  }else {
    NULL
  }
}

read_files_to_matrix <-function(input){
  regFile <- try(read.table(toString(input$regulatorGeneNamesFileName)))
  if(class(regFile)=='try-error'){ return(NULL)}
  regulators <- as.matrix(regFile)
  tarFile <- try(read.table(toString(input$targetGeneNamesFileName)))
  if(class(tarFile)=='try-error'){return(NULL)}
  targets <- as.matrix(tarFile)
  pert_file <- toString(input$perturbationMatrixFile)
  pertFile <- try(readLines(pert_file))
  if(class(pertFile)=='try-error') {return(NULL)}
  pert_input <- as.matrix(pertFile) 
  DEFILE <- try(read.table(input$differentialExpressionMatrixFile))
  if(class(DEFILE)=='try-error'){return(NULL)}
  de_component = as.matrix(DEFILE)
  ExpFile <- try(read.table(input$targetExpressionFile))
  if(class(ExpFile)== 'try-error'){return(NULL)}
  #check to see if tdata has header row(condition names)
  twoLines <- readLines(input$targetExpressionFile, n=2)
  firstLine <- strsplit(twoLines, " ")[[1]]
  hdr <- all(is.na(as.numeric(firstLine)))
  if(!hdr) { 
    
    tdata = as.matrix(read.table(input$targetExpressionFile, header=F))
    pert = t(sapply(targets, grepl, pert_input))
    
  }else{ 
    tdata = as.matrix(read.table(input$targetExpressionFile, header=T))
    collist <- firstLine
    rowlist <-as.character(unlist(targets))
    a = sapply(firstLine, grepl, pert_input)
    b = sapply(rowlist, grepl, pert_input)
    d=t(b) %*% a
    colnames(d) <- collist
    pert = d
    
  }
  rdata = tdata[match(regulators, targets),] 
  a = sapply(regulators, grepl, targets)    
  allowed = t(!a) 
  matrixInputFiles =list(tdata = tdata,hdr=hdr, pert_input=pert_input, pert = pert, de_component = de_component, regulatorsNames = regulators, targetsNames = targets,
                         rdata = rdata, allowed = allowed, inputFiles = input)
}

error_checking <- function(input){  
  if(ncol(input$regulatorsNames) != 1) return("Regulator names file should have 1 gene name per each row")
  
  if(ncol(input$targetsNames) != 1) return("Target names file should have 1 gene name per each row")
  #check # of experiments in tdata and DE
  if(nrow(input$tdata) != ncol(input$de_component)){ return("Expression Matrix and Differential Expression Matrix should have the same number of columns")}
  # check pert for file with header
  if(!input$hdr){
    
    #if no hdr, pert should have #cond rows
  if(nrow(input$pert_input) != ncol(input$tdata)){ return("your Expression file has no header. Therefore the perturbation file should have exactly n rows where n is the number of expression conditions. Please refer to the help page for more info")}
  }
  #check if pert has any true value
  if(!any(input$pert)){ return(" perterbation matrix is empty.")}
  
  NULL
  
}

shinyServer(function(input, output) {
    inputData <- reactive({
      validate(
        check_files_uploaded(input)
      )  
      targetExpressionFile = input$data$datapath
      differentialExpressionMatrixFile = input$DE$datapath
      perturbationMatrixFile = input$pert$datapath
      regulatorGeneNamesFileName = input$tf_orfs$datapath
      targetGeneNamesFileName = input$orfs$datapath
      inputFiles <- list(targetExpressionFile= targetExpressionFile, differentialExpressionMatrixFile=differentialExpressionMatrixFile,
                         perturbationMatrixFile=perturbationMatrixFile, regulatorGeneNamesFileName=regulatorGeneNamesFileName,
                         targetGeneNamesFileName=targetGeneNamesFileName)
      
      inputTables <- read_files_to_matrix(inputFiles)
      validate(
        need(!is.null(inputTables), "One of the input files is empty") %then%
        error_checking(inputTables)
      )
      inputTables
    })
    outputData <- reactive({
     
      input$getData
      if(input$getData==0) return(NULL)
  #     tempdir <- paste(tempdir(),"outputs", sep="/")
      results <- list()
      inputTables = isolate(inputData())
      source("./CODE/run_netprophet.r", local = TRUE)
      
      cat("I'm done")
      results

   })
  
  output$table <- renderDataTable({
    outputData()$combinedAdjLst
  }, options = list(aLengthMenu = c(10, 30, 50), iDisplayLength = 10))

  output$download_matrix <- downloadHandler(
    filename = function() { 
      validate(
        need(is.null(outputData())==FALSE, "No data to download yet")
      )
      paste('combined_model_matrix', '.txt', sep='') },
    content = function(file) {

    write.csv(outputData()$combinedAdjMtr, file)
  })
  output$download_list <- downloadHandler(
    filename = function() {
      validate(
        need(is.null(outputData())==FALSE, "No data to download yet")
      )
      paste('combined_model_list', '.txt', sep='') },
    content = function(file) {
      write.csv(outputData()$combinedAdjLst, file)
    })
  output$download_lasso <- downloadHandler(
    filename = function() { 
      validate(
        need(is.null(outputData())==FALSE, "No data to download yet")
      )
      paste('lasso_component', '.txt', sep='') },
    content = function(file) {
      write.csv(outputData()$lasso_component, file)
    })
  
#   output$downloadData <- downloadHandler(
#     filename = function() { 
# #       validate(
# #         need(input$getData>0, "No data to download yet. Please upload your data first and then press 'Calculate the results'"), %then%
# #         need()
# #           )
#       paste("netprophet_Output", '.tar', sep='') },
#     content = function(file) {
#       browser()
#       tar(file, results) 
#     }
#   )
#   
#   output$multiDownload <- downloadHandler(
#     filename = function(){
#       # Time-stamp tar filename
#       paste0("netprophet_output_files-", gsub("\\D", "_", Sys.time()), ".csv")
#     },
#     content = function(file){
#       browser()
#       # Loads the data. Could be reactive.
#       tempData <- inputDataset()
#       if(is.null(tempData)){
#         # Protect against empty data.
#         return(NULL)
#       }
#       
#         write.csv(tempData, file)
#       
# #       plateNames = unique(tempData$PlateID)
# #       tempdir = tempdir()
# #       dlist = plyr::dlply(tempData, "PlateID", function(x) x[, c("names", "well", "comments")])
# #       for(i in plateNames){
# #         write.csv(x = dlist[[i]], file = paste0(tempdir, "/", i, ".csv"), row.names = FALSE)
# #       }
# #       tar(tarfile = file, files = tempdir)
#     }
#   )
  
  output$mainnet <- reactiveAdjacencyMatrix(function() {
    
      
    mydata <- read.delim(".\\DATA\\YEAST_SUBNETWORK\\OUTPUT\\GLOBAL_SHRINKAGE\\combined_model.adjmtr", header=F, sep="\t")  
    gD <- simplify(graph.data.frame(mydata, directed=T))
    gAdj <- get.adjacency(gD,type="both", edges = T, sparse = FALSE)
    adjMatrix <- matrix(0, ncol=200, nrow=(200-24))
    adj <-rbind(adjMatrix, mydata)
    data <- read.delim(".\\DATA\\YEAST_SUBNETWORK\\INPUT\\data.withlbls.expr", header = T, sep = " ", quote = "\"'",dec = ".", row.names=1)
    names(adj) <-row.names(data)
    row.names(adj) <- row.names(data)
    net <- adj
    net <- abs(net)
    net[net < input$con_weight] <- 0
    net
    
   
    
  })
})
