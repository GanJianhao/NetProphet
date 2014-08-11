# library("shinyBS")
# library("knitr")
reactiveNetwork <- function (outputId) 
{
  HTML(paste("<div id=\"", outputId, "\" class=\"shiny-network-output\"><svg /></div>", sep=""))
}


googleAnalytics <- function(account="UA-36850640-1"){
  HTML(paste("<script type=\"text/javascript\">

    var _gaq = _gaq || [];
  _gaq.push(['_setAccount', '",account,"']);
  _gaq.push(['_setDomainName', 'rstudio.com']);
  _gaq.push(['_trackPageview']);

  (function() {
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  })();

  </script>", sep=""))
}

# Code to make a message that shiny is loading
# Make the loading bar
loadingBar <- tags$div(class="progress progress-striped active",
                       tags$div(class="bar", style="width: 100%;"))
# Code for loading message
loadingMsg <- tags$div(class="modal", tabindex="-1", role="dialog", 
                       "aria-labelledby"="myModalLabel", "aria-hidden"="true",
                       tags$div(class="modal-header",
                                tags$h3(id="myModalHeader", "Loading...")),
                       tags$div(class="modal-footer",
                                loadingBar))
# The conditional panel to show when shiny is busy
loadingPanel <- conditionalPanel(paste("input.getData > 0 &&", 
                                       "$('html').hasClass('shiny-busy')"),
                                 loadingMsg)

shinyUI( fluidPage(theme = "Cerulean.css",
  pageWithSidebar(
  headerPanel(""),
  sidebarPanel(
    
    h3("Upload Data Files"),
    
    tags$p("NetProphet needs 5 files as input. Please upload them properly"),
    fileInput('data', "Expression Data File",
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    bsPopover("data", 'Expression Data', " A space separated matrix of expression values of size # of genes x # of expression conditions.", trigger="hover", placement="right"),
        ######
    fileInput('DE', 'Differential Expression File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    bsPopover("DE", 'Differential Expression Matrix', "A space separated adjacency matrix of size # of regulator genes (rows) x # of target genes (columns).", trigger="hover", placement="right"),
    ######
    fileInput('tf_orfs', 'Regulator Genes Names File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    bsPopover("tf_orfs", 'Regulator Genes', "A file listing one regulator gene identifier per line.", trigger="hover", placement="right"),
    ######
    fileInput('orfs', 'Target Gene Names File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    bsPopover("orfs", 'Target Genes', "A file listing one target gene identifier per line.", trigger="hover", placement="right"),
    ######
    fileInput('pert','perturbation matrix File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    bsPopover("pert", 'Perturbation Matrix', "Please visit the 'help' tab for the format specifications", trigger="hover", placement="right"),
    
    actionButton('getData', 'calculate the result')
    
  ),
  mainPanel(
    
    tags$head(
      tags$style(type="text/css", ".jslider { max-width: 400px; }")
    ),
    tabsetPanel(id="tabSelected",
          tabPanel("Table", loadingPanel, dataTableOutput("table")),
          tabPanel("Download",
                   p("Below you can download the output files of NetProphet. Please refer to the help page for more details about each file"),
                   br(),
                   br(),
                   p("Download the results(combined lasso and differential expression) in the form of adjacency matrix:"),
                   downloadButton('download_matrix', 'Download adjacency matrix'),
                   br(),
                   br(),
                   
                   p("Download the results(combined lasso and differential expression) in the form of adjacency list:"),
                   downloadButton("download_list", 'Download adjacency list'),
                   br(),
                   br(),
                   p("download the Lasso adjacency matrix:"),
                   downloadButton('download_lasso', 'Download Lasso component matrix')),
          
        #  tabPanel("Plot",   
        #           sliderInput(inputId = "con_weight",
        #                       label = "Connection threshold:",
        #                       min = 0.0, max = 1.0, value = .06, step = 0.03),HTML("<hr />"), includeHTML("graph.js"),
        #                       reactiveNetwork(outputId = "mainnet") ),
         tabPanel("Help", includeRmd("help.Rmd"))
                   
          
    )  
  )
)))
