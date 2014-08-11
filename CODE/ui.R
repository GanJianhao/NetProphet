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
loadingPanel <- conditionalPanel(paste("input.get > 0 &&", 
                                       "$('html').hasClass('shiny-busy')"),
                                 loadingMsg)

shinyUI(pageWithSidebar(
  headerPanel("NetProphet"),
  sidebarPanel(
   
    h3("Upload Data Files"),
    tags$p("NetProphet needs 7 files as input. Please upload them properly"),
    fileInput('data', 'Choose target expression File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
#     fileInput('rdata', 'Choose regulator expression File',
#                accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
#     fileInput('allowed', 'Choose allowed matrix File',
#               accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    fileInput('pert', 'Choose perturbation matrix File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    fileInput('DE', 'Choose differential expression File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    fileInput('tf_orfs', 'Choose regulator genes names File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
    fileInput('orfs', 'Choose target gene name File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
     actionButton('get', 'calculate the result')
    
  ),
  mainPanel(
  #  tableOutput('contents')
    tabsetPanel(id="tabSelected",
          tabPanel("Table", loadingPanel, dataTableOutput("table")),
          tabPanel("Download", downloadButton('downloadData', 'Download')),
          tabPanel("Help", uiOutput("help")),
          tabPanel("Plot", plotOutput("basePlot"))
          
    )  
  )
))