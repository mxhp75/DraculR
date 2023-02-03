# library(BiocManager)
# options(repos = BiocManager::repositories())

library(dplyr)
library(plyr)
library(ggplot2)
library(patchwork)
library(shiny)
library(ggrepel)
library(scales)
library(tidyr)
library(magrittr)
library(reshape)
library(edgeR)
options(repos = BiocManager::repositories())
library(readr)
library(DT)
library(psych)
library(tools)
library(shinyhelper)

##### User defined functions #####

# negate %in%
`%notin%` <- Negate(`%in%`)

# function (x) to count the number of non-zero records in each column (ie per sample)
nonzero <- function(x) sum(x != 0)

##################################

# import our distribution difference dataframe
plotData_distDiff_dCq <- read_csv("www/plotData_distDiff_dCq.csv")

# import the rank and distribtuion difference data for GSE153813
rank_GSE153813 <- read_csv(file = "www/rank_GSE153813.csv")
unlistDistributionDifference_GSE153813 <- read_csv(file = "www/unlist_distributionDifference_GSE153813.csv")
GSE153813_info <- paste("The data containined in GSE153813 were obtained from",
                        " NCBI GEO. There is no associated publication" , sep = "")

countsRaw_GSE153813 <- read_csv(file = "www/counts_raw_GSE153813.csv")

rankDist_GSE153813 <- dplyr::full_join(rank_GSE153813, unlistDistributionDifference_GSE153813, by = "samplename") %>% 
  dplyr::mutate(., project = rep("GSE153813", nrow(.)))

# import the rank and distribtuion difference data for GSE118038
rank_GSE118038 <- read_csv(file = "www/rank_GSE118038.csv")
unlistDistributionDifference_GSE118038 <- read_csv(file = "www/unlist_distributionDifference_GSE118038.csv")
GSE118038_info <- paste("The data contained in GSE118038 were obtained from",
                        " NCBI GEO and originate from the article \"A preliminary study of micro-RNAs",
                        " as minimally invasive biomarkers for the diagnosis of prostate cancer patients.\"",
                        " J Exp Clin Cancer Res 2021 Feb 23;40(1):79.", sep = "")

countsRaw_GSE118038 <- read_csv(file = "www/counts_raw_GSE118038.csv")

rankDist_GSE118038 <- dplyr::full_join(rank_GSE118038, unlistDistributionDifference_GSE118038, by = "samplename") %>% 
  dplyr::mutate(., project = rep("GSE118038", nrow(.)))

# import the rank and distribtuion difference data for GSE105052
rank_GSE105052 <- read_csv(file = "www/rank_GSE105052.csv")
unlistDistributionDifference_GSE105052 <- read_csv(file = "www/unlist_distributionDifference_GSE105052.csv")
GSE105052_info <- paste("The data containined in GSE105052 were obtained from",
                        " NCBI GEO and originate from the article \"Small ",
                        " RNA-seq analysis of circulating miRNAs to identify phenotypic variability in",
                        " Friedreich's ataxia patients.\" Sci Data 2018 Mar 6;5:180021.", sep = "")

countsRaw_GSE105052 <- read_csv(file = "www/counts_raw_GSE105052.csv")

rankDist_GSE105052 <- dplyr::full_join(rank_GSE105052, unlistDistributionDifference_GSE105052, by = "samplename") %>% 
  dplyr::mutate(., project = rep("GSE105052", nrow(.)))

# import the rank and distribtuion difference data for GSE151341
rank_GSE151341 <- read_csv(file = "www/rank_GSE151341.csv")
unlistDistributionDifference_GSE151341 <- read_csv(file = "www/unlist_distributionDifference_GSE151341.csv")
GSE151341_info <- paste("The data containined in c were obtained from",
                        " NCBI GEO and originate from the article \"Sequencing",
                        " identifies a distinct signature of circulating microRNAs",
                        " in early radiographic knee osteoarthritis.\"",
                        " Osteoarthritis Cartilage 2020 Nov;28(11):1471-1481.", sep = "")

countsRaw_GSE151341 <- read_csv(file = "www/counts_raw_GSE151341.csv")

rankDist_GSE151341 <- dplyr::full_join(rank_GSE151341, unlistDistributionDifference_GSE151341, by = "samplename") %>% 
  dplyr::mutate(., project = rep("GSE151341", nrow(.)))

# global objects for imported data calculations

classifier_miRs <- data.frame(
  SYMBOL = c(
    "hsa-miR-106b-3p",
    "hsa-miR-140-3p",
    "hsa-miR-142-5p",
    "hsa-miR-532-5p",
    "hsa-miR-17-5p",
    "hsa-miR-19b-3p",
    "hsa-miR-30c-5p",
    "hsa-miR-324-5p",
    "hsa-miR-192-5p",
    "hsa-miR-660-5p",
    "hsa-miR-186-5p",
    "hsa-miR-425-5p",
    "hsa-miR-25-3p",
    "hsa-miR-363-3p",
    "hsa-miR-183-5p",
    "hsa-miR-451a",
    "hsa-miR-182-5p",
    "hsa-miR-191-5p",
    "hsa-miR-194-5p",
    "hsa-miR-20b-5p"
  ))

ui <- fluidPage(navbarPage(title = "DraculR",
                           
                           tabPanel("Methods",
                                    
                                    tags$h4("Background"),
                                    tags$br(),
                                    
                                    fluidRow(
                                      
                                      column(8,
                                             tags$h5("Welcome to DraculR, a Shiny App designed to help you detect red blood cell content contamination in miR-Seq data from human plasma."),
                                             tags$h5(HTML(paste(
                                               "All code used to calculate data shown here is available at the following",
                                               tags$a(href="https://github.com/mxhp75/DraculR.git", "git repository")
                                             )))),
                                      column(8,
                                             tags$h5(HTML(paste(
                                               "If you have questions",
                                               tags$a(href="mailto::melanie.smith@flinders.edu.au", "email me")
                                               )))),
                                      
                                      column(8,
                                             tags$h5(HTML(paste(
                                               "Note: If you plan to use DraculR regularly, we suggest downloading a local copy rather than using a web version. Instructions on how to do this can be found on the github repository, or at minute 18:39 here",
                                               tags$a(href="https://www.youtube.com/watch?v=vX3krP6JmOY&t=1119s", "NetworkChuck"),
                                               "or here",
                                               tags$a(href="https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository", "GitHub")
                                             )))),
                                      
                                      column(8,
                                             tags$h5(HTML(paste(
                                               "The Haemolysis metric described",
                                               tags$a(href="https://doi.org/10.3390/genes13071288", "here"),
                                               "and implemented in the DraculR ShinyR web-based application, performs an",
                                               tags$i("in silico"),
                                               "quality assessment to detect evidence of haemolysis contamination in the original plasma specimen, assigning each small RNA sequencing dataset into one of two categories."
                                             )))),
                                      
                                      column(8,
                                             tags$h5("The classification of ‘Clear’ or ‘Caution’ is designed to alert the user to potential quality control issues in the original plasma specimen."))
                                      
                                      ),
                                    
                                    tags$br(),
                                    tags$h4("Application"),
                                    
                                    fluidRow(

                                      column(8,
                                             tags$h5(HTML(paste(
                                               "DraculR provides a simple graphical user interface (GUI) that allows the user to upload small RNA sequencing data in the form of a raw counts table with names in mature miRNA format (e.g. hsa-miR-1-5p).",
                                               "The",
                                               tags$b("Import New Data"), 
                                               "tab provides the option to personalise table and figure titles, set filtering options and remove miRNA known to be differentially abundent between the user groups of interest",
                                               "from the calculation of the Haemolysis Metric.", "This final step ensures the calculation for haemolysis is not confounded in the event that one or more of the miRNA signature set is known to differ between groups.",
                                               sep = " "
                                               )))),
                                        
                                        column(8,
                                               tags$h5(HTML(paste(
                                                 "In Figure 1 we present an example distribution illustrative of a dataset classified as Clear (Figure 1a) and as one as Caution (Figure 1b).",
                                                 "In the first example (Figure 1a), the distance between the geometric mean of the background miRNA (blue) compared to that of the signature set miRNA (red) is very small whereas in the second example (Figure 1b) the distance between the geometric mean of the background miRNA compared to that of the signature set miRNA is larger.",
                                                 "Furthermore, this difference is greater than the Haemolysis Metric threshold of 1.9, established in",
                                                 tags$a(href="http://dx.doi.org/10.3390/genes13071288", "Smith et al."),
                                                 ", and thus the sample has been classified as ‘Caution’.",
                                                 "In Figure 1b the signature set distribution suggests these miRNAs are over-represented relative to Figure 1a",
                                                 "This suggests that RBC-associated miRNA have been added to the pool of miRNA isolated in the plasma.",
                                                 "In the haemolysed example (Figure 1b), we would recommend removing the sample data from further analysis.",
                                                 "However, where a decision is made to retain samples, the issue of haemolysis should be noted and may be a limitation to inference.",
                                                 sep = " "
                                               ))))
                                      ),
                                      
                                    fluidRow(
                                      column(12,
                                             p(
                                               tags$img(
                                                 src = "distribution.png",
                                                 alt = "Side by side distributions",
                                                 width = 800,
                                                 height = 500
                                               ),
                                               ""
                                             )),
                                      column(4,
                                             p(
                                               "Figure 1: In the first example (a) the sample is classified as ‘Clear’ indicating no evidence for haemolysis. The distance between the geometric mean of background and signature set miRNA is small. In the second example (b) the sample is classified as ‘Caution’ indicating that we found evidence suggestive of haemolysis. The geometric mean of background and signature set miRNA is further apart than that expected where no haemolysis is present."
                                               ))
                                      ),
                                    
                                    tags$br(),
                                    tags$h4("Data input"),
                                    
                                    fluidRow(
                                      column(8,
                                             tags$h5(HTML(paste(
                                               "DraculR allows the user to upload a raw, high throughput sequencing counts table for analysis.",
                                               "Normalisation is performed using the Trimmed Mean of M method (TMM) previously recommended in the",
                                               tags$a(href="http://dx.doi.org/10.1093/bioinformatics/btp616", "edgeR"),
                                               "workflow and the user controls features such as filtering for low expression (Step 1) and refining the haemolysis signature set based on",
                                               tags$i("a priori"),
                                               "knowledge of miRNAs that may be differentially expressed in the comparison of interest (Step 2).",
                                               "The purpose of removing miRNAs with a known association to the condition of interest is to help ensure any issues with haemolysis are not confounded with the research hypothesis.",
                                               "Note that samples with total miRNA read counts < 1 million are considered to be poorly sequenced and are recommended to be removed for quality control.",
                                               sep = " "
                                             ))))
                                    ),
                                    
                                    tags$br(),

                                    fluidRow(
                                      column(12,
                                             p(
                                               tags$img(
                                                 src = "step_1.png",
                                                 alt = "Step 1. Upload, filter and normalise",
                                                 width = 800,
                                                 height = 320
                                               ),
                                               ""
                                             )),
                                      column(4,
                                             p(
                                               "Step 1: Import a raw counts table generated by high throughput miRNA sequencing of human plasma libraries. These data will be filtered according to user specified requirements (n = number of samples in the smallest group of interest) and normalised using the Trimmed Mean of M (TMM) method."
                                             )),

                                      column(12,
                                             p(
                                               tags$img(
                                                 src = "step_2.png",
                                                 alt = "Step 2. Calculate the distribution difference",
                                                 width = 800,
                                                 height = 320
                                               ),
                                               ""
                                             )),
                                      column(4,
                                             p(
                                               HTML(paste("Step 2: The distribution difference between the background and signature miRNA counts is calculated on an individual sample basis allowing the user to upload one to many samples as required. In the case of",
                                                     tags$i("a priori"),
                                                     "knowledge of miRNA differentially abundant between a tested condition/control paradigm the user may choose to reduce the signature miRNA such that they do not include miRNA of interest (recommended).",
                                                     sep = " "
                                             ))))

                                      ),
                                    
                                    tags$br(),
                                    tags$h4("Visualisation"),
                                    
                                    fluidRow(
                                      
                                      column(8,
                                             tags$h5(
                                               "An essential feature of DraculR is that it allows users to visualise and assess the values obtained in the results, through sample specific and consolidated graphics including density plots, histograms and tables (Step 3). These features help the user decide on the level of haemolysis that may affect their analysis by providing a new quality metric. Using this metric the user may choose to remove samples from downstream analyses. However, irrespective of whether samples with a Haemolysis Metric above the suggested threshold are removed or retained, the new information may be important to the analysis of their miRNA sequencing data."
                                             ))
                                      ),
                                      
                                      tags$br(),
                                    
                                    fluidRow(
                                      
                                      column(12,
                                             p(
                                               tags$img(
                                                 src = "step_3.png",
                                                 alt = "Step 3. Visualise and interpret results",
                                                 width = 800,
                                                 height = 320
                                               ),
                                               ""
                                             )),
                                      column(6,
                                             p(
                                               "Step 3: Graphical results in the form of a density plot of individual distributions (i, ii) and a histogram of combined distribution differences (iv) are provided along with a combined table of results (iii). The user is provided with both a metric describing the amount of haemolysis and, if appropriate, a recommendation of caution (iii)."
                                             ))
                                      )
                                    ),

                           tabPanel("Instructions",

                                    tags$h4("Getting started"),
                                    
                                    tags$br(),
                                    
                                    fluidRow(
                                      
                                      column(8,
                                             tags$div(HTML(paste(
                                               "The",
                                               tags$b("Public Data Example"),
                                               "tab allows you to click through the plot output and raw data of four datasets available on NCBI GEO that we have run through our method. Samples we consider",
                                               tags$b("Clear"), "are seen in blue, those we consider should be used with",
                                               tags$b("Caution"), "are seen in red.",
                                               sep = " "
                                             ))))),
                                      
                                      tags$br(),
                                      
                                    fluidRow(
                                      column(8,
                                             tags$img(src = "GSE153813_instructionsTab.png",
                                                      height = 400,
                                                      width = 800)),
                                      
                                      column(4,
                                             tags$img(src = "drac.png",
                                                      height = 300,
                                                      width = 300))
                                      
                                      
                                    ),
                                    
                                    fluidRow(
                                      column(8,
                                             tags$h5(HTML(paste(
                                               "The histogram represents the queried dataset where three samples were flagged to be used with caution (red) and six were classified as clear of haemolysis (blue).",
                                               "As a reference, the barcode-style plots below the histogram display the Haemolysis Metric values of samples from an example dataset validated using the ΔCq method",
                                               tags$a(href="http://dx.doi.org/10.3390/genes13071288", "Smith et al."),
                                                "; red for those identified as haemolysed (ΔCq>7), blue for those classified as 'Clear' of haemolysis (ΔCq<7)."
                                             ))))
                                    ),
                                    
                                    fluidRow(
                                      
                                      column(12,
                                             tags$h5(HTML(paste("To import new data, move to the",
                                                           tags$b("Import new data"), "tab. From here click",
                                                           tags$b("Browse"),
                                                           "to access the miRNA expression value table (either as raw counts or normalised CPMs).")))),
                                      
                                      column(12,
                                             tags$h5("DraculR requires all input count files to be formatted as follows:")
                                      )
                                    ),
                                    
                                    
                                    fluidRow(
                                      column(12,
                                             tags$div(
                                               tags$ul(
                                                 tags$h4(
                                                   tags$li("Comma or tab delimited."),
                                                   tags$li("Samples in columns, miRNA expression observations in rows."),
                                                   tags$li("Column with miRNA names should have the column label “mirna_name”."),
                                                   tags$li("miRNA names should be in the miR-base format, e.g. “hsa-miR-123-3p”."),
                                                   tags$li("Sample names should be the other column labels"),
                                                   tags$li("Sample names should not include white space or special characters")
                                                 )
                                               )
                                             ),
                                             tags$h5(HTML(paste("An example input file can be downloaded",
                                                                tags$a(href="https://github.com/mxhp75/DraculR/tree/main/dataExample", "here"))
                                             )))),
                                    
                                    br(),
                                    
                                    fluidRow(
                                      
                                      column(8, 
                                             tags$h4("Importing a new file will populate the main page with a set of tabs that allow you to navigate through the new information.")
                                             )
                                      ),
                                    
                                    fluidRow(
                                      
                                      tags$img(src = "Picture_2.png",
                                               height = "40%",
                                               width = "40%")
                                    )),
                           
                           tabPanel("Public Data Example",
                                    
                                    tags$br(),
                                    
                                    tags$h5("Click one of these buttons to see the results"),
                                    
                                    fluidRow(
                                      column(12,
                                             actionButton(inputId = "GSE153813", label = "GSE153813"),
                                             actionButton(inputId = "GSE118038", label = "GSE118038"),
                                             actionButton(inputId = "GSE105052", label = "GSE105052"),
                                             actionButton(inputId = "GSE151341", label = "GSE151341"))
                                    ),
                                    
                                    br(),
                                    
                                    fluidRow(
                                      column = 6,
                                      offset = 0,
                                      textOutput("projectInfo")
                                    ),
                                    
                                    br(),
                                    
                                    fluidRow(column(8,
                                                    offset = 0.5,
                                                    plotOutput(outputId = "distributionDifference") %>% 
                                                      helper(type = "inline",
                                                             icon = "exclamation",
                                                             title = "Legend Help",
                                                             colour = "red",
                                                             size = "l",
                                                             content = c("The colours used here are consistent across all tabs",
                                                                         "<b>Caution GSE_____<b> = samples from the public data example with a Haemolysis Metric > 1.9 and should be used with caution",
                                                                         "<b>Clear GSE_____<b> = samples from the public data example with a Haemolysis Metric < 1.9 and are considered clear for use",
                                                                         "<b>Haemolysed (dCq)<b> = samples used in an example validation experiment with the qPCR method. These samples had dCq values > 7",
                                                                         "<b>Clear (dCq)<b> = samples used in an example validation experiment with the qPCR method. These samples had dCq values < 7")))
                                             ),
                                    br(),
                                    
                                    fluidRow(
                                      column(8, offset = 0,
                                             # DT::dataTableOutput("rawCounts"))
                                             DT::dataTableOutput("publicResults"))
                                    )),
                           
                           tabPanel("Import New Data",
                                    
                                    sidebarLayout(
                                      sidebarPanel(
                                        fileInput("rawDataFile","Upload the file"), # fileinput() function is used to get the file upload control option
                                        helpText("Max. file size is 5MB"),
                                        tags$hr(),
                                        fluidRow(
                                          column = 6,
                                          h5(helpText("Add a project title")),
                                          textInput(inputId = "project",
                                                    label = "Project",
                                                    "myProjectName"),
                                          verbatimTextOutput("value"),
                                          column = 6, offset = 6,
                                          
                                          h5(helpText("Apply your filtering value"))  %>% 
                                            
                                            helper(icon = "question",
                                                   colour = "green",
                                                   type = "markdown",
                                                   content = "filtering"),

                                          numericInput("filterNum", label = h5("Number in smallest group"), value = 1)),
 
                                        h5(helpText("Select the input file parameters below")),
                                        
                                        radioButtons(inputId = 'sep',
                                                     label = 'File separator',
                                                     choices = c(Comma = ',',
                                                                 Tab = '\t'),
                                                     selected = ','),
                                        
                                        h5(helpText("Select any miRNA that are differentially expressed between your groups ")) %>% 

                                          helper(icon = "question",
                                                 colour = "green",
                                                 type = "markdown",
                                                 content = "drop"),

                                        checkboxGroupInput(inputId = 'drop_miRs',
                                                     label = "Drop?",
                                                     choices = c("hsa-miR-106b-3p" = 'hsa-miR-106b-3p',
                                                                 "hsa-miR-140-3p" = 'hsa-mir-140-3p',
                                                                 "hsa-miR-186-5p" = 'hsa-miR-186-5p',
                                                                 "hsa-miR-425-5p" = 'hsa-miR-425-5p',
                                                                 "hsa-miR-142-5p" = 'hsa-miR-142-5p',
                                                                 "hsa-miR-532-5p" = 'hsa-miR-532-5p',
                                                                 "hsa-miR-17-5p" = 'hsa-miR-17-5p',
                                                                 "hsa-miR-25-3p" = 'hsa-miR-25-3p',
                                                                 "hsa-miR-363-3p" = 'hsa-miR-363-3p',
                                                                 "hsa-miR-183-5p" = 'hsa-miR-183-5p',
                                                                 "hsa-miR-660-5p" = 'hsa-miR-660-5p',
                                                                 "hsa-miR-451a" = 'hsa-miR-451a',
                                                                 "hsa-miR-19b-3p" = 'hsa-miR-19b-3p',
                                                                 "hsa-miR-182-5p" = 'hsa-miR-182-5p',
                                                                 "hsa-miR-30c-5p" = 'hsa-miR-30c-5p',
                                                                 "hsa-miR-324-5p" = 'hsa-miR-324-5p',
                                                                 "hsa-miR-191-5p" = 'hsa-miR-191-5p',
                                                                 "hsa-miR-192-5p" = 'hsa-miR-192-5p',
                                                                 "hsa-miR-194-5p" = "hsa-miR-194-5p",
                                                                 "hsa-miR-20b-5p" = 'hsa-miR-20b-5p'),
                                                     textOutput("txt")
                                                     )
                                        ),
                                      
                                      mainPanel(
                                        uiOutput("tb")
                                      )
                                    )

)))
  
# Input Validation Functions
is_csv_or_tsv <- function(input_file, radio_separator) {
  file1 <- input_file
  if(radio_separator == ",") { 
    sep_str = "csv"
  } else if(radio_separator == "\t") { 
    sep_str = "tsv"
  } else {
    print("Error: Unable to deduce file type extension")
  }
  
  if (identical(tolower(tools::file_ext(file1$datapath)), sep_str)) {
    NULL
  } else {
    print(paste0("Error: Uploaded file extension doesn't match the chosen delimiter, should be: \"", sep_str, "\""))
  }
}

has_correct_header <- function(df) {
  if("mir_name" %in% colnames(df)) {
    NULL
  } else {
    print("Input file must have a column named \"mir_name\" (containing miR IDs). Could not find a column by this name.")
  }
}

at_least_one_library_one_million_reads <- function(x) {
  
  colSumLibrary <- summarise_all(x, ~if(is.numeric(.)) sum(.) else "Total")

  if(any(colSumLibrary>1000000)) {
    NULL
  } else {
    print("Input file contains no libraries with greater than one million reads.")
  }
}

server <- function(input, output) {
  
  # make sure help functions are active
  observe_helpers(help_dir = "helpfiles",
                  withMathJax = TRUE)
  
  # This reactive function will take the inputs from UI.R and use them for read.table() to read the data from the file. It returns the dataset in the form of a dataframe.
  # file$datapath -> gives the path of the file
  uploadData <- reactive({
    file1 <- input$rawDataFile
    if(is.null(file1)){return()} 
    
    shiny::validate(is_csv_or_tsv(input$rawDataFile, input$sep))
    
    temp_df <- read.table(file = file1$datapath,
                          sep = input$sep,
                          header = TRUE,
                          stringsAsFactors = FALSE)
    
    shiny::validate(has_correct_header(temp_df))
    temp_df

    temp_counts <- read.table(file = file1$datapath,
                              sep = input$sep,
                              header = TRUE,
                              stringsAsFactors = FALSE) %>%
      dplyr::mutate_if(is.integer, as.numeric) %>%
      as.data.frame() %>%
      tibble::column_to_rownames("mir_name") %>%
      replace(is.na(.), 0)

    shiny::validate(at_least_one_library_one_million_reads(temp_counts))
    temp_df
    
  })

  # this reactive output will plot according to the public data reactive buttons
  plotDataPublic_miRNA <- reactiveValues(data = rankDist_GSE153813)
  
  observeEvent(input$GSE153813, { plotDataPublic_miRNA$data <- rankDist_GSE153813 })
  observeEvent(input$GSE118038, { plotDataPublic_miRNA$data <- rankDist_GSE118038 })
  observeEvent(input$GSE105052, { plotDataPublic_miRNA$data <- rankDist_GSE105052 })
  observeEvent(input$GSE151341, { plotDataPublic_miRNA$data <- rankDist_GSE151341 })
  
  output$projectInfo <- renderText({
    paste("You have selected", plotDataPublic_miRNA$data$project[1], ".",get(paste(plotDataPublic_miRNA$data$project[1],"_info", sep = "")))
  })
  
  # this reactive output will display the raw data according to the public data reactive buttons
  rawDataPublic_miRNA <- reactiveValues(data = countsRaw_GSE153813)
  
  observeEvent(input$GSE153813, { rawDataPublic_miRNA$data <- countsRaw_GSE153813 })
  observeEvent(input$GSE118038, { rawDataPublic_miRNA$data <- countsRaw_GSE118038 })
  observeEvent(input$GSE105052, { rawDataPublic_miRNA$data <- countsRaw_GSE105052 })
  observeEvent(input$GSE151341, { rawDataPublic_miRNA$data <- countsRaw_GSE151341 })
  
  # this reactive output will display the results table according to the public data reactive buttons
  tableDataPublic_miRNA <- reactiveValues(data = rankDist_GSE153813)
  
  observeEvent(input$GSE153813, { tableDataPublic_miRNA$data <- rankDist_GSE153813 })
  observeEvent(input$GSE118038, { tableDataPublic_miRNA$data <- rankDist_GSE118038 })
  observeEvent(input$GSE105052, { tableDataPublic_miRNA$data <- rankDist_GSE105052 })
  observeEvent(input$GSE151341, { tableDataPublic_miRNA$data <- rankDist_GSE151341 })
  

  # This reactive function will calculate the distribution difference and be 
  # available in all output tabs
  
  distDiff <- reactive({
    
    counts <- uploadData() %>% 
      dplyr::mutate_if(is.integer, as.numeric) %>%
      as.data.frame() %>% 
      tibble::column_to_rownames("mir_name") %>% 
      replace(is.na(.), 0)
    
    # identify samples with < 1 million reads
    lowCounts <- names(counts[, base::colSums(counts) < 1000000])
    
    # remove columns/samples with readcounts less than 1 million
    counts <- counts[, base::colSums(counts) > 1000000]
    
    # reduce any individual count less than five to zero
    counts[counts < 5] <- 0
    
    # remove miRNAs with zero counts in all samples
    counts <- counts[ base::rowSums(counts)!=0, ]
    
    # create the (super) minimal metadata table
    meta <- dplyr::data_frame(samplename = base::colnames(counts)) %>% 
      dplyr::mutate(., copy = samplename) %>% 
      tidyr::separate(., col = copy, into = c("ID", "condition"), sep = "_")
    
    # rank the samples by read counts and by unique miRs
    # this table will be joined downstream with the distribution difference table
    rank <- base::as.data.frame(base::colSums(counts)) %>%
      magrittr::set_colnames(., "readCounts") %>% 
      dplyr::arrange(., -(readCounts)) %>% 
      tibble::rownames_to_column("samplename") %>% 
      dplyr::left_join(., meta, by = "samplename") %>% 
      dplyr::select(., samplename, readCounts, condition) %>% 
      dplyr::mutate(., rank_readCounts = 1:nrow(.)) %>% 
      dplyr::full_join(.,
                       as.data.frame(t(numcolwise(nonzero)(as.data.frame(counts)))) %>%
                         tibble::rownames_to_column() %>%
                         magrittr::set_colnames(., c("samplename", "unique_miRs")) %>%
                         arrange(., desc(unique_miRs)) %>%
                         mutate(., rank_unique = 1:nrow(.)),
                       by = "samplename")
    
    
    # establish a DGEList object
    DGEList_public <- DGEList(counts = counts,
                              samples = meta)
    
    is.na(DGEList_public$counts) %>% table()
    
    # calculate normalisation factors and apply to the DGEList object
    DGEList_public <- calcNormFactors(DGEList_public, method = "TMM")
    
    # calculate the CPMs
    rawCPM <- cpm(DGEList_public, log = FALSE)
    # remove low expressed genes
    keep.exprs <- rowSums(rawCPM > 40) >= input$filterNum
    DGEList_public <- DGEList_public[keep.exprs,, keep.lib.sizes = FALSE]
    
    ## calculate the difference between the geometric mean of the distributions
    ### here we calculate the geometric mean of the classifier distribution and the
    ### "other" ensuring those taken from the classifier list are not included in
    ### other.
    # create a vector of sample names for use in the lapply
    varc <- dplyr::select(DGEList_public$samples, samplename) %>%
      tibble::remove_rownames() %>% 
      dplyr::pull(., samplename)
    
    # define the dropped classifiers as input from the groupCheckboxInput
    dropped <- subset(classifier_miRs, SYMBOL %in% input$drop_miRs)
    
    # define the final set of classifiers
    final_classifiers <- subset(classifier_miRs, SYMBOL %notin% input$drop_miRs)
    
    
    distributionDifference <- lapply(varc,function(x){
      # calculate the geometric mean of the two distribut ions (1 = classifier, 0 = other, 2 = dropped)
      dtmp <- dplyr::select(as.data.frame(edgeR::cpm(DGEList_public$counts, log = TRUE)), x) %>%
        tibble::rownames_to_column("mirna") %>% 
        mutate(., classifier = as.factor(ifelse(mirna %in% final_classifiers$SYMBOL, 1,
                                                ifelse(mirna %in% dropped$SYMBOL, 2,
                                                       ifelse(mirna %notin% classifier_miRs$SYMBOL, 0, NA)))))
      cdat_tmp <- with(dtmp,tapply(get(x),classifier,geometric.mean,na.rm=T))
      cdat <- data.frame("classifier"=rownames(cdat_tmp),"geometric.mean"=cdat_tmp)
      # calculate the difference between the two geometric means (classifier-other)  
      cdat_out <- dplyr::filter(cdat, classifier == 1)$geometric.mean - dplyr::filter(cdat, classifier == 0)$geometric.mean
      return(cdat_out)
    })
    
    names(distributionDifference) <- varc
    
    unlist_distributionDifference <- do.call(cbind.data.frame, distributionDifference) %>% 
      t() %>%
      magrittr::set_colnames("distributionDifference") %>%
      base::as.data.frame() %>% 
      tibble::rownames_to_column("samplename") %>% 
      dplyr::mutate(., haemoResult = ifelse(distributionDifference < 1.9, "Clear",
                                            ifelse(distributionDifference >= 1.9, "Caution", NA)))
    
    unlist_distributionDifference$haemoResult <- as.factor(unlist_distributionDifference$haemoResult)
    
    caution <- dim(filter(unlist_distributionDifference, haemoResult == "Caution"))
    
    dplyr::full_join(rank, unlist_distributionDifference, by = "samplename") %>%
      dplyr::mutate(., project = rep(input$project, nrow(.)))

  })
  
  
  # output$rawCounts <- DT::renderDataTable({
  #   
  #   rawDataPublic_miRNA$data
  #   
  # }, options = list(autoWidth = TRUE,
  #                   columnDefs = list(list(list(targets='_all',
  #                                               visible=TRUE,
  #                                               width='90') ))))
  
  output$publicResults <- DT::renderDataTable({
    
    tableDataPublic_miRNA$data %>% 
      dplyr::select(., Samplename = samplename,
                    `Haemolysis Metric` = distributionDifference,
                    `Haemolysis Result` = haemoResult,
                    Project = project)
    
  }, options = list(autoWidth = TRUE,
                    columnDefs = list(list(list(targets = '_all',
                                                visible = TRUE,
                                                width = '90') ))))
  
  output$distributionDifference <- renderPlot({
    
    # calculate the number of samples with a HM >= 1.9 Caution
    caution <- dim(filter(plotDataPublic_miRNA$data, haemoResult == "Caution"))
  
    # This code renders the barPlot/barcodePlot for public data examples

    # Background Barcode
    ## New facet label names for haemolysis variable
    ## needed to rename the facets
    # barcode.labs <- c("Haemolysed \n(dCq)", "Clear \n(dCq)")
    barcode.labs = c("", "")
    names(barcode.labs) <- c("haemolysed", "none")
    
    # Make the background Barcode plot
    
    geomTextInfo <- data.frame(
      label = c("Haemolysed \n(dCq)", "Clear \n(dCq)"),
      haemolysis = c("haemolysed", "none")
    )
    
    backgroundDataPlot <- ggplot() +
      geom_vline(show.legend = FALSE,
                 data = plotData_distDiff_dCq,
                 aes(xintercept = distributionDifference,
                     color = haemolysis),
                 size = 1) +
      
      coord_cartesian(xlim=c(0,5),
                      clip = "off") +
      
      xlab("Haemolysis Metric") +
      
      scale_colour_manual(values = c("#8B0000", "#29569D"),
                          name = "Haemolysis",
                          # labels = c("Haemolysed \n(dCq)", "Clear \n(dCq)")) +
      labels = c("", "")) +
      
      facet_grid(haemolysis ~ .,
                 # move the labels from right side to left side
                 switch = "y",
                 # swap out the old labels for the new ones
                 labeller = labeller(haemolysis = barcode.labs)) +
      
      theme_minimal() +
     
       # add the geom_text layer to replace the y-axis labels
      geom_text(
        data    = geomTextInfo,
        mapping = aes(x = -.5,
                      y = 0.1,
                      label = label)) +
      
      theme(strip.text.y = element_text(margin = margin(0),
                                        size = 14),
            strip.text.y.left = element_text(angle = 0),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(size = 16),
            axis.text.x = element_text(size = 13), 
            plot.margin = margin(t=0,r=0,b=1,l=2, "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())

    
    # public data example barplot
    newDataPlot <- ggplot() +
      
      geom_histogram(data = plotDataPublic_miRNA$data,
                     aes(x = distributionDifference,
                         fill = haemoResult,
                         colour = haemoResult),
                     breaks = seq(0,5,0.1),
                     alpha = 0.4, 
                     position = "identity",
                     lwd = 0.8) +
      
      coord_cartesian(xlim=c(0,5)) +
      
    scale_fill_manual(values = c("#8B0000", "#29569D",
                                 "#A9A9A9", "#CDCDCD"),
                      name = "Haemolysis",
                      labels = c(paste("Caution ", plotDataPublic_miRNA$data$project[1], sep = ""),
                                 paste("Clear ", plotDataPublic_miRNA$data$project[1], sep = ""),
                                 "Haemolysed (dCq)", "Clear (dCq)")) +
    scale_colour_manual(values = c("#8B0000", "#29569D",
                                   "#8B0000", "#29569D"),
                        name = "Haemolysis",
                        labels = c(paste("Caution ", plotDataPublic_miRNA$data$project[1], sep = ""),
                                   paste("Clear ", plotDataPublic_miRNA$data$project[1], sep = ""),
                                   "Haemolysed (dCq)", "Clear (dCq)")) +
      labs(
        x = "",
        y = "Number of samples",
        subtitle = paste("DraculR identified", caution[1], "samples to use with caution", sep = " ")
      ) +
      ggtitle(paste0(plotDataPublic_miRNA$data$project[1])) +

      # add a vertial dotted line to idicate the HM cut-off
      geom_vline(show.legend = FALSE,
                 xintercept = 1.9,
                 col = 2,
                 lty = 2) +
      
      theme_bw(base_size = 16) +
      theme(plot.subtitle = element_text(color="#8B0000"),
            legend.position = "top")
    
    # combine both plots (this makes sure the x-axix ticks are aligned)
    combinedPlot <- newDataPlot/backgroundDataPlot
    # layout with patchwork
    publicDataBarcodePlot <- combinedPlot +
      patchwork::plot_layout(heights = c(2,1))
    # print the combined plot to the screen
    publicDataBarcodePlot

  })
  
  # this reactive output contains the results table for the public data example
  output$tableDataPublic_miRNA <- renderTable({
    
    tableDataPublic_miRNA$data
    
  })
  
  ## Info for the Import New Data tab
  # this reactive output contains the dataset and display the file information
  output$about <- renderTable({
    if(is.null(uploadData())){return ()}
    
    # full table as uploaded by user
    input$rawDataFile
  })
  
  # this reactive output contains the summary of the dataset and display the summary in table format
  output$sum <- renderTable({
    if(is.null(uploadData())){return ()}
    
    # ordered summary table of results
    # cut the columns we don't need
    distDiff() %>% 
      dplyr::arrange(., samplename) %>% 
      dplyr::select(., Samplename = samplename,
                    `Haemolysis Metric` = distributionDifference,
                    `Haemolysis Result` = haemoResult,
                    Project = project)
    
  })
  
  DGEList_public <- reactive({
    
    counts <- uploadData() %>% 
      dplyr::mutate_if(is.integer, as.numeric) %>%
      as.data.frame() %>% 
      tibble::column_to_rownames("mir_name") %>% 
      replace(is.na(.), 0)
    
    # identify samples with < 1 million reads
    lowCounts <- names(counts[, base::colSums(counts) < 1000000])
    
    # remove columns/samples with readcouns less than 1 million
    counts <- counts[, base::colSums(counts) > 1000000]
    
    # reduce any individual count less than five to zero
    counts[counts < 5] <- 0
    
    # remove miRNAs with zero counts in all samples
    counts <- counts[ base::rowSums(counts)!=0, ]
    
    # create the (super) minimal metadata table
    meta <- dplyr::data_frame(samplename = base::colnames(counts)) %>% 
      dplyr::mutate(., copy = samplename) %>% 
      tidyr::separate(., col = copy, into = c("ID", "condition"), sep = "_")
    
    # rank the samples by read counts and by unique miRs
    # this table will be joined downstream with the distribution difference table
    rank <- base::as.data.frame(base::colSums(counts)) %>%
      magrittr::set_colnames(., "readCounts") %>% 
      dplyr::arrange(., -(readCounts)) %>% 
      tibble::rownames_to_column("samplename") %>% 
      dplyr::left_join(., meta, by = "samplename") %>% 
      dplyr::select(., samplename, readCounts, condition) %>% 
      dplyr::mutate(., rank_readCounts = 1:nrow(.)) %>% 
      dplyr::full_join(.,
                       as.data.frame(t(numcolwise(nonzero)(as.data.frame(counts)))) %>%
                         tibble::rownames_to_column() %>%
                         magrittr::set_colnames(., c("samplename", "unique_miRs")) %>%
                         arrange(., desc(unique_miRs)) %>%
                         mutate(., rank_unique = 1:nrow(.)),
                       by = "samplename")
    
    
    # establish a DGEList object
    DGEList_public <- DGEList(counts = counts,
                              samples = meta)
    
    is.na(DGEList_public$counts) %>% table()
    
    # calculate normalisation factors and apply to the DGEList object
    DGEList_public <- calcNormFactors(DGEList_public, method = "TMM")
    
    # calculate the CPMs
    rawCPM <- cpm(DGEList_public, log = FALSE)
    # remove low expressed genes
    keep.exprs <- rowSums(rawCPM > 40) >= input$filterNum
    DGEList_public <- DGEList_public[keep.exprs,, keep.lib.sizes = FALSE]
    
  })
  
  # This reactive output contains the user uploaded dataset and will display the dataset in table format
  output$table <- renderTable({
    if(is.null(uploadData())){return ()}
    uploadData() %>% 
      dplyr::mutate_if(is.integer, as.numeric) %>%
      replace(is.na(.), 0)
  })
  
  # This reactive output contains the dataset and will render the
  # mature miRNA as a function of read depth scatter plot
  output$distributions <- renderPlot({
    if(is.null(uploadData())){return ()}
    
    # define the final set of classifiers
    final_classifiers <- subset(classifier_miRs, SYMBOL %notin% input$drop_miRs)
    
    # subset the DGEList for one sample
    # temp <- as.data.frame(cpm(DGEList_public()$counts, log = TRUE)) %>%
    #   dplyr::select(., NPC0031_NPC)
    
    selectData <- reactive({
      
      as.data.frame(cpm(DGEList_public()$counts, log = TRUE)) %>%
        dplyr::select(., input$select)
      
    })
    

    # plot side by side densities
    ggplot() +
      geom_density(data = subset(selectData(), rownames(selectData()) %in% final_classifiers$SYMBOL) %>% 
                     dplyr::mutate(., colour = rep("Classifier")),
                   alpha = 0.5,
                   aes(x = !!sym(input$select),
                       fill = colour,
                       colour = colour)) +
      geom_density(data = subset(selectData(), rownames(selectData()) %notin% classifier_miRs) %>% 
                     dplyr::mutate(., colour = rep("Background")),
                   alpha = 0.5,
                   aes(x = !!sym(input$select),
                       fill = colour,
                       colour = colour)) +
      scale_fill_manual(values = c("#29569D", "#8B0000"),
                        name = "miRNA Set",
                        labels = c("Background", "Classifier")) +
      scale_colour_manual(values = c("#29569D", "#8B0000"),
                          name = "miRNA Set",
                          labels = c("Background", "Classifier")) +
      labs(title = paste(colnames(selectData()), " - ",
                         distDiff() %>%
                           dplyr::filter(., samplename == input$select) %>%
                           dplyr::select(., haemoResult) %>%
                           .[[1]],
                         sep = " "),
           x = "log2 CPM",
           y = "Density") +
      theme_bw(base_size = 16)
    
  })
  
  # This reactive output contains new objects and will display the histogram of 
  # distribution difference for uploaded counts above the background barcode
  output$dist_diff <- renderPlot({
    if(is.null(uploadData())){return ()}
    
    # This code renders the barPlot/barcodePlot for the user uploaded data
    
    # Background Barcode
    ## New facet label names for haemolysis variable
    ## needed to rename the facets
    # barcode.labs <- c("Haemolysed \n(dCq)", "Clear \n(dCq)")
    barcode.labs = c("", "")
    names(barcode.labs) <- c("haemolysed", "none")
    
    geomTextInfo <- data.frame(
      label = c("Haemolysed \n(dCq)", "Clear \n(dCq)"),
      haemolysis = c("haemolysed", "none")
    )
    
    # Make the background Barcode plot
    backgroundDataPlot <- ggplot() +
      
      geom_vline(show.legend = FALSE,
                 data = plotData_distDiff_dCq,
                 aes(xintercept = distributionDifference,
                     color = haemolysis),
                 size = 1) +
      
      coord_cartesian(xlim=c(0,5),
                      clip = "off") +
      
      xlab("Haemolysis Metric") +
      
      scale_colour_manual(values = c("#8B0000", "#29569D"),
                          name = "Haemolysis",
                          # labels = c("Haemolysed \n(dCq)", "Clear \n(dCq)")) +
                          labels = c("", "")) +
      
      facet_grid(haemolysis ~ .,
                 # move the labels from right side to left side
                 switch = "y",
                 # swap out the old labels for the new ones
                 labeller = labeller(haemolysis = barcode.labs)) +
      
      theme_minimal() +
      
      # add the geom_text layer to replace the y-axis labels
      geom_text(
        data    = geomTextInfo,
        mapping = aes(x = -.5,
                      y = 0.1,
                      label = label)) +
      
      # tweak labels, margins etc
      theme(strip.text.y = element_text(margin = margin(0),
                                        size = 14),
            strip.text.y.left = element_text(angle = 0),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(size = 16),
            axis.text.x = element_text(size = 13), 
            plot.margin = margin(t=0,r=0,b=1,l=2.5, "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    # make the upload data bar plot
    newDataPlot <- ggplot() +
      
      geom_histogram(data = distDiff(),
                     aes(x = distributionDifference,
                         fill = haemoResult,
                         colour = haemoResult),
                     breaks = seq(0,5,0.1),
                     alpha = 0.4, 
                     position = "identity",
                     lwd = 0.8) +
      
      coord_cartesian(xlim=c(0,5)) +
      
      scale_fill_manual(values = c("#8B0000", "#29569D",
                                   "#A9A9A9", "#CDCDCD"),
                        name = "Haemolysis",
                        labels = c("Caution", "Clear",
                                   "Haemolysed (dCq)", "Clear (dCq)")) +
      scale_colour_manual(values = c("#8B0000", "#29569D",
                                     "#8B0000", "#29569D"),
                          name = "Haemolysis",
                          labels = c("Caution", "Clear",
                                     "Haemolysed (dCq)", "Clear (dCq)")) +
      geom_vline(show.legend = FALSE,
                 xintercept = 1.9,
                 col = 2,
                 lty = 2) +

      labs(
        title = paste0("Haemolysis Metric: ", input$project),
        subtitle = paste("DraculR identified",
                         dim(filter(distDiff(),
                                    haemoResult == "Caution"))[1],
                         "samples to use with caution", sep = " "),
        x = "",
        y = "Number of samples"
      ) +
      
      theme_bw(base_size = 16) +
      theme(plot.subtitle = element_text(color="#8B0000"),
            legend.position = "top")
    
    # combine both plots (this makes sure the x-axix ticks are aligned)
    uploadCombinedPlot <- newDataPlot/backgroundDataPlot
    # print the combined plot to the screen
    uploadDataBarcodePlot <- uploadCombinedPlot + patchwork::plot_layout(heights = c(4,2))
    uploadDataBarcodePlot

  })
  
  # the following renderUI is used to dynamically generate the tabsets when the
  # file is loaded. Until the file is loaded, app will not show the tabset.
  output$tb <- renderUI({
    if(is.null(uploadData()))
      tags$h3("Test your own plasma miR-Seq data using ",
              tags$img(src = "drac.png",
                       heigth = 200,
                       width = 200))
    else
      tabsetPanel(tabPanel("About file",
                           tableOutput("about")),
                  
                  tabPanel("Data",
                           tableOutput("table")),
                  
                  tabPanel("Results Summary",
                           tableOutput("sum")),
                  
                  tabPanel("Distributions", 
                           tags$br(),
                           
                           fluidRow(
                             selectInput(inputId = "select",
                                         label = "Samplename",
                                         choices = colnames(DGEList_public()$counts),
                                         selected = NULL)),
                             fluidRow(
                               plotOutput("distributions")
                             )),
                  
                  tabPanel("Haemolysis Metric",
                           tags$br(),
                           
                           fluidRow(
                           plotOutput("dist_diff"))
                  ))
  })
  
  
}

shinyApp(ui = ui, server = server)

