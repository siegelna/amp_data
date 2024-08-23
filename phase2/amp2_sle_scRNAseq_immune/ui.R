library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")



sc2conf = readRDS("sc2conf.rds")
sc2def  = readRDS("sc2def.rds")



sc3conf = readRDS("sc3conf.rds")
sc3def  = readRDS("sc3def.rds")



sc4conf = readRDS("sc4conf.rds")
sc4def  = readRDS("sc4def.rds")



sc5conf = readRDS("sc5conf.rds")
sc5def  = readRDS("sc5def.rds")



### Start server code 
shinyUI(fluidPage( 
### HTML formatting of error messages 
 
tags$head(tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}"))), 
list(tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"))), 
 
   
### Page title 
titlePanel("AMP Phase II SLE"),  
navbarPage( 
  NULL,  
 navbarMenu("B Cell",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc1a1drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                           selected = sc1def$dimred[1]), 
            selectInput("sc1a1drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                        selected = sc1def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc1a1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc1a1togL % 2 == 1", 
          selectInput("sc1a1sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput("sc1a1sub1.ui"), 
          actionButton("sc1a1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1a1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc1a1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1a1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc1a1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc1a1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc1a1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc1a1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc1a1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a1inp1", "Cell information:", 
                           choices = sc1conf$UI, 
                           selected = sc1def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a1tog1 % 2 == 1", 
              radioButtons("sc1a1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc1a1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc1a1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a1oup1.ui"))), 
        downloadButton("sc1a1oup1.pdf", "Download PDF"), 
        downloadButton("sc1a1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("sc1a1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.sc1a1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("sc1a1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("sc1a1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a1tog2 % 2 == 1", 
              radioButtons("sc1a1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc1a1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("sc1a1oup2.ui"))), 
        downloadButton("sc1a1oup2.pdf", "Download PDF"), 
        downloadButton("sc1a1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc1a2drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                           selected = sc1def$dimred[1]), 
            selectInput("sc1a2drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                        selected = sc1def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc1a2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc1a2togL % 2 == 1", 
          selectInput("sc1a2sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput("sc1a2sub1.ui"), 
          actionButton("sc1a2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1a2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc1a2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1a2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc1a2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc1a2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc1a2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc1a2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc1a2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a2inp1", "Cell information:", 
                           choices = sc1conf$UI, 
                           selected = sc1def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a2tog1 % 2 == 1", 
              radioButtons("sc1a2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc1a2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc1a2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a2oup1.ui"))), 
        downloadButton("sc1a2oup1.pdf", "Download PDF"), 
        downloadButton("sc1a2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a2inp2", "Cell information:", 
                           choices = sc1conf$UI, 
                           selected = sc1def$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a2tog2 % 2 == 1", 
              radioButtons("sc1a2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc1a2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc1a2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a2oup2.ui"))), 
        downloadButton("sc1a2oup2.pdf", "Download PDF"), 
        downloadButton("sc1a2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc1a3drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                           selected = sc1def$dimred[1]), 
            selectInput("sc1a3drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                        selected = sc1def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc1a3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc1a3togL % 2 == 1", 
          selectInput("sc1a3sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput("sc1a3sub1.ui"), 
          actionButton("sc1a3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1a3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc1a3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1a3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc1a3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc1a3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc1a3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc1a3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc1a3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a3tog1 % 2 == 1", 
              radioButtons("sc1a3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc1a3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a3oup1.ui"))), 
        downloadButton("sc1a3oup1.pdf", "Download PDF"), 
        downloadButton("sc1a3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a3tog2 % 2 == 1", 
              radioButtons("sc1a3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc1a3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a3oup2.ui"))), 
        downloadButton("sc1a3oup2.pdf", "Download PDF"), 
        downloadButton("sc1a3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("sc1b2drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                           selected = sc1def$dimred[1]), 
           selectInput("sc1b2drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                       selected = sc1def$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("sc1b2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.sc1b2togL % 2 == 1", 
         selectInput("sc1b2sub1", "Cell information to subset:", 
                     choices = sc1conf[grp == TRUE]$UI, 
                    selected = sc1def$grp1), 
         uiOutput("sc1b2sub1.ui"), 
         actionButton("sc1b2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc1b2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("sc1b2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc1b2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("sc1b2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("sc1b2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("sc1b2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("sc1b2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("sc1b2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("sc1b2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("sc1b2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("sc1b2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.sc1b2tog1 % 2 == 1", 
         radioButtons("sc1b2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("sc1b2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("sc1b2oup1.ui"), 
       downloadButton("sc1b2oup1.pdf", "Download PDF"), 
       downloadButton("sc1b2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("sc1b2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("sc1b2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("sc1b2oup2.ui"), 
       downloadButton("sc1b2oup2.pdf", "Download PDF"), 
       downloadButton("sc1b2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("sc1b2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("sc1c1inp1", "Cell information (X-axis):", 
                   choices = sc1conf[grp == TRUE]$UI, 
                   selected = sc1def$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("sc1c1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("sc1c1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("sc1c1pts", "Show data points", value = FALSE), 
       actionButton("sc1c1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.sc1c1togL % 2 == 1", 
         selectInput("sc1c1sub1", "Cell information to subset:", 
                     choices = sc1conf[grp == TRUE]$UI, 
                     selected = sc1def$grp1), 
         uiOutput("sc1c1sub1.ui"), 
         actionButton("sc1c1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc1c1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("sc1c1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc1c1tog % 2 == 1", 
         sliderInput("sc1c1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("sc1c1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("sc1c1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("sc1c1oup.ui"),  
            downloadButton("sc1c1oup.pdf", "Download PDF"),  
            downloadButton("sc1c1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("sc1c1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("sc1c1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("sc1c2inp1", "Cell information to plot (X-axis):", 
                  choices = sc1conf[grp == TRUE]$UI, 
                  selected = sc1def$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("sc1c2inp2", "Cell information to group / colour by:", 
                  choices = sc1conf[grp == TRUE]$UI, 
                  selected = sc1def$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("sc1c2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("sc1c2flp", "Flip X/Y", value = FALSE), 
      actionButton("sc1c2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.sc1c2togL % 2 == 1", 
        selectInput("sc1c2sub1", "Cell information to subset:", 
                    choices = sc1conf[grp == TRUE]$UI, 
                    selected = sc1def$grp1), 
        uiOutput("sc1c2sub1.ui"), 
        actionButton("sc1c2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("sc1c2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("sc1c2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.sc1c2tog % 2 == 1", 
        radioButtons("sc1c2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("sc1c2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("sc1c2oup.ui"),  
           downloadButton("sc1c2oup.pdf", "Download PDF"),  
           downloadButton("sc1c2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("sc1c2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("sc1c2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("sc1d1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(sc1def$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("sc1d1grp", "Group by:", 
                    choices = sc1conf[grp == TRUE]$UI, 
                    selected = sc1conf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("sc1d1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("sc1d1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("sc1d1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("sc1d1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("sc1d1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc1d1togL % 2 == 1", 
          selectInput("sc1d1sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput("sc1d1sub1.ui"), 
          actionButton("sc1d1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1d1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("sc1d1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1d1tog % 2 == 1", 
          radioButtons("sc1d1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("sc1d1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("sc1d1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("sc1d1oupTxt")), 
             uiOutput("sc1d1oup.ui"), 
             downloadButton("sc1d1oup.pdf", "Download PDF"), 
             downloadButton("sc1d1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("sc1d1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("sc1d1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

navbarMenu("Myeloid",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc2a1drX", "X-axis:", choices = sc2conf[dimred == TRUE]$UI, 
                           selected = sc2def$dimred[1]), 
            selectInput("sc2a1drY", "Y-axis:", choices = sc2conf[dimred == TRUE]$UI, 
                        selected = sc2def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc2a1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc2a1togL % 2 == 1", 
          selectInput("sc2a1sub1", "Cell information to subset:", 
                      choices = sc2conf[grp == TRUE]$UI, 
                      selected = sc2def$grp1), 
          uiOutput("sc2a1sub1.ui"), 
          actionButton("sc2a1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc2a1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc2a1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc2a1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc2a1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc2a1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc2a1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc2a1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc2a1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("sc2a1inp1", "Cell information:", 
                           choices = sc2conf$UI, 
                           selected = sc2def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc2a1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc2a1tog1 % 2 == 1", 
              radioButtons("sc2a1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc2a1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc2a1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc2a1oup1.ui"))), 
        downloadButton("sc2a1oup1.pdf", "Download PDF"), 
        downloadButton("sc2a1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc2a1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc2a1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("sc2a1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.sc2a1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("sc2a1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("sc2a1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("sc2a1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc2a1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc2a1tog2 % 2 == 1", 
              radioButtons("sc2a1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc2a1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("sc2a1oup2.ui"))), 
        downloadButton("sc2a1oup2.pdf", "Download PDF"), 
        downloadButton("sc2a1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc2a1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc2a1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc2a2drX", "X-axis:", choices = sc2conf[dimred == TRUE]$UI, 
                           selected = sc2def$dimred[1]), 
            selectInput("sc2a2drY", "Y-axis:", choices = sc2conf[dimred == TRUE]$UI, 
                        selected = sc2def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc2a2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc2a2togL % 2 == 1", 
          selectInput("sc2a2sub1", "Cell information to subset:", 
                      choices = sc2conf[grp == TRUE]$UI, 
                      selected = sc2def$grp1), 
          uiOutput("sc2a2sub1.ui"), 
          actionButton("sc2a2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc2a2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc2a2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc2a2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc2a2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc2a2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc2a2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc2a2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc2a2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc2a2inp1", "Cell information:", 
                           choices = sc2conf$UI, 
                           selected = sc2def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc2a2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc2a2tog1 % 2 == 1", 
              radioButtons("sc2a2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc2a2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc2a2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc2a2oup1.ui"))), 
        downloadButton("sc2a2oup1.pdf", "Download PDF"), 
        downloadButton("sc2a2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc2a2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc2a2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc2a2inp2", "Cell information:", 
                           choices = sc2conf$UI, 
                           selected = sc2def$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc2a2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc2a2tog2 % 2 == 1", 
              radioButtons("sc2a2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc2a2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc2a2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc2a2oup2.ui"))), 
        downloadButton("sc2a2oup2.pdf", "Download PDF"), 
        downloadButton("sc2a2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc2a2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc2a2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc2a3drX", "X-axis:", choices = sc2conf[dimred == TRUE]$UI, 
                           selected = sc2def$dimred[1]), 
            selectInput("sc2a3drY", "Y-axis:", choices = sc2conf[dimred == TRUE]$UI, 
                        selected = sc2def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc2a3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc2a3togL % 2 == 1", 
          selectInput("sc2a3sub1", "Cell information to subset:", 
                      choices = sc2conf[grp == TRUE]$UI, 
                      selected = sc2def$grp1), 
          uiOutput("sc2a3sub1.ui"), 
          actionButton("sc2a3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc2a3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc2a3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc2a3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc2a3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc2a3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc2a3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc2a3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc2a3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc2a3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc2a3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc2a3tog1 % 2 == 1", 
              radioButtons("sc2a3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc2a3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc2a3oup1.ui"))), 
        downloadButton("sc2a3oup1.pdf", "Download PDF"), 
        downloadButton("sc2a3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc2a3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc2a3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc2a3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc2a3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc2a3tog2 % 2 == 1", 
              radioButtons("sc2a3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc2a3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc2a3oup2.ui"))), 
        downloadButton("sc2a3oup2.pdf", "Download PDF"), 
        downloadButton("sc2a3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc2a3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc2a3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("sc2b2drX", "X-axis:", choices = sc2conf[dimred == TRUE]$UI, 
                           selected = sc2def$dimred[1]), 
           selectInput("sc2b2drY", "Y-axis:", choices = sc2conf[dimred == TRUE]$UI, 
                       selected = sc2def$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("sc2b2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.sc2b2togL % 2 == 1", 
         selectInput("sc2b2sub1", "Cell information to subset:", 
                     choices = sc2conf[grp == TRUE]$UI, 
                    selected = sc2def$grp1), 
         uiOutput("sc2b2sub1.ui"), 
         actionButton("sc2b2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc2b2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("sc2b2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc2b2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("sc2b2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("sc2b2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("sc2b2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("sc2b2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("sc2b2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("sc2b2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("sc2b2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("sc2b2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.sc2b2tog1 % 2 == 1", 
         radioButtons("sc2b2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("sc2b2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("sc2b2oup1.ui"), 
       downloadButton("sc2b2oup1.pdf", "Download PDF"), 
       downloadButton("sc2b2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("sc2b2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("sc2b2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("sc2b2oup2.ui"), 
       downloadButton("sc2b2oup2.pdf", "Download PDF"), 
       downloadButton("sc2b2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("sc2b2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("sc2c1inp1", "Cell information (X-axis):", 
                   choices = sc2conf[grp == TRUE]$UI, 
                   selected = sc2def$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("sc2c1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("sc2c1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("sc2c1pts", "Show data points", value = FALSE), 
       actionButton("sc2c1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.sc2c1togL % 2 == 1", 
         selectInput("sc2c1sub1", "Cell information to subset:", 
                     choices = sc2conf[grp == TRUE]$UI, 
                     selected = sc2def$grp1), 
         uiOutput("sc2c1sub1.ui"), 
         actionButton("sc2c1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc2c1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("sc2c1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc2c1tog % 2 == 1", 
         sliderInput("sc2c1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("sc2c1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("sc2c1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("sc2c1oup.ui"),  
            downloadButton("sc2c1oup.pdf", "Download PDF"),  
            downloadButton("sc2c1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("sc2c1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("sc2c1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("sc2c2inp1", "Cell information to plot (X-axis):", 
                  choices = sc2conf[grp == TRUE]$UI, 
                  selected = sc2def$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("sc2c2inp2", "Cell information to group / colour by:", 
                  choices = sc2conf[grp == TRUE]$UI, 
                  selected = sc2def$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("sc2c2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("sc2c2flp", "Flip X/Y", value = FALSE), 
      actionButton("sc2c2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.sc2c2togL % 2 == 1", 
        selectInput("sc2c2sub1", "Cell information to subset:", 
                    choices = sc2conf[grp == TRUE]$UI, 
                    selected = sc2def$grp1), 
        uiOutput("sc2c2sub1.ui"), 
        actionButton("sc2c2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("sc2c2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("sc2c2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.sc2c2tog % 2 == 1", 
        radioButtons("sc2c2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("sc2c2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("sc2c2oup.ui"),  
           downloadButton("sc2c2oup.pdf", "Download PDF"),  
           downloadButton("sc2c2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("sc2c2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("sc2c2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("sc2d1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(sc2def$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("sc2d1grp", "Group by:", 
                    choices = sc2conf[grp == TRUE]$UI, 
                    selected = sc2conf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("sc2d1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("sc2d1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("sc2d1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("sc2d1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("sc2d1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc2d1togL % 2 == 1", 
          selectInput("sc2d1sub1", "Cell information to subset:", 
                      choices = sc2conf[grp == TRUE]$UI, 
                      selected = sc2def$grp1), 
          uiOutput("sc2d1sub1.ui"), 
          actionButton("sc2d1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc2d1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("sc2d1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc2d1tog % 2 == 1", 
          radioButtons("sc2d1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("sc2d1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("sc2d1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("sc2d1oupTxt")), 
             uiOutput("sc2d1oup.ui"), 
             downloadButton("sc2d1oup.pdf", "Download PDF"), 
             downloadButton("sc2d1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("sc2d1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("sc2d1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

navbarMenu("NK",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc3a1drX", "X-axis:", choices = sc3conf[dimred == TRUE]$UI, 
                           selected = sc3def$dimred[1]), 
            selectInput("sc3a1drY", "Y-axis:", choices = sc3conf[dimred == TRUE]$UI, 
                        selected = sc3def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc3a1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc3a1togL % 2 == 1", 
          selectInput("sc3a1sub1", "Cell information to subset:", 
                      choices = sc3conf[grp == TRUE]$UI, 
                      selected = sc3def$grp1), 
          uiOutput("sc3a1sub1.ui"), 
          actionButton("sc3a1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc3a1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc3a1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc3a1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc3a1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc3a1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc3a1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc3a1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc3a1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("sc3a1inp1", "Cell information:", 
                           choices = sc3conf$UI, 
                           selected = sc3def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc3a1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc3a1tog1 % 2 == 1", 
              radioButtons("sc3a1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc3a1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc3a1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc3a1oup1.ui"))), 
        downloadButton("sc3a1oup1.pdf", "Download PDF"), 
        downloadButton("sc3a1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc3a1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc3a1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("sc3a1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.sc3a1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("sc3a1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("sc3a1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("sc3a1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc3a1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc3a1tog2 % 2 == 1", 
              radioButtons("sc3a1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc3a1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("sc3a1oup2.ui"))), 
        downloadButton("sc3a1oup2.pdf", "Download PDF"), 
        downloadButton("sc3a1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc3a1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc3a1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc3a2drX", "X-axis:", choices = sc3conf[dimred == TRUE]$UI, 
                           selected = sc3def$dimred[1]), 
            selectInput("sc3a2drY", "Y-axis:", choices = sc3conf[dimred == TRUE]$UI, 
                        selected = sc3def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc3a2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc3a2togL % 2 == 1", 
          selectInput("sc3a2sub1", "Cell information to subset:", 
                      choices = sc3conf[grp == TRUE]$UI, 
                      selected = sc3def$grp1), 
          uiOutput("sc3a2sub1.ui"), 
          actionButton("sc3a2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc3a2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc3a2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc3a2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc3a2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc3a2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc3a2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc3a2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc3a2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc3a2inp1", "Cell information:", 
                           choices = sc3conf$UI, 
                           selected = sc3def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc3a2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc3a2tog1 % 2 == 1", 
              radioButtons("sc3a2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc3a2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc3a2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc3a2oup1.ui"))), 
        downloadButton("sc3a2oup1.pdf", "Download PDF"), 
        downloadButton("sc3a2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc3a2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc3a2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc3a2inp2", "Cell information:", 
                           choices = sc3conf$UI, 
                           selected = sc3def$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc3a2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc3a2tog2 % 2 == 1", 
              radioButtons("sc3a2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc3a2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc3a2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc3a2oup2.ui"))), 
        downloadButton("sc3a2oup2.pdf", "Download PDF"), 
        downloadButton("sc3a2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc3a2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc3a2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc3a3drX", "X-axis:", choices = sc3conf[dimred == TRUE]$UI, 
                           selected = sc3def$dimred[1]), 
            selectInput("sc3a3drY", "Y-axis:", choices = sc3conf[dimred == TRUE]$UI, 
                        selected = sc3def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc3a3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc3a3togL % 2 == 1", 
          selectInput("sc3a3sub1", "Cell information to subset:", 
                      choices = sc3conf[grp == TRUE]$UI, 
                      selected = sc3def$grp1), 
          uiOutput("sc3a3sub1.ui"), 
          actionButton("sc3a3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc3a3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc3a3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc3a3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc3a3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc3a3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc3a3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc3a3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc3a3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc3a3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc3a3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc3a3tog1 % 2 == 1", 
              radioButtons("sc3a3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc3a3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc3a3oup1.ui"))), 
        downloadButton("sc3a3oup1.pdf", "Download PDF"), 
        downloadButton("sc3a3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc3a3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc3a3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc3a3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc3a3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc3a3tog2 % 2 == 1", 
              radioButtons("sc3a3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc3a3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc3a3oup2.ui"))), 
        downloadButton("sc3a3oup2.pdf", "Download PDF"), 
        downloadButton("sc3a3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc3a3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc3a3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("sc3b2drX", "X-axis:", choices = sc3conf[dimred == TRUE]$UI, 
                           selected = sc3def$dimred[1]), 
           selectInput("sc3b2drY", "Y-axis:", choices = sc3conf[dimred == TRUE]$UI, 
                       selected = sc3def$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("sc3b2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.sc3b2togL % 2 == 1", 
         selectInput("sc3b2sub1", "Cell information to subset:", 
                     choices = sc3conf[grp == TRUE]$UI, 
                    selected = sc3def$grp1), 
         uiOutput("sc3b2sub1.ui"), 
         actionButton("sc3b2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc3b2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("sc3b2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc3b2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("sc3b2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("sc3b2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("sc3b2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("sc3b2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("sc3b2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("sc3b2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("sc3b2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("sc3b2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.sc3b2tog1 % 2 == 1", 
         radioButtons("sc3b2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("sc3b2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("sc3b2oup1.ui"), 
       downloadButton("sc3b2oup1.pdf", "Download PDF"), 
       downloadButton("sc3b2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("sc3b2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("sc3b2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("sc3b2oup2.ui"), 
       downloadButton("sc3b2oup2.pdf", "Download PDF"), 
       downloadButton("sc3b2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("sc3b2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("sc3c1inp1", "Cell information (X-axis):", 
                   choices = sc3conf[grp == TRUE]$UI, 
                   selected = sc3def$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("sc3c1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("sc3c1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("sc3c1pts", "Show data points", value = FALSE), 
       actionButton("sc3c1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.sc3c1togL % 2 == 1", 
         selectInput("sc3c1sub1", "Cell information to subset:", 
                     choices = sc3conf[grp == TRUE]$UI, 
                     selected = sc3def$grp1), 
         uiOutput("sc3c1sub1.ui"), 
         actionButton("sc3c1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc3c1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("sc3c1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc3c1tog % 2 == 1", 
         sliderInput("sc3c1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("sc3c1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("sc3c1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("sc3c1oup.ui"),  
            downloadButton("sc3c1oup.pdf", "Download PDF"),  
            downloadButton("sc3c1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("sc3c1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("sc3c1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("sc3c2inp1", "Cell information to plot (X-axis):", 
                  choices = sc3conf[grp == TRUE]$UI, 
                  selected = sc3def$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("sc3c2inp2", "Cell information to group / colour by:", 
                  choices = sc3conf[grp == TRUE]$UI, 
                  selected = sc3def$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("sc3c2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("sc3c2flp", "Flip X/Y", value = FALSE), 
      actionButton("sc3c2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.sc3c2togL % 2 == 1", 
        selectInput("sc3c2sub1", "Cell information to subset:", 
                    choices = sc3conf[grp == TRUE]$UI, 
                    selected = sc3def$grp1), 
        uiOutput("sc3c2sub1.ui"), 
        actionButton("sc3c2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("sc3c2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("sc3c2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.sc3c2tog % 2 == 1", 
        radioButtons("sc3c2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("sc3c2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("sc3c2oup.ui"),  
           downloadButton("sc3c2oup.pdf", "Download PDF"),  
           downloadButton("sc3c2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("sc3c2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("sc3c2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("sc3d1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(sc3def$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("sc3d1grp", "Group by:", 
                    choices = sc3conf[grp == TRUE]$UI, 
                    selected = sc3conf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("sc3d1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("sc3d1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("sc3d1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("sc3d1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("sc3d1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc3d1togL % 2 == 1", 
          selectInput("sc3d1sub1", "Cell information to subset:", 
                      choices = sc3conf[grp == TRUE]$UI, 
                      selected = sc3def$grp1), 
          uiOutput("sc3d1sub1.ui"), 
          actionButton("sc3d1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc3d1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("sc3d1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc3d1tog % 2 == 1", 
          radioButtons("sc3d1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("sc3d1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("sc3d1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("sc3d1oupTxt")), 
             uiOutput("sc3d1oup.ui"), 
             downloadButton("sc3d1oup.pdf", "Download PDF"), 
             downloadButton("sc3d1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("sc3d1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("sc3d1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

navbarMenu("Plasma Cell",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc4a1drX", "X-axis:", choices = sc4conf[dimred == TRUE]$UI, 
                           selected = sc4def$dimred[1]), 
            selectInput("sc4a1drY", "Y-axis:", choices = sc4conf[dimred == TRUE]$UI, 
                        selected = sc4def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc4a1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc4a1togL % 2 == 1", 
          selectInput("sc4a1sub1", "Cell information to subset:", 
                      choices = sc4conf[grp == TRUE]$UI, 
                      selected = sc4def$grp1), 
          uiOutput("sc4a1sub1.ui"), 
          actionButton("sc4a1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc4a1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc4a1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc4a1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc4a1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc4a1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc4a1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc4a1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc4a1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("sc4a1inp1", "Cell information:", 
                           choices = sc4conf$UI, 
                           selected = sc4def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc4a1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc4a1tog1 % 2 == 1", 
              radioButtons("sc4a1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc4a1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc4a1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc4a1oup1.ui"))), 
        downloadButton("sc4a1oup1.pdf", "Download PDF"), 
        downloadButton("sc4a1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc4a1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc4a1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("sc4a1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.sc4a1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("sc4a1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("sc4a1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("sc4a1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc4a1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc4a1tog2 % 2 == 1", 
              radioButtons("sc4a1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc4a1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("sc4a1oup2.ui"))), 
        downloadButton("sc4a1oup2.pdf", "Download PDF"), 
        downloadButton("sc4a1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc4a1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc4a1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc4a2drX", "X-axis:", choices = sc4conf[dimred == TRUE]$UI, 
                           selected = sc4def$dimred[1]), 
            selectInput("sc4a2drY", "Y-axis:", choices = sc4conf[dimred == TRUE]$UI, 
                        selected = sc4def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc4a2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc4a2togL % 2 == 1", 
          selectInput("sc4a2sub1", "Cell information to subset:", 
                      choices = sc4conf[grp == TRUE]$UI, 
                      selected = sc4def$grp1), 
          uiOutput("sc4a2sub1.ui"), 
          actionButton("sc4a2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc4a2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc4a2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc4a2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc4a2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc4a2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc4a2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc4a2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc4a2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc4a2inp1", "Cell information:", 
                           choices = sc4conf$UI, 
                           selected = sc4def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc4a2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc4a2tog1 % 2 == 1", 
              radioButtons("sc4a2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc4a2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc4a2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc4a2oup1.ui"))), 
        downloadButton("sc4a2oup1.pdf", "Download PDF"), 
        downloadButton("sc4a2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc4a2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc4a2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc4a2inp2", "Cell information:", 
                           choices = sc4conf$UI, 
                           selected = sc4def$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc4a2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc4a2tog2 % 2 == 1", 
              radioButtons("sc4a2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc4a2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc4a2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc4a2oup2.ui"))), 
        downloadButton("sc4a2oup2.pdf", "Download PDF"), 
        downloadButton("sc4a2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc4a2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc4a2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc4a3drX", "X-axis:", choices = sc4conf[dimred == TRUE]$UI, 
                           selected = sc4def$dimred[1]), 
            selectInput("sc4a3drY", "Y-axis:", choices = sc4conf[dimred == TRUE]$UI, 
                        selected = sc4def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc4a3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc4a3togL % 2 == 1", 
          selectInput("sc4a3sub1", "Cell information to subset:", 
                      choices = sc4conf[grp == TRUE]$UI, 
                      selected = sc4def$grp1), 
          uiOutput("sc4a3sub1.ui"), 
          actionButton("sc4a3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc4a3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc4a3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc4a3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc4a3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc4a3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc4a3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc4a3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc4a3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc4a3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc4a3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc4a3tog1 % 2 == 1", 
              radioButtons("sc4a3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc4a3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc4a3oup1.ui"))), 
        downloadButton("sc4a3oup1.pdf", "Download PDF"), 
        downloadButton("sc4a3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc4a3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc4a3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc4a3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc4a3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc4a3tog2 % 2 == 1", 
              radioButtons("sc4a3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc4a3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc4a3oup2.ui"))), 
        downloadButton("sc4a3oup2.pdf", "Download PDF"), 
        downloadButton("sc4a3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc4a3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc4a3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("sc4b2drX", "X-axis:", choices = sc4conf[dimred == TRUE]$UI, 
                           selected = sc4def$dimred[1]), 
           selectInput("sc4b2drY", "Y-axis:", choices = sc4conf[dimred == TRUE]$UI, 
                       selected = sc4def$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("sc4b2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.sc4b2togL % 2 == 1", 
         selectInput("sc4b2sub1", "Cell information to subset:", 
                     choices = sc4conf[grp == TRUE]$UI, 
                    selected = sc4def$grp1), 
         uiOutput("sc4b2sub1.ui"), 
         actionButton("sc4b2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc4b2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("sc4b2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc4b2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("sc4b2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("sc4b2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("sc4b2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("sc4b2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("sc4b2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("sc4b2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("sc4b2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("sc4b2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.sc4b2tog1 % 2 == 1", 
         radioButtons("sc4b2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("sc4b2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("sc4b2oup1.ui"), 
       downloadButton("sc4b2oup1.pdf", "Download PDF"), 
       downloadButton("sc4b2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("sc4b2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("sc4b2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("sc4b2oup2.ui"), 
       downloadButton("sc4b2oup2.pdf", "Download PDF"), 
       downloadButton("sc4b2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("sc4b2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("sc4c1inp1", "Cell information (X-axis):", 
                   choices = sc4conf[grp == TRUE]$UI, 
                   selected = sc4def$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("sc4c1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("sc4c1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("sc4c1pts", "Show data points", value = FALSE), 
       actionButton("sc4c1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.sc4c1togL % 2 == 1", 
         selectInput("sc4c1sub1", "Cell information to subset:", 
                     choices = sc4conf[grp == TRUE]$UI, 
                     selected = sc4def$grp1), 
         uiOutput("sc4c1sub1.ui"), 
         actionButton("sc4c1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc4c1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("sc4c1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc4c1tog % 2 == 1", 
         sliderInput("sc4c1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("sc4c1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("sc4c1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("sc4c1oup.ui"),  
            downloadButton("sc4c1oup.pdf", "Download PDF"),  
            downloadButton("sc4c1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("sc4c1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("sc4c1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("sc4c2inp1", "Cell information to plot (X-axis):", 
                  choices = sc4conf[grp == TRUE]$UI, 
                  selected = sc4def$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("sc4c2inp2", "Cell information to group / colour by:", 
                  choices = sc4conf[grp == TRUE]$UI, 
                  selected = sc4def$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("sc4c2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("sc4c2flp", "Flip X/Y", value = FALSE), 
      actionButton("sc4c2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.sc4c2togL % 2 == 1", 
        selectInput("sc4c2sub1", "Cell information to subset:", 
                    choices = sc4conf[grp == TRUE]$UI, 
                    selected = sc4def$grp1), 
        uiOutput("sc4c2sub1.ui"), 
        actionButton("sc4c2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("sc4c2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("sc4c2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.sc4c2tog % 2 == 1", 
        radioButtons("sc4c2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("sc4c2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("sc4c2oup.ui"),  
           downloadButton("sc4c2oup.pdf", "Download PDF"),  
           downloadButton("sc4c2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("sc4c2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("sc4c2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("sc4d1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(sc4def$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("sc4d1grp", "Group by:", 
                    choices = sc4conf[grp == TRUE]$UI, 
                    selected = sc4conf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("sc4d1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("sc4d1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("sc4d1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("sc4d1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("sc4d1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc4d1togL % 2 == 1", 
          selectInput("sc4d1sub1", "Cell information to subset:", 
                      choices = sc4conf[grp == TRUE]$UI, 
                      selected = sc4def$grp1), 
          uiOutput("sc4d1sub1.ui"), 
          actionButton("sc4d1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc4d1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("sc4d1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc4d1tog % 2 == 1", 
          radioButtons("sc4d1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("sc4d1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("sc4d1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("sc4d1oupTxt")), 
             uiOutput("sc4d1oup.ui"), 
             downloadButton("sc4d1oup.pdf", "Download PDF"), 
             downloadButton("sc4d1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("sc4d1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("sc4d1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

navbarMenu("T Cell",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc5a1drX", "X-axis:", choices = sc5conf[dimred == TRUE]$UI, 
                           selected = sc5def$dimred[1]), 
            selectInput("sc5a1drY", "Y-axis:", choices = sc5conf[dimred == TRUE]$UI, 
                        selected = sc5def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc5a1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc5a1togL % 2 == 1", 
          selectInput("sc5a1sub1", "Cell information to subset:", 
                      choices = sc5conf[grp == TRUE]$UI, 
                      selected = sc5def$grp1), 
          uiOutput("sc5a1sub1.ui"), 
          actionButton("sc5a1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc5a1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc5a1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc5a1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc5a1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc5a1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc5a1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc5a1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc5a1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("sc5a1inp1", "Cell information:", 
                           choices = sc5conf$UI, 
                           selected = sc5def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc5a1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc5a1tog1 % 2 == 1", 
              radioButtons("sc5a1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc5a1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc5a1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc5a1oup1.ui"))), 
        downloadButton("sc5a1oup1.pdf", "Download PDF"), 
        downloadButton("sc5a1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc5a1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc5a1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("sc5a1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.sc5a1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("sc5a1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("sc5a1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("sc5a1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc5a1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc5a1tog2 % 2 == 1", 
              radioButtons("sc5a1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc5a1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("sc5a1oup2.ui"))), 
        downloadButton("sc5a1oup2.pdf", "Download PDF"), 
        downloadButton("sc5a1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc5a1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc5a1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc5a2drX", "X-axis:", choices = sc5conf[dimred == TRUE]$UI, 
                           selected = sc5def$dimred[1]), 
            selectInput("sc5a2drY", "Y-axis:", choices = sc5conf[dimred == TRUE]$UI, 
                        selected = sc5def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc5a2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc5a2togL % 2 == 1", 
          selectInput("sc5a2sub1", "Cell information to subset:", 
                      choices = sc5conf[grp == TRUE]$UI, 
                      selected = sc5def$grp1), 
          uiOutput("sc5a2sub1.ui"), 
          actionButton("sc5a2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc5a2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc5a2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc5a2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc5a2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc5a2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc5a2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc5a2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc5a2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc5a2inp1", "Cell information:", 
                           choices = sc5conf$UI, 
                           selected = sc5def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc5a2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc5a2tog1 % 2 == 1", 
              radioButtons("sc5a2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc5a2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc5a2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc5a2oup1.ui"))), 
        downloadButton("sc5a2oup1.pdf", "Download PDF"), 
        downloadButton("sc5a2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc5a2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc5a2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc5a2inp2", "Cell information:", 
                           choices = sc5conf$UI, 
                           selected = sc5def$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc5a2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc5a2tog2 % 2 == 1", 
              radioButtons("sc5a2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc5a2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc5a2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc5a2oup2.ui"))), 
        downloadButton("sc5a2oup2.pdf", "Download PDF"), 
        downloadButton("sc5a2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc5a2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc5a2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc5a3drX", "X-axis:", choices = sc5conf[dimred == TRUE]$UI, 
                           selected = sc5def$dimred[1]), 
            selectInput("sc5a3drY", "Y-axis:", choices = sc5conf[dimred == TRUE]$UI, 
                        selected = sc5def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("sc5a3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc5a3togL % 2 == 1", 
          selectInput("sc5a3sub1", "Cell information to subset:", 
                      choices = sc5conf[grp == TRUE]$UI, 
                      selected = sc5def$grp1), 
          uiOutput("sc5a3sub1.ui"), 
          actionButton("sc5a3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc5a3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("sc5a3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc5a3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc5a3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc5a3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc5a3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc5a3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc5a3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("sc5a3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc5a3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc5a3tog1 % 2 == 1", 
              radioButtons("sc5a3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc5a3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc5a3oup1.ui"))), 
        downloadButton("sc5a3oup1.pdf", "Download PDF"), 
        downloadButton("sc5a3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc5a3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc5a3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("sc5a3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc5a3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc5a3tog2 % 2 == 1", 
              radioButtons("sc5a3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc5a3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc5a3oup2.ui"))), 
        downloadButton("sc5a3oup2.pdf", "Download PDF"), 
        downloadButton("sc5a3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc5a3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc5a3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("sc5b2drX", "X-axis:", choices = sc5conf[dimred == TRUE]$UI, 
                           selected = sc5def$dimred[1]), 
           selectInput("sc5b2drY", "Y-axis:", choices = sc5conf[dimred == TRUE]$UI, 
                       selected = sc5def$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("sc5b2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.sc5b2togL % 2 == 1", 
         selectInput("sc5b2sub1", "Cell information to subset:", 
                     choices = sc5conf[grp == TRUE]$UI, 
                    selected = sc5def$grp1), 
         uiOutput("sc5b2sub1.ui"), 
         actionButton("sc5b2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc5b2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("sc5b2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc5b2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("sc5b2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("sc5b2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("sc5b2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("sc5b2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("sc5b2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("sc5b2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("sc5b2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("sc5b2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.sc5b2tog1 % 2 == 1", 
         radioButtons("sc5b2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("sc5b2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("sc5b2oup1.ui"), 
       downloadButton("sc5b2oup1.pdf", "Download PDF"), 
       downloadButton("sc5b2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("sc5b2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("sc5b2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("sc5b2oup2.ui"), 
       downloadButton("sc5b2oup2.pdf", "Download PDF"), 
       downloadButton("sc5b2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("sc5b2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("sc5c1inp1", "Cell information (X-axis):", 
                   choices = sc5conf[grp == TRUE]$UI, 
                   selected = sc5def$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("sc5c1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("sc5c1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("sc5c1pts", "Show data points", value = FALSE), 
       actionButton("sc5c1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.sc5c1togL % 2 == 1", 
         selectInput("sc5c1sub1", "Cell information to subset:", 
                     choices = sc5conf[grp == TRUE]$UI, 
                     selected = sc5def$grp1), 
         uiOutput("sc5c1sub1.ui"), 
         actionButton("sc5c1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc5c1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("sc5c1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc5c1tog % 2 == 1", 
         sliderInput("sc5c1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("sc5c1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("sc5c1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("sc5c1oup.ui"),  
            downloadButton("sc5c1oup.pdf", "Download PDF"),  
            downloadButton("sc5c1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("sc5c1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("sc5c1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("sc5c2inp1", "Cell information to plot (X-axis):", 
                  choices = sc5conf[grp == TRUE]$UI, 
                  selected = sc5def$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("sc5c2inp2", "Cell information to group / colour by:", 
                  choices = sc5conf[grp == TRUE]$UI, 
                  selected = sc5def$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("sc5c2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("sc5c2flp", "Flip X/Y", value = FALSE), 
      actionButton("sc5c2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.sc5c2togL % 2 == 1", 
        selectInput("sc5c2sub1", "Cell information to subset:", 
                    choices = sc5conf[grp == TRUE]$UI, 
                    selected = sc5def$grp1), 
        uiOutput("sc5c2sub1.ui"), 
        actionButton("sc5c2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("sc5c2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("sc5c2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.sc5c2tog % 2 == 1", 
        radioButtons("sc5c2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("sc5c2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("sc5c2oup.ui"),  
           downloadButton("sc5c2oup.pdf", "Download PDF"),  
           downloadButton("sc5c2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("sc5c2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("sc5c2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("sc5d1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(sc5def$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("sc5d1grp", "Group by:", 
                    choices = sc5conf[grp == TRUE]$UI, 
                    selected = sc5conf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("sc5d1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("sc5d1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("sc5d1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("sc5d1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("sc5d1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.sc5d1togL % 2 == 1", 
          selectInput("sc5d1sub1", "Cell information to subset:", 
                      choices = sc5conf[grp == TRUE]$UI, 
                      selected = sc5def$grp1), 
          uiOutput("sc5d1sub1.ui"), 
          actionButton("sc5d1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc5d1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("sc5d1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc5d1tog % 2 == 1", 
          radioButtons("sc5d1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("sc5d1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("sc5d1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("sc5d1oupTxt")), 
             uiOutput("sc5d1oup.ui"), 
             downloadButton("sc5d1oup.pdf", "Download PDF"), 
             downloadButton("sc5d1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("sc5d1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("sc5d1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

   
br(), 
p("", style = "font-size: 125%;"), 
p(em("This webpage was made using "), a("ShinyCell", 
  href = "https://github.com/SGDDNB/ShinyCell",target="_blank")), 
br(),br(),br(),br(),br() 
))) 
 
 
 
 