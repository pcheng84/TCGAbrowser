dashboardPage(
  dashboardHeader(title = "Menu"
    #dropdownMenu(type = "notifications",
    #messageItem(
    #             from = "Dear User",
    #          message = "Please email phil.cheng@usz.ch for questions",
    #        icon = icon("question")))
  ),

  dashboardSidebar(
    sidebarMenu(id = "sidebarmenu",
                menuItem("Cancer", tabName ="cancer", icon = icon("stethoscope", lib="font-awesome")),
                conditionalPanel("input.sidebarmenu === 'cancer'",
                                 selectizeInput("cancername",
                                                "Cancer",
                                                choices = "",
                                                options = list(maxOptions = 35)
                                                )
                #actionButton("goButton2", "Select Cancer")),
                ),

                menuItem("RNAseq explorer", tabName = "a",  icon = icon("bar-chart", lib="font-awesome")),
                conditionalPanel("input.sidebarmenu === 'a'",
                                 selectizeInput("genename",
                                                "Gene (Type gene below)",
                                                choices = "",
                                                options = list(maxOptions = 10)
                                                ),
                actionButton("goButton", "submit"),
                sliderInput(inputId="quantile", label = "Percentile %:", min = 1, max = 50, value = 10, step = 1),
                          radioButtons("Survival", "Survival Time (only overall available)",
                                       c("Overall" = "overall")
                                       #"PFS" = "pfs")),
                                       )
                          #checkboxGroupInput("TNM", "TNM Stage",
                          #                   choices = unique(pat$TNM),
                          #                  selected = unique(pat$TNM))
                )
    )
  ),
  dashboardBody(
    fluidRow(
      tabItems(
        tabItem(tabName = "cancer",
                h2("Welcome to the TCGA Browser v0.9 (beta-testing version)"),
                #h3("Maintenance on Monday January 25 2016 13:00 - 17:00 (UTC +1) "),
                h4("Please email Phil Cheng for any errors and questions."),
                a(href="mailto://phil.cheng@usz.ch", "phil.cheng@usz.ch"),
                h4("Please select a cancer on the left side, then click on RNAseq"),
                tabBox(title="", width=12,
                       id="cancerbox",
                       tabPanel("Datasets loaded",
                                valueBoxOutput("RNAnum", width = 3),
                                valueBoxOutput("mutnum", width = 3),
                                valueBoxOutput("cpnum", width = 3),
                                valueBoxOutput("patnum", width = 3),
                                downloadButton("RNA", "RNAseq table"),
                                downloadButton("Mut", "Mutation matrix"),
                                downloadButton("Cop", "Copy number matrix"),
                                downloadButton("Pat", "Patient table")
                                ),
                                #verbatimTextOutput("test")),
                       tabPanel("Patient data table", dataTableOutput("pattable"))
                      #tabPanel("rna table", dataTableOutput("d1table")),
                      #tabPanel("mutation table", dataTableOutput("m1table")),
                      #tabPanel("copy number table", dataTableOutput("cp1table"))
                )
        ),
        tabItem(tabName = "a",
                infoBoxOutput("cancerselect", width=12),
                fluidRow(
                  box(title = "RNA expression", solidHeader=F, status = "warning", width=4, collapsible = F,
                      plotlyOutput("geneplot", height = 800)),
                  tabBox(title = "", width=8, height = 800,
                         id = "tabbox1",
                         tabPanel("Survival plot", plotOutput("plot", height= 800)),
                         tabPanel("Exome plot", plotOutput("exome", height= 800)),
                         tabPanel("Copy Number Plot", plotOutput("copyplot", height= 800)),
                         #tabPanel("Univariate Cox", plotOutput("coxplot"), verbatimTextOutput("coxinfo")),
                         #tabPanel("Multivariate Cox", verbatimTextOutput("mcoxinfo")),
                         tabPanel("Heatmap", downloadLink("downloadheat", "Download"), plotOutput("heatmap", height=800)),
                         tabPanel("Differential expression", dataTableOutput("DEG")),
                         tabPanel("GSVA Analysis", htmlOutput("GSVA")),
                         tabPanel("GSVA Heatmap", plotOutput("gsvaheatmap", height = 800)),
                         tabPanel("Reactome Plots", selectInput("graph", "Plot",
                                                                choices =c("Dotplot", "Enrichment", "Cnet")),
                                  plotOutput("react", height = 800)),
                         tabPanel("Bar graphs",
                                  selectInput("patgene",
                                              "Choose a factor:",
                                              choices = c("Age",
                                                          "Gender",
                                                          "Stage")
                                              ),
                                  showOutput("bargraph", "nvd3")
                                )
                       )
                )

        )
      )
    )
  )
)
