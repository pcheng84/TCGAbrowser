#' tcgashiny function
#'
#' Launches local Shiny app to explore your dataset
#'
#' @param pat data frame containing clinical data
#' @param rna data frame containing RNAseq counts
#' @param mut data frame containing binary matrix of mutated genes
#' @param cp data frame containing copy number values
#' @param rppa data frame containing RPPA values
#'
#' @import data.table
#' @import shiny
#' @import shinydashboard
#' @import plotly
#'
#'
#' @return shiny app that explores the 5 datasets
#'
#' @examples
#' data(skcm)
#' tcgashiny(pat, rna, mut, cp, rppa)
#'
#' @export
#'

tcgashiny <- function(pat, rna, mut, cp, rppa) {
  shinyApp(ui = dashboardPage(
    dashboardHeader(title = "Menu"
    ),
    dashboardSidebar(
      sidebarMenu(id = "sidebarmenu",
                  menuItem("RNAseq explorer", tabName = "a",  icon = icon("bar-chart", lib="font-awesome")),
                  conditionalPanel("input.sidebarmenu === 'a'",
                                   selectizeInput("genename",
                                                  "Gene (Type gene below)",
                                                  choices = "",
                                                  options = list(maxOptions = 30)
                                   ),
                                   actionButton("goButton", "submit"),
                                   sliderInput(inputId="quantile",
                                               label = "Percentile %:",
                                               min = 1,
                                               max = 50,
                                               value = 10,
                                               step = 1),
                                   radioButtons("Survival",
                                                "Survival Time (only overall available)",
                                                c("Overall" = "overall")
                                                #"PFS" = "pfs")),
                                   )
                                   #checkboxGroupInput("TNM", "TNM Stage",
                                   #choices = unique(pat$TNM),
                                   #selected = unique(pat$TNM))
                  ),
                  menuItem("Mutation explorer", tabName = "b",  icon = icon("link", lib="font-awesome")),
                  conditionalPanel("input.sidebarmenu === 'b'",
                                   selectizeInput("genename2",
                                                  "Gene (Type gene below)",
                                                  choices = "",
                                                  options = list(maxOptions = 30)
                                   ),
                                   actionButton("goButton2", "submit"),
                                   radioButtons("Survival2",
                                                "Survival Time (only overall available)",
                                                c("Overall" = "overall")
                                                #"PFS" = "pfs")),
                                   )
                                   #checkboxGroupInput("TNM", "TNM Stage",
                                   #choices = unique(pat$TNM),
                                   #selected = unique(pat$TNM))
                  )
      )
    ),
    dashboardBody(
      fluidRow(
        tabItems(
          tabItem(tabName = "cancer",
                  h2("Welcome to the TCGABrowser Shiny App"),
                  #h3("Maintenance on Monday January 25 2016 13:00 - 17:00 (UTC +1) "),
                  tabBox(title="", width=12,
                         id="cancerbox",
                         tabPanel("Datasets loaded",
                                  fluidRow(
                                    valueBoxOutput("RNAnum", width = 2),
                                    valueBoxOutput("mutnum", width = 2),
                                    valueBoxOutput("cpnum", width = 2),
                                    valueBoxOutput("patnum", width = 2),
                                    valueBoxOutput("rppanum", width = 2)
                                  )

                         )

          )),
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
                           tabPanel("Heatmap", plotOutput("heatmap", height=800)),
                           tabPanel("Differential expression", dataTableOutput("DEG")),
                           tabPanel("GSVA Analysis", htmlOutput("GSVA")),
                           tabPanel("GSVA Heatmap", plotOutput("gsvaheatmap", height = 800)),
                           tabPanel("Reactome Plots", selectInput("graph", "Plot",
                                                                  choices =c("Dotplot", "Enrichment", "Cnet")),
                                    plotOutput("react", height = 800)),
                           tabPanel("RPPA Heatmap", plotOutput("rppaheat", height = 800)),
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
          ),

          tabItem(tabName = "b",
                  infoBoxOutput("cancerselect2", width=12),
                  fluidRow(
                    box(title = "RNA expression", solidHeader=F, status = "warning", width=4, collapsible = F,
                        plotlyOutput("geneplot2", height = 800)),
                    tabBox(title = "", width=8, height = 800,
                           id = "tabbox2",
                           tabPanel("Survival plot", plotOutput("plot2", height= 800)),
                           tabPanel("Exome plot", plotOutput("exome2", height= 800)),
                           tabPanel("Copy Number Plot", plotOutput("copyplot2", height= 800)),
                           #tabPanel("Univariate Cox", plotOutput("coxplot"), verbatimTextOutput("coxinfo")),
                           #tabPanel("Multivariate Cox", verbatimTextOutput("mcoxinfo")),
                           tabPanel("Heatmap", plotOutput("heatmap2", height=800)),
                           tabPanel("Differential expression", dataTableOutput("DEG2")),
                           tabPanel("GSVA Analysis", htmlOutput("GSVA2")),
                           tabPanel("GSVA Heatmap", plotOutput("gsvaheatmap2", height = 800)),
                           tabPanel("Reactome Plots", selectInput("graph2", "Plot",
                                                                  choices =c("Dotplot", "Enrichment", "Cnet")),
                                    plotOutput("react2", height = 800)),
                           tabPanel("Bar graphs",
                                    selectInput("patgene2",
                                                "Choose a factor:",
                                                choices = c("Age",
                                                            "Gender",
                                                            "Stage")
                                    ),
                                    showOutput("bargraph2", "nvd3")
                           )
                    )
                  )
          )
        )
      )
    )
  ),
  server = function(input, output, session) {
    library(data.table)
    d1 <- rna
    m1 <- mut
    cp1 <- cp
    p1 <- rppa
    pat <- pat
    gene.name <- unique(d1$Gene)
    gene.name2 <- unique(m1$Gene)

    output$RNAnum <- renderValueBox({
      valueBox(paste(dim(d1)[2], subtitle = "RNAseq samples", color = "blue"))
    })


    output$mutnum <- renderValueBox({
      if(length(m1()) == 0) {
        valueBox(paste("0"), subtitle = "Exome seq samples", color = "red")
      } else {
        valueBox(paste(dim(m1)[2], subtitle = "Exome seq samples", color = "red"))
      }
    })

    output$cpnum <- renderValueBox({
      valueBox(paste(dim(cp1)[2], subtitle = "Copy Number samples", color = "purple"))
    })

    output$rppanum <- renderValueBox({
      valueBox(paste(dim(p1)[2], subtitle = "RPPA samples", color = "orange"))
    })

    output$patnum <- renderValueBox({
      valueBox(paste(dim(pat)[1], subtitle = "Clinical data samples", color = "black"))
    })


    testpat <- reactive({
      pat <- pat
      setkey(pat, bcr_patient_barcode, name)
      #pat[, TNM := toupper(paste(pat$pathologyTstage, pat$pathologyNstage, pat$pathologyMstage))]
      pat$vitalstatus <- as.numeric(pat$vitalstatus)
      pat$yearstobirth <- as.numeric(pat$yearstobirth)
      pat <- pat[!(is.na(pat$name))]
      pat$gender <- factor(pat$gender, levels=c("male", "female"))
      setkey(pat, name)
      pat
    })


    pat <<- reactive({
      pat
    })

    d1 <<- reactive({
        d1 <- rna
        setkey(d1, Gene)
        d1
    })


    cp1 <<- reactive({
        cp1 <- cp
        setkey(cp1, Gene)
        cp1
    })


    m1 <<- reactive({
          m1 <- m1
          setkey(m1, Gene)
          m1
    })

    p1 <<- reactive({
          p1 <- data.table(rppa)
          p1.names <- gsub("(TCGA-.*?-.*?-.*?)-.*", "\\1", colnames(p1))
          setnames(p1, p1.names)
          setkey(p1, Gene)
          p1
    })

    updateSelectizeInput(session, "genename", choices = gene.name, server = T, selected = "EGFR")

    gene <- reactive({
      if(is.null(input$genename))
        return("EGFR")
      else {
        toupper(input$genename)
      }
    })


    pat.d1.gene <- reactive({
      input$goButton
      rnasubset(pat(), d1(), isolate(gene()), isolate(input$quantile))
    })

    output$geneplot <- renderPlotly({
      input$goButton
      plotlygenelevel(isolate(pat.d1.gene()))


    })

    output$plot <- renderPlot({
      input$goButton
      genesurv(pat.d1.gene(), isolate(gene()))
    })

    cox.surv <- reactive({
      gene <- gene()
      pat.gene <- quan()
      time <- time()
      cox <- if(time == "Overall") coxph(Surv(days, vitalstatus) ~ gene, data = pat.gene) else coxph(Surv(pfs_days, pfs) ~ gene, data = pat.gene)
      cox
    })

    mcox.surv <- reactive({
      gene <- gene()
      pat.gene <- quan()
      time <- time()
      mcox <- if(time == "Overall") coxph(Surv(days, vitalstatus) ~(gene + yearstobirth), data = pat.gene) else coxph(Surv(pfs_days, pfs) ~(gene + yearstobirth + gleasonscore), data = pat.gene)
      mcox
    })

    output$coxplot <- renderPlot({
      cox <- cox.surv()
      plot(survfit(cox), xlab = "Survival Time (Days)", ylab = "Survival")
    })

    output$coxinfo <- renderPrint({
      cox <- cox.surv()
      print(summary(cox))
    })

    output$mcoxinfo <- renderPrint({
      mcox <- mcox.surv()
      print(summary(mcox))
    })

    output$exome <- renderPlot({
      gene.mut <- diffmut(pat.d1.gene(), m1())
      plotlymut(pat.d1.gene(), m1(), gene.mut, gene())
    })

    output$copyplot <- renderPlot({
      gene.cp <- diffcp(pat.d1.gene(), cp1())
      plotlycp(pat.d1.gene(), cp1(), gene.cp, gene())
    })


    hmap <- reactive({
      withProgress(message = "Calculating Differential Expression", value=0.1, {
        # Number of times we'll go through the loop
        n <- 10

        for (i in 1:n) {
          # Each time through the loop, add another row of data. This is
          # a stand-in for a long-running computation.

          Sys.sleep(0.25)
          # Increment the progress bar, and update the detail text.
          incProgress(1/n, detail = paste("Calculating...", i))
        }
      })
      fit3 <- rnadeg(pat.d1.gene(), d1())
      fit3
    })

    #use orior count for log transformation
    output$heatmap <- renderPlot({
      input$goButton
      rnaheat(pat.d1.gene(), d1(), hmap(), isolate(gene()))

    })

    output$DEG <- renderDataTable(hmap())

    gsva <- reactive({
      gene.gsva <- rnagsva(pat.d1.gene(), d1())
    })

    gsvasig <- reactive({
      gene.gsva.sig <- rnagsvasig(pat.d1.gene(), gsva())
    })

    output$GSVA <- renderGvis({

      gvisTable(gsvasig(),options=list(page='enable',
                                       pageSize= 20,
                                       height = "automatic",
                                       width = "automatic"))
    })

    output$gsvaheatmap <- renderPlot({
      rnagsvaheat(pat.d1.gene(), gsva(), gsvasig(), gene())

    })

    output$react <- renderPlot({
      gene.react <- rnareact(hmap(), lookup)
      graph <- switch(input$graph,
                      Dotplot = "dot",
                      Enrichment = "map",
                      Cnet = "cnet")
      plotreact(gene.react, hmap(), graph, 15)
    })


    output$rppaheat <- renderPlot({
      rppaheat(pat.d1.gene(), p1(), gene())
    })

    graphfactor <- reactive({
      switch(input$patgene,
             Age = "age_group",
             Gender ="gender",
             Stage = "pathologicstage"
             #Gleason = "something"
      )
    })

    output$bargraph <- renderChart({
      pat.gene <- pat.d1.gene()
      group <- graphfactor()

      gr <- pat.gene[, .N , by=.(gene2, with(pat.gene, get(group)))][order(with)]
      setkey(gr, gene2, with)
      gr2 <- gr[gr$gene2 !="middle"]
      setnames(gr2, 2, group)

      n1 <- nPlot(N ~ gene2, group = group, data = gr2, type = "multiBarChart")
      n1$chart(margin = list(left = 100))
      n1$yAxis(axisLabel  ="Number of patients")
      n1$addParams(dom = "bargraph")
      return(n1)
    })

    #Exome explorer section

    updateSelectizeInput(session, "genename2", choices = gene.name2, server = T, selected = "EGFR")

    gene2 <- reactive({
      if(is.null(input$genename2))
        return("EGFR")
      else {
        toupper(input$genename2)
      }
    })


    pat.d1.gene2 <- reactive({
      input$goButton2
      mutsubset(pat(), d1(), m1(), isolate(gene2()))
    })

    output$geneplot2 <- renderPlotly({
      input$goButton2
      plotlygenelevel(isolate(pat.d1.gene2()))


    })

    output$plot2 <- renderPlot({
      input$goButton2
      genesurv(pat.d1.gene2(), isolate(gene2()))
    })

    output$exome2 <- renderPlot({
      gene.mut <- diffmut(pat.d1.gene2(), m1())
      plotlymut(pat.d1.gene2(), m1(), gene.mut, gene2())
    })

    output$copyplot2 <- renderPlot({
      gene.cp <- diffcp(pat.d1.gene2(), cp1())
      plotlycp(pat.d1.gene2(), cp1(), gene.cp, gene2())
    })


    hmap2 <- reactive({
      withProgress(message = "Calculating Differential Expression", value=0.1, {
        # Number of times we'll go through the loop
        n <- 10

        for (i in 1:n) {
          # Each time through the loop, add another row of data. This is
          # a stand-in for a long-running computation.

          Sys.sleep(0.25)
          # Increment the progress bar, and update the detail text.
          incProgress(1/n, detail = paste("Calculating...", i))
        }
      })
      fit3 <- rnadeg(pat.d1.gene2(), d1())
      fit3
    })

    #use orior count for log transformation
    output$heatmap2 <- renderPlot({
      input$goButton2
      rnaheat(pat.d1.gene2(), d1(), hmap2(), isolate(gene2()))

    })


    output$DEG2 <- renderDataTable(hmap2())

    gsva2 <- reactive({
      gene.gsva <- rnagsva(pat.d1.gene2(), d1())
    })

    gsvasig2 <- reactive({
      gene.gsva.sig <- rnagsvasig(pat.d1.gene2(), gsva2())
    })

    output$GSVA2 <- renderGvis({

      gvisTable(gsvasig2(),options=list(page='enable',
                                        pageSize= 20,
                                        height = "automatic",
                                        width = "automatic"))
    })

    output$gsvaheatmap2 <- renderPlot({
      rnagsvaheat(pat.d1.gene2(), gsva2(), gsvasig2(), gene2())

    })

    output$react2 <- renderPlot({
      gene.react <- rnareact(hmap2(), lookup)
      graph <- switch(input$graph2,
                      Dotplot = "dot",
                      Enrichment = "map",
                      Cnet = "cnet")
      plotreact(gene.react, hmap2(), graph, 15)
    })

    graphfactor2 <- reactive({
      switch(input$patgene2,
             Age = "age_group",
             Gender ="gender",
             Stage = "pathologicstage"
             #Gleason = "something"
      )
    })

    output$bargraph2 <- renderChart({
      pat.gene <- pat.d1.gene2()
      group <- graphfactor2()

      gr <- pat.gene[, .N , by=.(gene2, with(pat.gene, get(group)))][order(with)]
      setkey(gr, gene2, with)
      gr2 <- gr[gr$gene2 !="middle"]
      setnames(gr2, 2, group)

      n2 <- nPlot(N ~ gene2, group = group, data = gr2, type = "multiBarChart")
      n2$chart(margin = list(left = 100))
      n2$yAxis(axisLabel  ="Number of patients")
      n2$addParams(dom = "bargraph")
      return(n2)
    })

  }



  )
}
