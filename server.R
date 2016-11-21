function(input, output, session) {

  updateSelectizeInput(session, "cancername", choices = cancers$Cancer_name, server = T, selected = "Skin Cutaneous Melanoma")



    datasets <- reactive({
    if(is.na(input$cancername)) input$cancername <- "SKCM" else{

    sessionEnvir <- sys.frame()
    if (!is.na(input$cancername)){
     load(paste0("./cancers/", cancers[input$cancername, Cancer],"/multi.RData"), sessionEnvir)
      }}
    })

  output$cancerselect <- renderInfoBox({
    infoBox(title = input$cancername, value = cancers[input$cancername, Cancer], subtitle = "RNAseq explorer", icon = icon("area-chart"), color = "aqua", width = 6)
    })

  output$RNAnum <- renderValueBox({
    valueBox(paste(dim(d1())[2]), subtitle = "RNAseq samples", color = "blue")
  })


  output$mutnum <- renderValueBox({
    if(length(m1()) == 0) {
    valueBox(paste("0"), subtitle = "Exome seq samples", color = "red")
    } else {
    valueBox(paste(dim(m1())[2]), subtitle = "Exome seq samples", color = "red")
    }
  })

  output$cpnum <- renderValueBox({
    valueBox(paste(dim(cp1())[2]), subtitle = "Copy Number samples", color = "purple")
  })

  output$rppanum <- renderValueBox({
    valueBox(paste(dim(p1())[2]), subtitle = "RPPA samples", color = "orange")
  })

  output$patnum <- renderValueBox({
    valueBox(paste(dim(pat())[1]), subtitle = "Clinical data samples", color = "black")
  })


  output$RNA <- downloadHandler(
    filename = function() {
      paste("RNAseq_", Sys.Date(), ".txt", sep = "")
      },
    content = function(con) {
      write.table(d1(), con, sep = "\t", row.names=F)
    }
  )

  output$Mut <- downloadHandler(
    filename = function() {
      paste("Mutation_", Sys.Date(), ".txt", sep = "")
    },
    content = function(con) {
      write.table(m1(), con, sep = "\t", row.names=FALSE)
    }
  )

  output$Cop <- downloadHandler(
    filename = function() {
      paste("CopyNumber_", Sys.Date(), ".txt", sep="")
    },
    content = function(con) {
      write.table(cp1(), con, sep = "\t", row.names = FALSE)
    }
  )

  output$Pat <- downloadHandler(
    filename = function() {
      paste("Patient_", Sys.Date(), ".txt", sep = "")
    },
    content = function(con) {
      write.table(pat(), con, sep = "\t", row.names = FALSE)
    }
  )

  output$RPPA <- downloadHandler(
    filename = function() {
      paste("RPPA_", Sys.Date(), ".txt", sep="")
    },
    content = function(con) {
      write.table(p1(), con, sep = "\t", row.names = FALSE)
    }
  )

  testpat <- reactive({
    datasets()
    pat <- combi[[1]]
    setkey(pat, bcr_patient_barcode, name)
    #pat[, TNM := toupper(paste(pat$pathologyTstage, pat$pathologyNstage, pat$pathologyMstage))]
    pat$vitalstatus <- as.numeric(pat$vitalstatus)
    pat$yearstobirth <- as.numeric(pat$yearstobirth)
    pat <- pat[!(is.na(pat$name))]
    pat$gender <- factor(pat$gender, levels=c("male", "female"))
    setkey(pat, name)
    pat
  })

  output$pattable <- renderDataTable(testpat(), extensions = c('ColVis', 'Responsive'), options = list(dom = 'C<"clear">lfrtip', pageLength = 10))


  pat <<- reactive({
    if(!is.null(datasets())) {

      datasets()
      pat <- combi[[1]]
      setkey(pat, bcr_patient_barcode, name)
      pat$vitalstatus <- as.numeric(pat$vitalstatus)
      pat$yearstobirth <- as.numeric(pat$yearstobirth)
      pat$daystosubmittedspecimendx <- as.numeric(pat$daystosubmittedspecimendx)
      pat <- pat[!(is.na(pat$name))]
      pat <- pat[!(is.na(pat$days))]
      pat[, years := round(days/365.25, 2)]
      pat[, TCGA_day := days - daystosubmittedspecimendx]
      pat[, TCGA_year := round(TCGA_day/365.25,2)]
      pat <- pat[!is.na(pat$TCGA_day)]
      pat$gender <- factor(pat$gender, levels=c("male", "female"))
      setkey(pat, name)
      pat$age_group <- cut(pat$yearstobirth + pat$daystosubmittedspecimendx/365.25, br=seq(from = 9, to = 99, by = 10), labels = c("10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99"))
      pat$age_group <- factor(pat$age_group, levels=(sort(unique(pat$age_group))))
      pat
    }
    pat
  })

  d1 <<- reactive({
    if(!is.null(datasets())) {
      datasets()
      d1 <- combi[[2]]
      setkey(d1, Gene)
      d1
    }
  })

  output$d1table <- renderDataTable(d1(), extensions = c('ColVis', 'Responsive'), options = list(dom = 'C<"clear">lfrtip', pageLength = 10))

  cp1 <<- reactive({
    if(!is.null(datasets())) {
      datasets()
      cp1 <- combi[[3]]
      setkey(cp1, Gene)
    }
  })

  output$cp1table <- renderDataTable(cp1(), extensions = c('ColVis', 'Responsive'), options = list(dom = 'C<"clear">lfrtip', pageLength = 10))

  m1 <<- reactive({
    if(!is.null(datasets())) {
      datasets()
      if(length(combi[[4]]) == 0) {
        m1 <- NULL
      } else {
      m1 <- combi[[4]]
      setkey(m1, Gene)
      }
      m1
    }
  })

  p1 <<- reactive({
    if(!is.null(datasets())) {
      datasets()
      if(length(combi[[5]]) == 0) {
        p1 <- NULL
      } else {
        p1 <- combi[[5]]
        setkey(p1, Gene)
      }
      p1
      }
    })

  output$m1table <- renderDataTable(m1(),
                                    extensions = c('ColVis', 'Responsive'),
                                    options = list(dom = 'C<"clear">lfrtip', pageLength = 10)
                                    )

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

  output$downloadheat <- downloadHandler(
       filename = function() {
         gene <- gene()
         paste(Sys.Date(), cancers[input$cancername, Cancer], gene, input$quantile, '.tiff', sep="_")
       },
       content = function(con) {
         tiff(con, width = 1200, height = 1600, units = "px", pointsize=18)
         par(mar=c(7,4,10,2))
         rnaheat(pat.d1.gene(), d1(), hmap(), gene())
         dev.off()
       }
     )

  output$DEG <- renderDataTable(hmap(),
                                extensions = c('TableTools'),
                                options = list(dom = 'T<"clear">lfrtip',
                                               pageLength = 10,
                                               lengthMenu = c(10,25,100,1000, 2000),
                                               tableTools = list(sSwfPath = copySWF("www")))
  )

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

  output$cancerselect2 <- renderInfoBox({
    infoBox(title = input$cancername, value = cancers[input$cancername, Cancer], subtitle = "Mutation explorer", icon = icon("area-chart"), color = "aqua", width = 6)
  })


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

  output$downloadheat2 <- downloadHandler(
    filename = function() {
      gene <- gene()
      paste(Sys.Date(), cancers[input$cancername, Cancer], gene2, '.tiff', sep="_")
    },
    content = function(con) {
      tiff(con, width = 1200, height = 1600, units = "px", pointsize=18)
      par(mar=c(7,4,10,2))
      rnaheat(pat.d1.gene2(), d1(), hmap2(), gene2())
      dev.off()
    }
  )

  output$DEG2 <- renderDataTable(hmap2(),
                                extensions = c('TableTools'),
                                options = list(dom = 'T<"clear">lfrtip',
                                               pageLength = 10,
                                               lengthMenu = c(10,25,100,1000, 2000),
                                               tableTools = list(sSwfPath = copySWF("www")))
  )

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
