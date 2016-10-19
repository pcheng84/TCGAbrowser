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

  output$patnum <- renderValueBox({
    valueBox(paste(dim(pat())[1]), subtitle = "Clinical data samples", color = "black")
  })


  output$RNA <- downloadHandler(
    filename = function() {
      paste("RNAseq_", Sys.Date(), ".txt", sep="")
      },
    content = function(con) {
      write.table(d1(), con, sep="\t", row.names=F)
    }
  )

  output$Mut <- downloadHandler(
    filename = function() {
      paste("Mutation_", Sys.Date(), ".txt", sep="")
    },
    content = function(con) {
      write.table(m1(), con, sep="\t", row.names=F)
    }
  )

  output$Cop <- downloadHandler(
    filename = function() {
      paste("CopyNumber_", Sys.Date(), ".txt", sep="")
    },
    content = function(con) {
      write.table(cp1(), con, sep="\t", row.names=F)
    }
  )

  output$Pat <- downloadHandler(
    filename = function() {
      paste("Patient_", Sys.Date(), ".txt", sep="")
    },
    content = function(con) {
      write.table(pat(), con, sep="\t", row.names=F)
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

  output$m1table <- renderDataTable(m1(),
                                    extensions = c('ColVis', 'Responsive'),
                                    options = list(dom = 'C<"clear">lfrtip', pageLength = 10)
                                    )

  updateSelectizeInput(session, "genename", choices = gene.name, server = T, selected = "EGFR")

  gene <- reactive({
    if(is.null(input$genename))
      return("EGFR")
    else {
      input$goButton
      toupper(isolate(input$genename))
    }
  })

  pat.d1 <- reactive({
    d1 <- d1()
    pat <- pat()
    #tumor <- input$TNM
    #d2 <- d1[,c("Gene", na.omit(pat$name[which(pat$TNM %in% tumor)])), with=FALSE]
    d2 <- d1[, c("Gene", colnames(d1)[colnames(d1) %in% pat$name]), with=F]
    setkey(d2, Gene)
    d2
  })

  pat.d1.gene <- reactive({
    patsubset(pat(), d1(), gene(), input$quantile)
  })

  output$geneplot <- renderPlotly({
    plotlygenelevel(pat.d1.gene())


  })

  time <- reactive({
    switch(input$Survival,
           "overall" = "Overall"
           #"pfs" = "PFS"
           )
  })

  output$plot <- renderPlot({
    genesurv(pat.d1.gene(), gene())
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

  hgenes <- reactive({
    pheno <- quan()
    pat.d1 <- pat.d1()
    d2<- pat.d1[, pheno[, name], with=F]
    d2[, Gene := gene.name]
    setkey(d2, Gene)
    d2
  })

  hmap <- reactive({
    withProgress(message = "Calculating Differential Expression", value=0.1, {
      Sys.sleep(0.25)

      phenos <- quan()
      d2 <- hgenes()
      d3 <- DGEList(counts=d2[,setdiff(colnames(d2), "Gene"), with=F], genes=gene.name, group=phenos$gene)
      isexpr <- rowSums(cpm(d3)>1) >= (ncol(d3)/2) #only keeps genes with at least 1 count-per-million in at least half the samples
      d3 <- d3[isexpr,]
      design <- model.matrix(~phenos$gene)
      v2 <- voom(d3, design, plot=F)
      fit <- lmFit(v2, design)
      fit2 <- eBayes(fit)
      fit3 <- topTable(fit2, coef=2, n=Inf, adjust.method="BH", p.value=0.05, lfc=1,sort="p")
      jobLength = 10
      for (i in 1:jobLength) {
        # Do work
        incProgress(0.1, detail = paste("part", i))
      }
      fit3
    })
  })

  #use orior count for log transformation
  output$heatmap <- renderPlot({
    pat.gene <- quan()
    genes <- hgenes()
    deg <- hmap()
    setkey(genes, Gene)
    genes2 <- cpm(genes[, setdiff(colnames(genes), "Gene"), with=F], log=T, normalized.lib.sizes=F)
    rownames(genes2) <- genes$Gene
    map <- genes2[match(na.omit(deg$genes[1:100]), rownames(genes2)),]
    colsidecolors <- matrix((pat.gene$gene * (-1)) + 3)
    colnames(colsidecolors) <- gene()
    #heatmap(map, col=colorRampPalette(c("blue", "white", "red"))(100), labCol=F, scale="row", cexRow= 0.3, ColSideColors = as.character(pat.gene$gene+1))
    heatmap3(map, ColSideColors = colsidecolors, cexRow=0.2, cexCol=0.2,
             legendfun=function() showLegend(legend=c("High", "Low"), col=c("black","firebrick")),verbose=TRUE)
  })

  output$downloadheat <- downloadHandler(


       filename = function() {
         gene <- gene()
         paste(Sys.Date(), cancers[input$cancername, Cancer], gene, input$quantile, '.tiff', sep="_")
       },
       content = function(con) {
         pat.gene <- quan()
         genes <- hgenes()
         deg <- hmap()
         setkey(genes, Gene)
         genes2 <- cpm(genes[, setdiff(colnames(genes), "Gene"), with=F], log=TRUE, normalized.lib.sizes=F)
         rownames(genes2) <- genes$Gene
         map <- genes2[match(na.omit(deg$genes[1:100]), rownames(genes2)),]
         tiff(con, width = 1200, height = 1600, units = "px", pointsize=18)
         par(mar=c(7,4,10,2))
         colsidecolors <- matrix((pat.gene$gene * (-1)) + 3)
         colnames(colsidecolors) <- gene()
         #names(colsidescolors) <- gene()
         #heatmap(map, col=colorRampPalette(c("blue", "white", "red"))(100), labCol=F, scale="row", cexRow= 0.3, ColSideColors = as.character(pat.gene$gene+1))
         heatmap3(map, ColSideColors = colsidecolors, cexRow=1, cexCol=1,
                  legendfun=function() showLegend(legend=c("High", "Low"), col=c("black","firebrick"), cex=3),verbose=TRUE)
         dev.off()

       }
     )


  keggfactor <- reactive({
    switch(input$highlow,
           Upregulated = "greater",
           Downregulated = "less")
  })

  output$Kegg <- renderGvis({
    deg <- hmap()
    limma.fc <- deg$logFC
    names(limma.fc) <- lookup$entrez[match(deg$genes, lookup$gene)]
    kf <- keggfactor()
    fc.kegg.p <- gage(limma.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
    sel <- fc.kegg.p$greater[, "p.val"] < 0.05 & !is.na(fc.kegg.p$greater[,"p.val"])
    greater <- data.frame(cbind(Pathway = rownames(fc.kegg.p$greater[sel,]),round(fc.kegg.p$greater[sel,1:5],5)))
    sel.1 <- fc.kegg.p$less[,"p.val"] < 0.05 & !is.na(fc.kegg.p$less[,"p.val"])
    less <-data.frame(cbind(Pathway = rownames(fc.kegg.p$less[sel.1,]), round(fc.kegg.p$less[sel.1,1:5],5)))
    if(kf == "greater") gvisTable(greater,options=list(page='enable', pageSize= 20, width=950))
    else gvisTable( less,options=list(page='enable', pageSize= 20, width=950))

  })

  output$download_KEGG <- renderUI({

  })

  output$DEG <- renderDataTable(hmap(),
                                extensions = c('TableTools'),
                                options = list(dom = 'T<"clear">lfrtip',
                                               pageLength = 10,
                                               lengthMenu = c(10,25,100,1000, 2000),
                                               tableTools = list(sSwfPath = copySWF("www")))
                                )


  graphfactor <- reactive({

    switch(input$patgene,
           Age = "age_group",
           Gender ="gender",
           Stage = "pathologicstage"
           #Gleason = "something"
    )
  })

  output$bargraph <- renderChart({
    gene <- gene()
    pat.gene <- quan()
    group <- graphfactor()

    gr <- pat.gene[, .N , by=.(gene, with(pat.gene, get(group)))][order(with)]
    setkey(gr, gene, with)
    setnames(gr, 2, group)
    gr[gene == 1, gene2 := "low"]
    gr[gene == 2, gene2 := "high"]
    n1 <- nPlot(N ~ gene2, group = group, data = gr, type = "multiBarChart")
    n1$chart(margin = list(left = 100))
    n1$yAxis(axisLabel  ="Number of patients")
    n1$addParams(dom = "bargraph")
    return(n1)
  })

}
