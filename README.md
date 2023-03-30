# GSEA_app

A shiny app which takes RNA-seq count matrix and accompanying clinical file, does DE analysis followed by gene-set enrichmenet analysis (GSEA) and plot the results.

```R
library(shiny)
library(DESeq2)
library(tibble)
library(dplyr)
library(fgsea)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)



# UI function
ui <- fluidPage(
  titlePanel("DE Analysis App"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Data Input", 
                 fileInput("rna", "Select RNA file"),
                 fileInput("uromol_clin", "Select clinical data file"),
                 fileInput("pathways", "Select the gene set file (gmt file)"),
                 actionButton("run_analysis", "Run Analysis")
        ),
        tabPanel("DE result",
                 downloadButton("download_results", "Download Results")
        ),
        tabPanel("GSEA result", 
                 downloadButton("downloadGSEA", "Download Results")
        ),
        tabPanel("Plots",
                 plotOutput("plot1"),
                 plotOutput("plot2"),
                 plotOutput("plot3")
        )
      )
    ),
    mainPanel(
      verbatimTextOutput("summary")
    )
  )
)

# Update server function
server <- function(input, output) {
  # Set maximum file size to 300MB
  options(shiny.maxRequestSize = 300 * 1024^2)
  
  # Reactive expression for data input
  data_input <- reactive({
    req(input$rna)
    req(input$uromol_clin)
    req(input$pathways)
    rna <- read.csv(input$rna$datapath, header = T, sep = ",", row.names = 1)
    uromol_clin <- read.csv(input$uromol_clin$datapath, sep = ",", header = T, row.names = 1)
    pathways <- gmtPathways(input$pathways$datapath)
    
    # Ensure that row names/column names in both dataframes are the same 
    uromol_clin <- uromol_clin[colnames(rna),]
    
    # Set grouping variable
    uromol_clin$isTa <- as.factor(ifelse(uromol_clin$Stage == "Ta", "Ta", "non_Ta"))
    
    # Return data frames
    list(rna = rna, uromol_clin = uromol_clin, pathways =pathways)
  })
  
  # Reactive expression for DE analysis
  de_analysis <- reactive({
    req(input$run_analysis)
    rna <- data_input()$rna
    uromol_clin <- data_input()$uromol_clin
    
    # Make DESeqDataSet object
    dds <- DESeqDataSetFromMatrix(countData = rna,
                                  colData = uromol_clin,
                                  design = ~ isTa)
    # Perform DE analysis
    dds <- DESeq(dds)
    
    # Extract results
    res <- data.frame(results(dds))
    #    res <- cbind(as.data.frame(rowData(dds)), res)
    #    res <- res[, -duplicated(colnames(res))]
    #    res <- res %>% 
    #      filter(baseMean > 1) %>%
    #      mutate(log2FoldChange = log2FoldChange,
    #             adj_pval = padj)
    
    res$row <- rownames(res)
    # Map Ensembl gene IDs to symbol. First create a mapping table.
    ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                        key=res$row, 
                                        columns="SYMBOL",
                                        keytype="ENSEMBL")
    names(ens2symbol)[1] <- "row"
    
    ens2symbol <- as_tibble(ens2symbol)
    
    # Joining
    res <- merge(data.frame(res), ens2symbol, by=c("row"))
    # Remove the NAs, averaging statitics for a multi-hit symbols
    res2 <- res %>% 
      dplyr::select(SYMBOL, stat) %>% 
      na.omit() %>% 
      distinct() %>% 
      group_by(SYMBOL) %>% 
      summarize(stat=mean(stat))
    
    # Creating a named vector [ranked genes]
    ranks <- res2$stat
    names(ranks) <- res2$SYMBOL
    #Running fgsea algorithm:
    fgseaRes <- fgseaMultilevel(pathways=data_input()$pathways, stats=ranks)
    
    # Tidy the results:
    fgseaResTidy <<- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES)) # order by normalized enrichment score (NES)
    fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
    #Return data frame
    list(res, fgseaResTidy, fgseaRes, ranks)
  })
  # plot 1 output
  output$plot1 <- renderPlot({
    fgseaResTidy <- de_analysis()[[2]]
    # Add code to generate the plot here
    cols <- c("non-significant" = "grey", "significant" = "red")
    ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title="Hallmark pathways Enrichment Score from GSEA")
    
  })
  #plot2 output
  output$plot2 <- renderPlot({
    ranks <- de_analysis()[[4]]
    pathways <- data_input()$pathways
    fgseaRes <- de_analysis()[[3]]
    
    plotGseaTable(pathways[fgseaRes$pathway[fgseaRes$padj < 0.05]], ranks, fgseaRes, 
                  gseaParam=0.5)
    
  })
  
  #plot3 output
  output$plot3 <- renderPlot({
    # defining variables
    res <- de_analysis()[[1]]
    pathways <- data_input()$pathways
    fgseaResTidy <- de_analysis()[[2]]
    fgseaRes <- de_analysis()[[3]]
    # To see what genes are in each of these pathways:
    gene.in.pathway <- pathways %>% 
      enframe("pathway", "SYMBOL") %>% 
      unnest(cols = c(SYMBOL)) %>% 
      inner_join(res, by="SYMBOL")
    
    # pathways with significant enrichment score
    sig.path <- fgseaResTidy$pathway[fgseaResTidy$adjPvalue == "significant"]
    sig.gen <- unique(na.omit(gene.in.pathway$SYMBOL[gene.in.pathway$pathway %in% sig.path]))
    
    ### create a new data-frame that has '1' for when a gene is part of a term, and '0' when not
    h.dat <- dcast(gene.in.pathway[, c(1,2)], SYMBOL~pathway)
    rownames(h.dat) <- h.dat$SYMBOL
    h.dat <- h.dat[, -1]
    
    h.dat <- h.dat[rownames(h.dat) %in% sig.gen, ]
    h.dat <- h.dat[, colnames(h.dat) %in% sig.path]
    
    # keep those genes with 3  or more occurrences
    #table(data.frame(rowSums(h.dat)))
    
    # 1       2    3    4    5    6 
    # 1604  282   65   11    1    1 
    #h.dat <- h.dat[data.frame(rowSums(h.dat)) >= 3, ]
    num_non_na <- apply(h.dat, 1, function(x) sum(!is.na(x)))
    h.dat <- h.dat[num_non_na >=3, ]
    
    #
    topTable <- res[res$SYMBOL %in% rownames(h.dat), ]
    rownames(topTable) <- topTable$SYMBOL
    
    # match the order of rownames in toptable with that of h.dat
    topTableAligned <- topTable[which(rownames(topTable) %in% rownames(h.dat)),]
    topTableAligned <- topTableAligned[match(rownames(h.dat), rownames(topTableAligned)),]
    #all(rownames(topTableAligned) == rownames(h.dat))
    
    # colour bar for -log10(adjusted p-value) for sig.genes
    dfMinusLog10FDRGenes <- data.frame(-log10(
      topTableAligned[which(rownames(topTableAligned) %in% rownames(h.dat)), 'padj']))
    dfMinusLog10FDRGenes[dfMinusLog10FDRGenes == 'Inf'] <- 0
    
    # colour bar for fold changes for sigGenes
    dfFoldChangeGenes <- data.frame(
      topTableAligned[which(rownames(topTableAligned) %in% rownames(h.dat)), 'log2FoldChange'])
    
    # merge both
    dfGeneAnno <- data.frame(dfMinusLog10FDRGenes, dfFoldChangeGenes)
    colnames(dfGeneAnno) <- c('Gene score', 'Log2FC')
    
    dfGeneAnno[,2] <- ifelse(dfGeneAnno$Log2FC > 0, 'Up-regulated',
                             ifelse(dfGeneAnno$Log2FC < 0, 'Down-regulated', 'Unchanged'))
    colours <- list(
      'Log2FC' = c('Up-regulated' = 'royalblue', 'Down-regulated' = 'yellow'))
    haGenes <- rowAnnotation(
      df = dfGeneAnno,
      col = colours,
      width = unit(1,'cm'),
      annotation_name_side = 'top')
    
    # Now a separate colour bar for the GSEA enrichment padj. This will 
    # also contain the enriched term names via annot_text()
    
    # colour bar for enrichment score from fgsea results
    dfEnrichment <- fgseaRes[, c("pathway", "NES")]
    dfEnrichment <- dfEnrichment[dfEnrichment$pathway %in% colnames(h.dat)]
    dd <- dfEnrichment$pathway
    dfEnrichment <- dfEnrichment[, -1]
    rownames(dfEnrichment) <- dd
    
    colnames(dfEnrichment) <- 'Normalized\n Enrichment score'
    haTerms <- HeatmapAnnotation(
      df = dfEnrichment,
      Term = anno_text(
        colnames(h.dat),
        rot = 45,
        just = 'right',
        gp = gpar(fontsize = 12)),
      annotation_height = unit.c(unit(1, 'cm'), unit(8, 'cm')),
      annotation_name_side = 'left')
    
    # now generate the heatmap
    hmapGSEA <- Heatmap(h.dat,
                        name = 'GSEA hallmark pathways enrichment',
                        split = dfGeneAnno[,2],
                        col = c('0' = 'white', '1' = 'forestgreen'),
                        rect_gp = gpar(col = 'grey85'),
                        cluster_rows = TRUE,
                        show_row_dend = TRUE,
                        row_title = 'Top Genes',
                        row_title_side = 'left',
                        row_title_gp = gpar(fontsize = 11, fontface = 'bold'),
                        row_title_rot = 90,
                        show_row_names = TRUE,
                        row_names_gp = gpar(fontsize = 11, fontface = 'bold'),
                        row_names_side = 'left',
                        row_dend_width = unit(35, 'mm'),
                        cluster_columns = TRUE,
                        show_column_dend = TRUE,
                        column_title = 'Enriched terms',
                        column_title_side = 'top',
                        column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                        column_title_rot = 0,
                        show_column_names = FALSE,
                        show_heatmap_legend = FALSE,
                        clustering_distance_columns = 'euclidean',
                        clustering_method_columns = 'ward.D2',
                        clustering_distance_rows = 'euclidean',
                        clustering_method_rows = 'ward.D2',
                        bottom_annotation = haTerms)
    
  })
  
  # Output DE analysis results
  output$summary <- renderPrint({
    req(input$run_analysis)
    de_analysis()[[2]]
  })
  
  # Download results
  output$download_results <- downloadHandler(
    filename = function() {
      paste("DE_results_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(de_analysis()[[1]], file, row.names = F)
    }
  )
  output$downloadGSEA <- downloadHandler(
    filename = function() {
      paste("fgsea_results", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(apply(de_analysis()[[2]],2,as.character), file, row.names = FALSE)
    }
  )
  
}


# Run the app
shinyApp(ui, server)

```