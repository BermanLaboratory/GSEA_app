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
                 plotOutput("myPlot")
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
    
    #Return data frame
    list(res, fgseaResTidy)
  })
  
  output$myPlot <- renderPlot({
    
    # Add code to generate the plot here
    fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
    cols <- c("non-significant" = "grey", "significant" = "red")
    ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title="Hallmark pathways Enrichment Score from GSEA")
    
  })
  
  # Output DE analysis results
  output$summary <- renderPrint({
    req(input$run_analysis)
    de_analysis()
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