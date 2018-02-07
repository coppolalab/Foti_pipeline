library(RnBeads)
library(biomaRt)
library(parallel)
library(sva)
library(WGCNA)
library(PMCMR)
library(BayesFactor)
library(limma)

library(Cairo)
library(UpSetR)
library(openxlsx)

library(R.utils)
library(broom)
library(magrittr)
library(stringr)
library(readxl)
library(tidyverse)

#THIS MUST BE FIXED IN PACKAGE
rnb.execute.na.removal.fix <- function (rnb.set, threshold = rnb.getOption("filtering.missing.value.quantile")) {
    if (!inherits(rnb.set, "RnBSet")) {
        stop("invalid value for rnb.set")
    }
    if (!(is.double(threshold) && length(threshold) == 1 && (!is.na(threshold)))) {
        stop("invalid value for threshold")
    }
    if (!(0 <= threshold && threshold <= 1)) {
        stop("invalid value for threshold; expected a value between 0 and 1")
    }
    filterRes <- RnBeads:::rnb.execute.na.removal.internal(rnb.set, NULL,
        threshold)
    list(dataset.before = rnb.set, dataset = remove.sites(rnb.set,
        filterRes$filtered), filtered = filterRes$filtered, threshold = threshold,
        naCounts = filterRes$naCounts)
}

UniqueNames <- function(data_list) {
    comparison_name <- names(data_list)
    df_extract <- data_list[[1]]
    colnames(df_extract)[5] <- str_c(comparison_name, "_log2FC")
    colnames(df_extract)[6] <- str_c(comparison_name, "_p_val")
    colnames(df_extract)[7] <- str_c(comparison_name, "_p_val_adj")
    colnames(df_extract)[8] <- str_c(comparison_name, "_combinedRank")
    list(df_extract)
}

UniqueNames3 <- function(data_list) {
    comparison_name <- names(data_list)
    df_extract <- data_list[[1]]
    colnames(df_extract)[6] <- str_c(comparison_name, "_log2FC")
    colnames(df_extract)[7] <- str_c(comparison_name, "_p_val")
    colnames(df_extract)[8] <- str_c(comparison_name, "_p_val_adj")
    colnames(df_extract)[9] <- str_c(comparison_name, "_combinedRank")
    list(df_extract)
}

UniqueNames2 <- function(data_list) {
    comparison_name <- names(data_list)
    df_extract <- data_list[[1]]
    colnames(df_extract)[6] <- str_c(comparison_name, "_log2FC")
    colnames(df_extract)[7] <- str_c(comparison_name, "_p_val")
    colnames(df_extract)[8] <- str_c(comparison_name, "_p_val_adj")
    colnames(df_extract)[9] <- str_c(comparison_name, "_combinedRank")
    list(df_extract)
}

DMWorkbook <- function(dataset, filename) {
    pval_cols <- colnames(dataset) %>% str_detect("p_val$") %>% which
    adj_pval_cols <- colnames(dataset) %>% str_detect("p_val_adj") %>% which
    logfc_cols <- colnames(dataset) %>% str_detect("log2FC") %>% which
    description_cols <- colnames(dataset) %>% str_detect("Description") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset)
    sig_pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval_cols, rows = 1:nrow(dataset), 
                          rule = "<0.005", style = sig_pvalues)
    conditionalFormatting(wb, 1, cols = adj_pval_cols, rows = 1:nrow(dataset), 
                          rule = "<0.05", style = sig_pvalues)
    conditionalFormatting(wb, 1, cols = logfc_cols, rows = 1:nrow(dataset), 
                          style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = c(1,4), widths = 20)
    setColWidths(wb, 1, cols = 2, widths = "auto")
    setColWidths(wb, 1, cols = 5:ncol(dataset), widths = "auto")
    setColWidths(wb, 1, cols = description_cols, widths = 45)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

DMWorkbook2 <- function(dataset, filename) {
    pval_cols <- colnames(dataset) %>% str_detect("p_val$") %>% which
    adj_pval_cols <- colnames(dataset) %>% str_detect("p_val_adj") %>% which
    logfc_cols <- colnames(dataset) %>% str_detect("log2FC") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset)
    sig_pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval_cols, rows = 1:nrow(dataset), 
                          rule = "<0.005", style = sig_pvalues)
    conditionalFormatting(wb, 1, cols = adj_pval_cols, rows = 1:nrow(dataset), 
                          rule = "<0.05", style = sig_pvalues)
    conditionalFormatting(wb, 1, cols = logfc_cols, rows = 1:nrow(dataset), 
                          style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1:5, widths = 20)
    setColWidths(wb, 1, cols = 6:ncol(dataset), widths = "auto")
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

VolcanoPlot <- function(top_table, filename, plot_name, cutoff_column = "Adj_P_Value", 
                        plot_column = "P_Value", log_column = "logFC", cutoff = 0.05, 
                        xlabel = "Log Fold Change", ylabel = "Log P-value") {
    top_table$Significant <- factor(top_table[[cutoff_column]] < cutoff)
    top_table$Log.Pvalue <- -log10(top_table[[plot_column]])
    p <- ggplot(top_table, aes_string(x = log_column, y = "Log.Pvalue")) + 
        geom_point(aes(color = Significant)) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position = "none", 
              plot.background = element_blank(),
              panel.border = element_rect(size = 1, color = "black"),
              plot.title = element_text(hjust = 0.5)) + 
        xlab(xlabel) + 
        ylab(ylabel) + 
        ggtitle(plot_name) 

    CairoPDF(filename, width = 6, height = 6, bg = "transparent")
    print(p)
    dev.off()
}

MapVolcanoPlots <- function(comparison, prefix, tab_list) {
    tab_reduce <- select(as.data.frame(tab_list[[comparison]]), mean.mean.quot.log2, comb.p.val, comb.p.adj.fdr)
    colnames(tab_reduce) <- c("logFC", "P_Value", "Adj_P_Value")
    VolcanoPlot(tab_reduce, str_c(prefix, comparison, sep = "_"), comparison)
}

PCAPlot <- function(plot_df, color_column, filename){
    p <- ggplot(data = plot_df, aes_string(x = "PC1", y = "PC2", col = color_column)) + 
        geom_point() + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              plot.background = element_blank(), 
              legend.background = element_blank(),
              panel.border = element_rect(color = "black", size = 1)) 

    CairoPDF(filename, width = 8, height = 6, bg = "transparent")
    print(p)
    dev.off()
}

PCAPlot2 <- function(plot_df, label_column, filename, xlims = NULL){
    p <- ggplot(data = plot_df, aes_string(x = "PC1", y = "PC2", label = label_column)) + 
        geom_text() + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
            plot.background = element_blank(), 
              legend.background = element_blank(),
              panel.border = element_rect(color = "black", size = 1)) 
    if (!is.null(xlims)) {
        p <- p + xlim(xlims)
    }

    CairoPDF(filename, width = 8, height = 6, bg = "transparent")
    print(p)
    dev.off()
}

CollapseGenes <- function(gene_df, col_name) {
    unique_vector <- unique(as.character(gene_df[[col_name]]))
    if (length(unique_vector) != 1) {
        return_vector <- str_c(unique_vector, collapse = ",")
    } else {
        return_vector <- unique_vector
    }
    return(return_vector)
}

MakeCompColumn <- function(filter_string, group_vector) {
    filter_match <- str_detect(group_vector, filter_string)
    new_vector <- group_vector
    new_vector[!filter_match] <- NA
    new_vector
}

rnb.options(disk.dump.big.matrices = FALSE) #Set to FALSE because data set is small
rnb.options(enforce.memory.management = TRUE)
rnb.options(logging.disk = TRUE)
rnb.options(import.gender.prediction = FALSE)

source("../../code/common_functions.R") #only needed for parallelized ReadRDSgz and SaveRDSgz

#Read in sample annotation - there must be a column called barcode that matches the array IDs in the raw data exactly
sample_data <- read_csv("../2017-9059-2 Sample Sheet.csv")
sample_data_split <- str_split_fixed(sample_data$`External Sample ID`, " ", 3)
sample_data$Line <- sample_data_split[,1]
sample_data$Cell_Type <- sample_data_split[,2]
sample_data$Replicate <- sample_data_split[,3]
sample_data$barcode <- str_c(sample_data$`chip ID`, sample_data$stripe, sep = "_")
sample_data$Sample_Name <- str_c(sample_data$Cell_Type, sample_data$Replicate, sep = "_")
sample_data_final <- select(sample_data, barcode, Line, Cell_Type, Replicate, Sample_Name)
fixed_samplenames <- filter(sample_data_final, Sample_Name == "_") %>% extract2("Line")
sample_data_final$Sample_Name[sample_data_final$Sample_Name == "_"] <- fixed_samplenames
write_csv(sample_data_final, "../sample_annotation.csv")

#Tell RnBeads where to find directory with IDAT files, sample annotation, and where to write reports
idat_dir <- "../Raw_Data_HD"
sample_annotation <- "../sample_annotation.csv"
report_dir <- "./reports"

#Initialize RnBeads
rnb.initialize.reports(report_dir)
logger.start(fname = NA)
parallel.setup(8)
data_source <- c(idat_dir, sample_annotation)

#Import raw data from IDAT dirs - this process is slow!
result <- rnb.run.import(data.source = data_source, data.type = "infinium.idat.dir", dir.reports = report_dir)
rnb_set <- result$rnb.set
SaveRDSgz(rnb_set, "./save/rnb.set.rda")

#Generate QC report
rnb.run.qc(rnb_set, report_dir)

#For methylation clock - note that the preprocessing steps which remove probes are skipped
rnb_norm_bmiq <- rnb.execute.normalization(rnb_set, method = "bmiq", bgcorr.method = "methylumi.noob") #Run noob background correction and BMIQ normalization - this is very slow!
rnb_narm_bmiq <- rnb.execute.na.removal.fix(rnb_norm_bmiq, 0)$dataset #Drop missing probes
site_bmiq <- meth(rnb_narm_bmiq) #Extract normalized beta values
rownames(site_bmiq) <- rownames(annotation(rnb_narm_bmiq, type = "sites")) #Set rownames to probe IDs

#Remove probes with missing values and replace with rows of NAs - this is only necessary if you want to upload to the webtool.  Otherwise you can skip!
all_probes <- rownames(annotation(rnb_norm_bmiq, type = "sites")) #Get rownames of full probe set
missing_probes <- all_probes[!(all_probes %in% rownames(site_bmiq))] #Find probes which were removed because they were missing
missing_matrix <- matrix(rep(NA, length(missing_probes) * ncol(site_bmiq)), nrow = length(missing_probes), ncol = ncol(site_bmiq)) %>% set_rownames(missing_probes) #make a matrix of NAs corresponding to missing probes
site_complete <- rbind(site_bmiq, missing_matrix) #append row of NAs to other probes
colnames(site_complete) <- sample_data_final$barcode
write.csv(site_complete, "./site_complete.csv") #export normalized data to csv

#Unsuccessful attempt to estimate DNA methylation age
#TransformAge <- function(current_age, adult_age = 20) { 
    #if(current_age < 0) {
        #new_age <- (1 + adult_age) * exp(current_age) - 1
        #return(new_age)
    #} else {
        #new_age <- (1 + adult_age) * current_age + adult_age 
        #return(new_age)
    #}
#} 

#dnam_model <- read_csv("./AdditionalFile3.csv")
#site_clock <- site_bmiq[dnam_model$CpGmarker[-1],]
#predicted_age <- dnam_model$CoefficientTraining[1] + (t(site_clock) %*% dnam_model$CoefficientTraining[-1]) %>%
    #map_dbl(TransformAge)

#DNA methylation age
dnam_import <- read_csv("./site_complete.output.csv")
dnam_import$SampleID %<>% str_replace("^X", "")
dnam_plot <- select(dnam_import, SampleID, DNAmAge) %>% 
    cbind(sample_data_final) %>% as_tibble %>%
    mutate_at(vars(DNAmAge), signif, digits = 3) %>%
    filter(str_detect("G19|G20", Line)) #%>%
    #filter(SampleID != "201496940015_R06C01")
dnam_plot$Sample_Name <- str_c(dnam_plot$Line, dnam_plot$Replicate, sep = "_") %>% 
    factor(levels = c("G19_CD140a+", "G20_CD140a+", "G19_CD44+", "G20_CD44+"))

p <- ggplot(dnam_plot, aes(Sample_Name, DNAmAge, col = Sample_Name)) + 
    geom_point() + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none", 
          plot.background = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(size = 1, color = "black")) + 
    xlab("Cell Type") + 
    ylab("DNA Methylation Age") + 
    scale_color_discrete(name = "Cell Type") 

CairoPDF("dnamage.plot", width = 7, height = 6, bg = "transparent")
print(p)
dev.off()

#Non-parametric ANOVA and post-hoc test
kw_statistic <- kruskal.test(DNAmAge ~ Sample_Name, dnam_plot) %>% tidy
dunn_statistic <- posthoc.kruskal.dunn.test(DNAmAge ~ Sample_Name, dnam_plot, p.adjust.method = "fdr")

#Bayesian ANOVA
anova_age <- lm(DNAmAge ~ Sample_Name, dnam_plot) %>% anova %>% tidy
anova_aov <- aov(DNAmAge ~ Sample_Name, dnam_plot) %>% TukeyHSD %>% tidy

dnam_bf <- anovaBF(DNAmAge ~ Sample_Name, data.frame(dnam_plot)) 
dnam_posterior <- posterior(dnam_bf, iterations = 10000) %>% 
    data.frame %>% select(mu, matches("^Sample_Name")) 
colnames(dnam_posterior) %<>% str_replace_all("Sample_Name\\.", "") %>% str_replace_all("\\.$", "") %>% str_replace_all("\\.", "_")
dnam_posterior_mu <- dnam_posterior$mu
dnam_posterior_predictive <- sweep(select(dnam_posterior, -mu), 1, dnam_posterior_mu, "+")
dnam_posterior_plot <- gather(dnam_posterior_predictive, Sample_Name, Estimate)
dnam_posterior_plot$Sample_Name %<>% factor(levels = c("G19_CD140a", "G20_CD140a", "G19_CD44", "G20_CD44"))

#Compute posterior probability of non-zero difference between contrasts
dnam_posterior_44 <- dnam_posterior$G20_CD44 - dnam_posterior$G19_CD44 
dnam_posterior_44_prob <- which(dnam_posterior_44 < 0) %>% length %>% divide_by(10000)

dnam_posterior_140 <- dnam_posterior$G20_CD140 - dnam_posterior$G19_CD140 
dnam_posterior_140_prob <- which(dnam_posterior_140 < 0) %>% length %>% divide_by(10000)

#Plot posterior estimate of methylation age for each cell type
p <- ggplot(dnam_posterior_plot, aes(x = Sample_Name, y = Estimate, fill = Sample_Name)) + 
    geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) + 
    theme_bw() + 
    theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_rect(color = "black", size = 1),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_blank()) + 
    ylab("Methylation Age")

CairoPDF("estimate_age", width = 7, height = 6, bg = "transparent")
print(p)
dev.off()

#For differential methylation - this includes the full set of filtering steps that remove probes
rnb_filter <- rnb.execute.context.removal(rnb_set)$dataset #Remove non CpG probes
rnb_filter <- rnb.execute.snp.removal(rnb_filter, snp = "any")$dataset #Remove SNPs included for identifying sample swaps
rnb_filter <- rnb.execute.sex.removal(rnb_filter)$dataset #Remove sex chromosomes
SaveRDSgz(rnb_filter, "./save/rnb_filter.rda")

rnb_greedy <- rnb.execute.greedycut(rnb_filter) #Identify probes with low detection score p-value
filter_sites <- rnb_greedy$sites
rnb_filter <- remove.sites(rnb_filter, filter_sites) #Drop probes with low detection score p-value
rnb_filter <- rnb.execute.na.removal.fix(rnb_filter, 0)$dataset #Remove probes with missing values
rnb_filter <- rnb.execute.variability.removal(rnb_filter, 0.005)$dataset #Drop probes with very low variability
rnb_norm <- rnb.execute.normalization(rnb_filter, method = "bmiq", bgcorr.method = "methylumi.noob") #Run noob background correction and BMIQ normalization - this is very slow!
SaveRDSgz(rnb_norm, "./save/rnb_norm.rda")

#Extract beta values and phenotype
promoters_beta <- meth(rnb_norm, type = "promoters")
promoters_annot <- annotation(rnb_norm, type = "promoters") 
promoters_annot$Ensembl_ID <- rownames(promoters_annot)
rownames(promoters_beta) <- rownames(promoters_annot)
write.csv(promoters_beta, "promoters.csv") 
SaveRDSgz(promoters_beta, "./save/promoters_beta.rda")

genes_beta <- meth(rnb_norm, type = "genes")
genes_annot <- annotation(rnb_norm, type = "genes") 
genes_annot$Ensembl_ID <- rownames(genes_annot)
rownames(genes_beta) <- rownames(genes_annot)
write.csv(genes_beta, "genes.csv") 

sites_beta <- meth(rnb_norm, type = "sites")
sites_annot <- annotation(rnb_norm, type = "sites") 
sites_annot$Illumina_ID <- rownames(sites_annot)
rownames(sites_beta) <- rownames(sites_annot)
write.csv(sites_beta, "sites.csv") 

cpgislands_beta <- meth(rnb_norm, type = "cpgislands")
annot_cpgislands_orig <- annotation(rnb_norm, type = "cpgislands")
annot_cpgislands_orig$CpG_Island <- str_c("CpG_", rownames(annot_cpgislands_orig)) 
rownames(cpgislands_beta) <- str_c("CpG_", rownames(annot_cpgislands_orig))
write.csv(cpgislands_beta, "cpgislands.csv") 

pheno_export <- pheno(rnb_norm)
SaveRDSgz(pheno_export, "./save/pheno_export.rda")

##### I have not run the code past this point, so it will not work, at least for the differential methylation step #####

#Principal components analysis
dred_sites <- rnb.execute.dreduction(rnb_norm)
dred_promoters <- rnb.execute.dreduction(rnb_norm, target = "promoters")
dred_genes <- rnb.execute.dreduction(rnb_norm, target = "genes")
dred_cpgislands <- rnb.execute.dreduction(rnb_norm, target = "cpgislands")

#PCA on all sites
dred_sites_plot <- data.frame(dred_sites$mds$manhattan[,1:2])
colnames(dred_sites_plot) <- c("PC1", "PC2")
dred_sites_plot %<>% mutate(Cell_Type = pheno(rnb_norm)$Cell_Type, Sample_Name = pheno(rnb_norm)$Sample_Name)
PCAPlot(dred_sites_plot, "Cell_Type", "sites_pca.pdf")
#PCAPlot2(dred_sites_plot, "Sample_Name", "sites_pca_label.pdf", c(-30000,16000))

#PCA on promoters
dred_promoters_plot <- data.frame(dred_promoters$mds$manhattan[,1:2])
colnames(dred_promoters_plot) <- c("PC1", "PC2")
dred_promoters_plot %<>% mutate(Cell_Type = pheno(rnb_norm)$Cell_Type, Sample_Name = pheno(rnb_norm)$Sample_Name)
PCAPlot(dred_promoters_plot, "Cell_Type", "promoters_pca.pdf")
#PCAPlot2(dred_promoters_plot, "Sample_Name", "promoters_pca_label.pdf", c(-2000,1500))

#PCA on gene bodies
dred_genes_plot <- data.frame(dred_genes$mds$manhattan[,1:2])
colnames(dred_genes_plot) <- c("PC1", "PC2")
dred_genes_plot %<>% mutate(Cell_Type = pheno(rnb_norm)$Cell_Type, Sample_Name = pheno(rnb_norm)$Sample_Name)
PCAPlot(dred_genes_plot, "Cell_Type", "genes_pca.pdf")
#PCAPlot2(dred_genes_plot, "Sample_Name", "genes_pca_label.pdf", c(-1600,1000))

#PCA on CpG Islands
dred_cpgislands_plot <- data.frame(dred_cpgislands$mds$manhattan[,1:2])
colnames(dred_cpgislands_plot) <- c("PC1", "PC2")
dred_cpgislands_plot %<>% mutate(Cell_Type = pheno(rnb_norm)$Cell_Type, Sample_Name = pheno(rnb_norm)$Sample_Name)
PCAPlot(dred_cpgislands_plot, "Cell_Type", "cpgislands_pca.pdf")
#PCAPlot2(dred_cpgislands_plot, "Sample_Name", "cpgislands_pca_label.pdf", c(-1500,1000))

#Annotate CpG Islands
promoters_intervals <- IRanges(start = promoters_annot$Start, end = promoters_annot$End)
promoters_grange <- GRanges(seqnames = promoters_annot$Chromosome, ranges = promoters_intervals, strand = promoters_annot$Strand, mcols = data.frame(Ensembl_ID = rownames(promoters_annot)))
genes_intervals <- IRanges(start = genes_annot$Start, end = genes_annot$End)
genes_grange <- GRanges(seqnames = genes_annot$Chromosome, ranges = genes_intervals, strand = genes_annot$Strand, mcols = data.frame(Ensembl_ID = rownames(genes_annot)))
islands_intervals <- IRanges(start = annot_cpgislands_orig$Start, end = annot_cpgislands_orig$End)
islands_grange <- GRanges(seqnames = annot_cpgislands_orig$Chromosome, ranges = islands_intervals, strand = annot_cpgislands_orig$Strand, mcols = data.frame(CpG_Island = annot_cpgislands_orig$CpG_Island))

merge_genes <- mergeByOverlaps(islands_grange, genes_grange, type = "any")
merge_genes_df <- data.frame(merge_genes) %>% 
    select(mcols.CpG_Island, genes_grange.strand, mcols.Ensembl_ID)
colnames(merge_genes_df) <- c("CpG_Island", "Strand", "Ensembl_ID")
merge_genes_df$Region <- "Gene_Body"

merge_promoters <- mergeByOverlaps(islands_grange, promoters_grange, type = "any")
merge_promoters_df <- data.frame(merge_promoters) %>% 
    select(mcols.CpG_Island, promoters_grange.strand, mcols.Ensembl_ID)
colnames(merge_promoters_df) <- c("CpG_Island", "Strand", "Ensembl_ID")
merge_promoters_df$Region <- "Promoter"
merge_all_df <- rbind(merge_genes_df, merge_promoters_df)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
cpgislands_bm_table <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = as.character(unique(merge_all_df$Ensembl_ID)), mart = ensembl)
colnames(cpgislands_bm_table) <- c("Ensembl_ID", "Symbol")

merge_all_df %<>% left_join(cpgislands_bm_table)
merge_all_df$Symbol[is.na(merge_all_df$Symbol)] <- "Unannotated"
merge_all_df$Symbol[nchar(merge_all_df$Symbol) == 0] <- "Unannotated"

merge_all_df$Combined <- str_c(merge_all_df$CpG_Island, ";", merge_all_df$Region, ":", merge_all_df$Strand)
merge_all_collapse_ensembl <- split(merge_all_df, merge_all_df$Combined) %>% 
    map_chr(CollapseGenes, "Ensembl_ID")
merge_all_collapse_ensembl_df <- tibble(Combined = names(merge_all_collapse_ensembl), all = merge_all_collapse_ensembl)
merge_all_collapse_ensembl_df$CpG_Island <- str_split_fixed(merge_all_collapse_ensembl_df$Combined, ";", 2) %>% extract(TRUE,1)
merge_all_collapse_ensembl_df$Location <- str_split_fixed(merge_all_collapse_ensembl_df$Combined, ";", 2) %>% 
    extract(TRUE,2) %>% str_c("Ensembl_ID_", .)
merge_all_spread_ensembl <- select(merge_all_collapse_ensembl_df, -Combined) %>% 
    spread(Location, all)

merge_all_collapse <- split(merge_all_df, merge_all_df$Combined) %>% 
    map_chr(CollapseGenes, "Symbol")
merge_all_collapse_df <- tibble(Combined = names(merge_all_collapse), all = merge_all_collapse)
merge_all_collapse_df$CpG_Island <- str_split_fixed(merge_all_collapse_df$Combined, ";", 2) %>% extract(TRUE,1)
merge_all_collapse_df$Location <- str_split_fixed(merge_all_collapse_df$Combined, ";", 2) %>% extract(TRUE,2)
merge_all_spread <- select(merge_all_collapse_df, -Combined) %>% 
    spread(Location, all)

merge_all_final <- left_join(merge_all_spread, merge_all_spread_ensembl) %>% left_join(annot_cpgislands_orig)
write.xlsx(merge_all_final, "cpgislands_annot.xlsx")

#Gene Annotation
promoters_bm_table <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description', 'gene_biotype'), filters = 'ensembl_gene_id', values = rownames(promoters_annot), mart = ensembl)
promoters_bm_table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(promoters_bm_table) <- c("Ensembl_ID", "Symbol", "Description", "Gene_Type")

genes_bm_table <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol','description', 'gene_biotype'), filters = 'ensembl_gene_id', values = rownames(genes_annot), mart = ensembl)
genes_bm_table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(genes_bm_table) <- c("Ensembl_ID", "Symbol", "Description", "Gene_Type")

#Differential methylation
contrasts_list <- list(NPC_vs_ES = "NPC|ES",
                       Astro_DIV180_vs_ES = "ASTRO_DIV180|ES",
                       Astro_DIV180_vs_NPC = "ASTRO_DIV180|NPC",
                       Astro_DIV210_vs_Astro_DIV180 = "ASTRO_DIV210|ASTRO_DIV180",
                       Astro_DIV240_vs_Astro_DIV210 = "ASTRO_DIV240|ASTRO_DIV210",
                       Astro_DIV240_vs_Astro_DIV180 = "ASTRO_DIV240|ASTRO_DIV180",
                       OPC_DIV180_vs_ES = "OPC_DIV180|ES",
                       OPC_DIV180_vs_NPC = "OPC_DIV180|NPC",
                       OPC_DIV210_vs_OPC_DIV180 = "OPC_DIV210|OPC_DIV180",
                       OPC_DIV240_vs_OPC_DIV210 = "OPC_DIV240|OPC_DIV210",
                       OPC_DIV240_vs_OPC_DIV180 = "OPC_DIV240|OPC_DIV180",
                       Astro_DIV180_vs_OPC_DIV180 = "ASTRO_DIV180|OPC_DIV180",
                       Astro_DIV210_vs_OPC_DIV210 = "ASTRO_DIV210|OPC_DIV210",
                       Astro_DIV240_vs_OPC_DIV240 = "ASTRO_DIV240|OPC_DIV240")

cell_type_reformat <- str_replace_all(pheno(rnb_norm)$Cell_Type, "\\(", "_") %>% 
    str_replace_all("\\)", "") %>%
    factor(levels = c("ES", "NPC", "OPC_DIV180", "OPC_DIV210", "OPC_DIV240", "ASTRO_DIV180", "ASTRO_DIV210", "ASTRO_DIV240"))

contrast_cols <- map(contrasts_list, MakeCompColumn, cell_type_reformat) %>% reduce(data.frame) %>% set_colnames(names(contrasts_list))

rnb_norm <- addPheno(rnb_norm, contrast_cols$NPC_vs_ES, "NPC_vs_ES") %>%
    addPheno(contrast_cols$Astro_DIV180_vs_ES, "Astro_DIV180_vs_ES") %>%
    addPheno(contrast_cols$Astro_DIV180_vs_NPC, "Astro_DIV180_vs_NPC") %>%
    addPheno(contrast_cols$Astro_DIV210_vs_Astro_DIV180, "Astro_DIV210_vs_Astro_DIV180") %>%
    addPheno(contrast_cols$Astro_DIV240_vs_Astro_DIV210, "Astro_DIV240_vs_Astro_DIV210") %>%
    addPheno(contrast_cols$Astro_DIV240_vs_Astro_DIV180, "Astro_DIV240_vs_Astro_DIV180") %>%
    addPheno(contrast_cols$OPC_DIV180_vs_ES, "OPC_DIV180_vs_ES") %>%
    addPheno(contrast_cols$OPC_DIV180_vs_NPC, "OPC_DIV180_vs_NPC") %>%
    addPheno(contrast_cols$OPC_DIV210_vs_OPC_DIV180, "OPC_DIV210_vs_OPC_DIV180") %>%
    addPheno(contrast_cols$OPC_DIV240_vs_OPC_DIV210, "OPC_DIV240_vs_OPC_DIV210") %>%
    addPheno(contrast_cols$OPC_DIV240_vs_OPC_DIV180, "OPC_DIV240_vs_OPC_DIV180") %>%
    addPheno(contrast_cols$Astro_DIV180_vs_OPC_DIV180, "Astro_DIV180_vs_OPC_DIV180") %>%
    addPheno(contrast_cols$Astro_DIV210_vs_OPC_DIV210, "Astro_DIV210_vs_OPC_DIV210") %>%
    addPheno(contrast_cols$Astro_DIV240_vs_OPC_DIV240, "Astro_DIV240_vs_OPC_DIV240") 

comp_cols <- colnames(contrast_cols)
reg_types <- c("genes", "promoters", "cpgislands")

diffmeth <- rnb.execute.computeDiffMeth(rnb_norm, pheno.cols = comp_cols, region.types = reg_types)
SaveRDSgz(diffmeth.adj, "./save/diffmeth.adj.rda")

comparisons <- get.comparisons(diffmeth)
comparisons_format <- str_replace(comparisons, " \\(.*$", "") %>% 
    str_replace("\\.", "") %>%
    str_replace_all(" ", "_")

tab_sites <- map(comparisons, get.table, object = diffmeth, region.type = "sites", return.data.frame = TRUE) 
names(tab_sites) <- comparisons_format
tab_promoters <- map(comparisons, get.table, object = diffmeth, region.type = "promoters", return.data.frame = TRUE) 
names(tab_promoters) <- comparisons_format
tab_genes <- map(comparisons, get.table, object = diffmeth, region.type = "genes", return.data.frame = TRUE) 
names(tab_genes) <- comparisons_format
tab_cpgislands <- map(comparisons, get.table, object = diffmeth, region.type = "cpgislands", return.data.frame = TRUE) 
names(tab_cpgislands) <- comparisons_format

SaveRDSgz(tab_sites, "./save/tab_sites.rda")

map(comparisons_format, MapVolcanoPlots, "promoters", tab_promoters)
map(comparisons_format, MapVolcanoPlots, "genes", tab_genes)
#map(comparisons_format, MapVolcanoPlots, "sites", tab_sites)
map(comparisons_format, MapVolcanoPlots, "cpgislands", tab_cpgislands)

#Add annotation to DE tables
tab_promoters_annot <- map(tab_promoters, mutate, Ensembl_ID = rownames(promoters_annot)) %>% 
    map(join, promoters_annot) %>%
    map(join, promoters_bm_table) %>%
    map(dplyr::select, Ensembl_ID, Symbol:Gene_Type, mean.mean.quot.log2, comb.p.val, comb.p.adj.fdr, combinedRank, Chromosome:Strand, num.sites, CpG:G, entrezID, mean.mean.g1:mean.mean.diff, mean.num.na.g1:mean.nsamples.covg.thresh.g2) 
SaveRDSgz(tab_promoters_annot, "./save/tab_promoters_annot.rda")

tab_genes_annot <- map(tab_genes, mutate, Ensembl_ID = rownames(genes_annot)) %>% 
    map(join, genes_annot) %>%
    map(join, genes_bm_table) %>%
    map(dplyr::select, Ensembl_ID, Symbol:Gene_Type, mean.mean.quot.log2, comb.p.val, comb.p.adj.fdr, combinedRank, Chromosome:Strand, num.sites, CpG:G, entrezID, mean.mean.g1:mean.mean.diff, mean.num.na.g1:mean.nsamples.covg.thresh.g2) 
SaveRDSgz(tab_genes_annot, "./save/tab_genes_annot.rda")

tab_sites_annot <- map(tab_sites, mutate, Illumina_ID = sites_annot$Illumina_ID) %>% 
    map(join, sites_annot) %>%
    map(dplyr::select, Illumina_ID, Chromosome:Strand, mean.quot.log2, diffmeth.p.val, diffmeth.p.adj.fdr, combinedRank, `CGI Relation`:GC, AddressA:`Mismatches B`, `SNPs 3`:`Cross-reactive`, mean.g1:mean.diff, max.g1:min.diff, num.na.g1:covg.thresh.nsamples.g2) 
SaveRDSgz(tab_sites_annot, "./save/tab_sites_annot.rda")

tab_cpgislands_annot <- map(tab_cpgislands, mutate, CpG_Island = str_c("CpG", rownames(annot_cpgislands_orig), sep = "_")) %>% 
    map(join, merge_all_final) %>%
    map(select, CpG_Island, `Gene_Body:-`:`Promoter:+`, mean.mean.quot.log2, comb.p.val, comb.p.adj.fdr, combinedRank, Chromosome:Strand, num.sites, `Ensembl_ID_Gene_Body:-`:`Ensembl_ID_Promoter:+`, CpG:G, mean.mean.g1:mean.mean.diff, mean.num.na.g1:mean.nsamples.covg.thresh.g2) 
SaveRDSgz(tab_cpgislands_annot, "./save/tab_cpgislands_annot.rda")

names(tab_promoters_annot) <- comparisons_format
names(tab_genes_annot) <- comparisons_format
#names(tab_sites_annot) <- comparisons_format
names(tab_cpgislands_annot) <- comparisons_format

##Subset tables for export
promoters_table_reduce <- map(tab_promoters_annot, dplyr::select, Ensembl_ID:G) %>% 
    lmap(UniqueNames) %>% 
    reduce(join) %>%
    dplyr::select(Ensembl_ID:Gene_Type, dplyr::contains("log2FC"), 
                  dplyr::matches("p_val$"), dplyr::contains("p_val_adj"), 
                  dplyr::contains("combinedRank"), Chromosome:G) 
promoters_table_reduce[,grepl("p_val|log2FC", colnames(promoters_table_reduce))] %<>% signif(3)
DMWorkbook(promoters_table_reduce, "./promoters.xlsx")

genes_table_reduce <- map(tab_genes_annot, dplyr::select, Ensembl_ID:G) %>% 
    lmap(UniqueNames) %>% 
    reduce(join) %>%
    dplyr::select(Ensembl_ID:Gene_Type, dplyr::contains("log2FC"), 
                  dplyr::matches("p_val$"), dplyr::contains("p_val_adj"), 
                  dplyr::contains("combinedRank"), Chromosome:G) 
genes_table_reduce[,grepl("p_val|log2FC", colnames(genes_table_reduce))] %<>% signif(3)
DMWorkbook(genes_table_reduce, "./genes.xlsx")

sites_table_reduce <- map(tab_sites_annot, dplyr::select, Illumina_ID:Context) %>%
    lmap(UniqueNames3) %>% 
    reduce(join) %>%
    dplyr::select(Illumina_ID:Strand, dplyr::contains("log2FC"), dplyr::matches("p_val$"), dplyr::contains("p_val_adj"), dplyr::contains("combinedRank"), `CGI Relation`:Context) 
sites_table_reduce[,grepl("p_val|log2FC", colnames(sites_table_reduce))] %<>% signif(3)
SaveRDSgz(sites_table_reduce, "./save/sites_table_reduce.rda")
write.csv(sites_table_reduce, "./sites_combined.csv", row.names = FALSE)

cpgislands_table_reduce <- map(tab_cpgislands_annot, dplyr::select, CpG_Island:`Ensembl_ID_Promoter:+`) %>% 
    lmap(UniqueNames2) %>% 
    reduce(left_join) %>%
    dplyr::select(CpG_Island:`Promoter:+`, dplyr::contains("log2FC"), 
                  dplyr::matches("p_val$"), dplyr::contains("p_val_adj"), 
                  dplyr::contains("combinedRank"), Chromosome:`Ensembl_ID_Promoter:+`) 
cpgislands_table_reduce[,grepl("p_val|log2FC", colnames(cpgislands_table_reduce))] %<>% signif(3)
DMWorkbook2(cpgislands_table_reduce, "./cpgislands.xlsx")

#Upset
comp_upset <- str_subset(colnames(promoters_table_reduce), "p_val_adj")[c(1:4,6:8,13:14)] 
comp_filters <- str_c(comp_upset, " < 0.05")
comp_format <- str_replace_all(comp_upset, "_p_val_adj", "")

top_upset_promoters <- select(promoters_table_reduce, Ensembl_ID, dplyr::contains("p_val_adj")) 
top_list_promoters <- map(comp_filters, filter_, .data = top_upset_promoters) %>% map(extract2, "Ensembl_ID") 
names(top_list_promoters) <- comp_format
top_list_promoters_set1 <- top_list_promoters[c(1:3,6:7)]
top_list_promoters_set2 <- top_list_promoters[c(4:5,8:9)]

top_upset_genes <- select(genes_table_reduce, Ensembl_ID, dplyr::contains("p_val_adj")) 
top_list_genes <- map(comp_filters, filter_, .data = top_upset_genes) %>% map(extract2, "Ensembl_ID") 
names(top_list_genes) <- comp_format
top_list_genes_set1 <- top_list_genes[c(1:3,6:7)]
top_list_genes_set2 <- top_list_genes[c(4:5,8:9)]

top_upset_cpgislands <- select(cpgislands_table_reduce, CpG_Island, dplyr::contains("p_val_adj")) 
top_list_cpgislands <- map(comp_filters, filter_, .data = top_upset_cpgislands) %>% map(extract2, "CpG_Island") 
names(top_list_cpgislands) <- comp_format
top_list_cpgislands_set1 <- top_list_cpgislands[c(1:3,6:7)]
top_list_cpgislands_set2 <- top_list_cpgislands[c(4:5,8:9)]

CairoPDF("upset_promoters_es_npc", width = 12, height = 6, bg = "transparent")
upset(fromList(top_list_promoters_set1), nsets = 5, nintersects = NA, order.by = "freq")
dev.off()

CairoPDF("upset_promoters_opc_astro", width = 10, height = 6, bg = "transparent")
upset(fromList(top_list_promoters_set2), nsets = 4, nintersects = NA, order.by = "freq")
dev.off()

CairoPDF("upset_genes_es_npc", width = 12, height = 6, bg = "transparent")
upset(fromList(top_list_genes_set1), nsets = 5, nintersects = NA, order.by = "freq")
dev.off()

CairoPDF("upset_genes_opc_astro", width = 10, height = 6, bg = "transparent")
upset(fromList(top_list_genes_set2), nsets = 4, nintersects = NA, order.by = "freq")
dev.off()

CairoPDF("upset_cpgislands_es_npc", width = 12, height = 6, bg = "transparent")
upset(fromList(top_list_cpgislands_set1), nsets = 5, nintersects = NA, order.by = "freq")
dev.off()

CairoPDF("upset_cpgislands_opc_astro", width = 10, height = 6, bg = "transparent")
upset(fromList(top_list_cpgislands_set2), nsets = 4, nintersects = NA, order.by = "freq")
dev.off()

