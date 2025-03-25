# Initial Proteomics Analysis file
# practice working with and analysing proteomics data

# install package from BioConductor BiocManager::install("msqrob2", force = TRUE)

library(adegenet)
library(QFeatures)
library(plotly)
library(msqrob2)
library(ggrepel)
library(limma)
library(lfproQC)
library(tidyverse)
library(readxl)
library(reshape2)
library(stringr)
library(ComplexHeatmap)

# read Excel sheets into a list with named elements then map to data frame ------
dir_path <- paste0(getwd(), "/data/")
re_file <- "combined_annotated_proteins.xlsx"

read_sheets <- function(dir_path, file){
  xlsx_file <- paste0(dir_path, file)
  xlsx_file %>%
    excel_sheets() %>%
    set_names() %>%
    map_df(read_excel, path = xlsx_file, .id = 'sheet_name') %>%
    mutate(file_name = file) %>%
    select(file_name, sheet_name, everything())
}

dat.df <- list.files(dir_path, re_file) %>%
  map_df(~ read_sheets(dir_path, .))

# Add in batch info using which genome sample was used as FASTA file for IDing -----
batch1 <- c('092VI10', '092VI11', '092VI12', '092VI13', '092VI14')
batch2 <- c('092VI15', '092VI16')
batch3 <- c('092VI17', '092VI18')
batches <- data.frame(sheet_name = c(batch1, batch2, batch3), 
                      batch_genome = c(rep('091VI-19', 5), 
                                       rep('091VI-20', 2), 
                                       rep('091VI-21', 2)))

prot.df <- dat.df %>% 
  dplyr::select(-c(1, 55, 4, 5, 10)) %>% 
  filter(Contaminant == FALSE) %>% 
  mutate(eggnog_COG = str_extract(eggNOG_OGs, '^.+?(?=\\@)')) %>% 
  left_join(batches, by = join_by(sheet_name == sheet_name)) %>% 
  dplyr::select(-Contaminant, -eggNOG_OGs)

# reshape data so that abundances by COG not accession are in wide format -------
# if accession number needed, just add it back in where COG is and put COG on end
m.prot <- melt(prot.df, 
               id.vars = c('sheet_name', 'batch_genome', 'eggnog_COG', 
                           'Sequence', 'max_annot_lvl', 'Description_y', 
                           'PFAMs'), 
               measure.vars = 'Abundance: F1: Sample')

# keep missing values as missing for now, next line subs in 0 for them if needed
# put accession at start so formula is accession + batch_genome + eggnog_COG ~
prot.mat <- dcast(data = m.prot, 
                  eggnog_COG + batch_genome ~ sheet_name, 
                  fun.aggregate = sum)
adj.df <- prot.df %>% 
  dplyr::select(-sheet_name) %>% 
  left_join(prot.mat, 
            by = join_by(eggnog_COG == eggnog_COG, 
                         batch_genome == batch_genome))
adj.df <- adj.df[, -c(35:36, 44:46, 28, 3)] # 35:36 may need adjust to 36:37

# create QFeature object --------
coldata <- DataFrame(
  quantCols = seq(1, 9, 1), 
  sample = batches$sheet_name, 
  batch_genome = batches$batch_genome
)

q1 <- readQFeatures(assayData = adj.df, 
                    quantCols = 43:51, 
                    colData = coldata, 
                    name = 'proteins')

head(assay(q1[['proteins']]))
colData(q1)
rowData(q1[['proteins']])$nNonZero <- rowSums(assay(q1[['proteins']]) > 0)
head(q1[[1]])

# assess number of missing values per sample - don't want to delete them just yet
q1 <- zeroIsNA(q1, i = seq_along(q1))
nNA(q1, i = seq_along(q1))

MSnbase::plotNA(assay(q1[['proteins']])) + 
  xlab('Protein index (ordered by data completeness)')

# add log2-transformed data as 2nd assay to not overwrite original
q1 <- addAssay(q1, 
               logTransform(q1[[1]]), 
               name = 'protein_log')

plotDensities(assay(q1[[1]]))
plotDensities(assay(q1[[2]]))

q1 <- filterFeatures(q1, ~ nNonZero >= 2)
nrow(q1[['protein_log']])
# normalise using centered median on log2 transformed data -----
q1 <- addAssay(q1, 
               normalize(q1[['protein_log']], 
                         method = 'center.median'), 
               name = 'protein_norm')

plotDensities(assay(q1[['protein_log']]))
plotDensities(assay(q1[['protein_norm']]))

plotDensities(assay(q1[['protein_log']]), group = colData(q1)$genome_batch)
plotDensities(assay(q1[['protein_norm']]), group = colData(q1)$genome_batch)

# Missing value imputation as methods do not allow for missing values ------
q1.imp <- QFeatures::impute(q1, 
                            method = 'MinProb', 
                            i = 'protein_norm', 
                            name = 'norm_imp')


# Plots and further visuals using normalized data -----
boxplot(assay(q1.imp[['norm_imp']]), col = palette()[-1],
        main = "Peptide distribtutions after normalisation", ylab = "intensity")

plotMA(assay(q1.imp[['norm_imp']]), 
       array = 1)

limma::plotMDS(assay(q1.imp[['norm_imp']]), 
               top = 1000, 
               col = as.numeric(colData(q1.imp)$batch_genome), 
               labels = colnames(assay(q1.imp[['norm_imp']])))

mds_check <- limma::plotMDS(assay(q1.imp[['norm_imp']]), 
                            top = 1000, 
                            plot = FALSE, 
                            col = as.numeric(colData(q1.imp)$batch_genome), 
                            labels = colnames(assay(q1.imp[['norm_imp']])))
mds_check$x

mds_gg <- data.frame(label = rownames(mds_check$distance.matrix.squared), 
                     x = mds_check$x, 
                     y = mds_check$y, 
                     fill = batches$batch_genome)

ggplot(data = mds_gg, aes(x = x, y = y, fill = fill)) + 
  geom_point(size = 3, shape = 21, color = 'black') +
  geom_text_repel(aes(label = label), data = mds_gg, max.overlaps = 10) +
  xlab(paste0("Leading logFoldChange dimension 1 (",
              round(mds_check$var.explained[1] * 100, ), "%)")) +
  ylab(paste0("Leading logFoldChange dimension 2 (",
              round(mds_check$var.explained[2] * 100, ), "%)")) + 
  labs(title = 'Proteome MDS Using log2 Fold Changes Between Samples') + 
  guides(fill = guide_legend(title = 'MI ID WGS Metagenome Group')) + 
  theme(axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom', 
        plot.title.position = 'plot', 
        panel.background = element_rect(fill = 'white'))





#Heatmap(matrix = assay(q1.imp, 'norm_imp'), 
#        show_row_names = FALSE, 
#        clustering_method_rows = 'ward.D2', 
#        clustering_method_columns = 'ward.D2', 
#        row_dend_width = unit(10, 'mm'))





