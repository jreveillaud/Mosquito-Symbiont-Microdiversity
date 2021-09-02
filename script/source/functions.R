make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))

# PREPARE DATA FOR MED 

load_sequence_file <- function(path_input, filename){
  # FICHIER COMBINED_SEQUENCES.FASTA
  
  # chemin du fichier fasta
  setwd(path_input)
  
  # import du fichier fasta
  #fasta <- readDNAStringSet("combined_sequences.fasta")
  fasta <- readDNAStringSet(filename)
  
  # stocker les headers des séquences
  seq_name = names(fasta)
  
  # garder uniquement les noms d'échantillons et le numéro de read
  seq_name = sub("\\|.*", "", seq_name) 
  
  # stocker les séquences
  sequence = paste(fasta)
  
  # créer un df avec le nouveau header et la séquence
  fasta_df <- data.frame(seq_name, sequence)
  colnames(fasta_df) <- c("Sample", "Seq")
  
  rm(seq_name, sequence, fasta)
  
  fasta_df$Only_sample <- fasta_df$Sample %>% gsub(pattern="_.*", replacement="")
  
  samples <- fasta_df$Only_sample %>% as.factor() %>% levels()
  
  nread <- c()
  for(i in 1:length(samples)){
    #paste0(samples[i], " : ", fasta_df[fasta_df$Only_sample==samples[i],] %>% nrow(), "\n") %>% cat()
    n <- fasta_df[fasta_df$Only_sample==samples[i],] %>% nrow()
    nread <- c(nread, n)
  }
  
  count_table <- data.frame(samples, nread)
  
  results <- list(count_table, fasta_df)
  return(results)
}


# PREPARE DATA FOR OLIGOTYPING 

prepare_data_for_oligotyping <- function(path_RDATA, path_MED, path_NODES, path_metadata, metadata, path_output, genus){
  
  # Préparation des tables
  
  ## Taxonomy
  setwd(path_RDATA)
  taxo <- read.table("1B_MED_nodes_taxo.tsv")
  
  ## Metadata
  metadata <- read.csv(paste0(path_metadata, metadata), sep=";")

  ## Seqtab
  
  # files import
  setwd(path_MED) # matrix count from MED
  matrix_count <- read.table("MATRIX-COUNT.txt", header = TRUE) %>% t()
  colnames(matrix_count) <- matrix_count[1,]
  matrix_count <- matrix_count[-1,]
  matrix_count <- matrix_count %>% as.data.frame()
  
  # sequences of representative nodes from MED
  fasta <- readDNAStringSet("NODE-REPRESENTATIVES.fasta")
  
  # df sequence format
  fasta <- fasta %>% as.data.frame()
  colnames(fasta) <- "seq"
  fasta$nodes <- rownames(fasta)
  
  # remove gaps
  fasta$seq <- fasta$seq %>% gsub(pattern = "-", replacement = "")
  
  # change row and colnames for merge sequences with correspondant nodes
  rownames(matrix_count) <- gsub(x=rownames(matrix_count), pattern="X", replacement="")
  colnames(matrix_count) <- gsub(x=colnames(matrix_count), pattern="Sample-", replacement="")
  fasta$nodes <- gsub(x=fasta$nodes, pattern="\\|.*", replacement="")
  
  # add sequences to matrix count
  matrix_count$nodes <- rownames(matrix_count)
  matrix_count <- matrix_count %>% merge(fasta %>% dplyr::select(c(nodes, seq)), by="nodes")
  rownames(matrix_count) <- matrix_count$seq
  
  
  matrix_count_nodes <- matrix_count
  rownames(matrix_count_nodes) <- matrix_count_nodes$nodes
  
  matrix_count <- matrix_count %>% dplyr::select(-c(seq, nodes))
  matrix_count_nodes <- matrix_count_nodes %>% dplyr::select(-c(seq, nodes))
  
  # reference table nodes/sequences
  taxo$seq <- rownames(taxo)
  taxo <- taxo %>% merge(fasta, by="seq")
  taxo_nodes <- taxo 
  rownames(taxo_nodes) <- taxo_nodes$nodes
  taxo_nodes <- taxo_nodes %>% dplyr::select(-c(seq, nodes))
  
  
  
  ## Select nodes depending on Genus
  
  taxo_nodes_genus <- taxo_nodes[taxo_nodes$Genus==genus,]
  taxo_nodes_genus <- taxo_nodes_genus %>% dplyr::select(-c(Species)) %>% na.omit()
  nodes_genus <- taxo_nodes_genus %>% rownames()
  nodes_genus
  
  setwd(path_NODES)
  df <- data.frame()
  for(i in 1:length(nodes_genus)){
    fastaFile <- readDNAStringSet(paste0(nodes_genus[i], ".fa"))
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)
    node = nodes_genus[i]
    df <- df %>% rbind(data.frame(seq_name, sequence, node))
  }
  
  colnames(df) <- c("Sample", "Seq", "Node")
  
  df %>%
    dplyr::select(-c(Seq)) %>%
    reshape::melt(id.vars="Node") %>%
    dplyr::count(Node) %>%
    dplyr::arrange(-n) -> count_nodes
  
  colnames(count_nodes) <- c("Node", "Count")
  
  ## Export wolbachia sequences
  setwd(path_output)
  write.fasta(sequences=as.list(df$Seq), names=df$Sample, file.out=paste0("oligotyping_", genus, "_sequences.fasta"))
  write.table(count_nodes, paste0("MED_nodes_counts.tsv"), row.names=FALSE, quote=FALSE, sep="\t")
}


# COMPUTE STATS ON PHYLOSEQ OBJECTS

compute_sample_counts <- function(ps){
  sample_counts <- ps %>% sample_sums() %>% length()
  return(sample_counts)
}

compute_taxa_counts <- function(ps){
  taxa_counts <- ps %>% taxa_sums() %>% length()
  return(taxa_counts)
}

compute_read_counts <- function(ps){
  read_counts <- ps %>% taxa_sums() %>% sum()
  return(read_counts)
}

check_ps <- function(ps){
  # check and remove ASV < 1 in each sample
  ps <- prune_taxa(taxa_sums(ps) >= 1, ps)
  
  # check and remove samples with 0 reads
  ps <- prune_samples(sample_sums(ps) >= 1, ps)
}


# DECONTAM

# make_decontam <- function(ps, title, path, threshold){
#   
#   setwd(path)
#   plots <- list()
#   
#   # Preprocess
#   df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
#   df$LibrarySize <- sample_sums(ps)
#   df <- df[order(df$LibrarySize),]
#   df$Index <- seq(nrow(df))
#   
#   p <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=Dna)) + geom_point()
#   plots[[1]] <- p 
#   
#   as.numeric(get_variable(ps, "Dna"))
#   get_variable(ps, "Dna")
#   sample_data(ps)
#   
#   sample_data(ps)$Dna <- as.numeric(get_variable(ps, "Dna"))
#   
#   
#   # Identifier les contaminants - méthode de PREVALENCE
#   sample_data(ps)$is.neg <- sample_data(ps)$Control == "Control sample"
#   contamdf.prev01 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=threshold)
#   
#   table(contamdf.prev01$contaminant)
#   contam.asv.prev01 <- row.names(contamdf.prev01[contamdf.prev01$contaminant == TRUE, ])
#   
#   ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
#   ps.pa.neg <- prune_samples(sample_data(ps.pa)$Control == "Control sample", ps.pa)
#   ps.pa.pos <- prune_samples(sample_data(ps.pa)$Control == "True sample", ps.pa)
#   
#   ## Créer un dataframe de prévalence chez les échantillons positifs et négatifs
#   df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
#                       contaminant=contamdf.prev01$contaminant)
#   
#   p2 <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
#     xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
#   
#   plots[[2]] <- p2
#   
#   ## Ecrire la table d'ASV contaminants détectés par prévalence
#   tax <- tax_table(ps) %>% as.matrix()
#   contam <- tax[row.names(tax) %in% contam.asv.prev01, ]
#   
#   write.table(contam, paste0(title, "_decontam_contaminants.tsv"),
#               sep="\t", quote=F, col.names=NA)
#   
#   # Return
#   output <- list(contam.asv.prev01, contam, plots, contamdf.prev01)
#   return(output)
# }


# CHANGE ORDER OF LEVELS 

change_organ <- function(ps){
  sample_data(ps)$Organ <- factor(sample_data(ps)$Organ, levels=c("Whole", "Pool", "Intestine", "Ovary", "Salivary gland"))
}

change_species <- function(ps){
  sample_data(ps)$Species <- factor(sample_data(ps)$Species, levels=c("Culex pipiens", "Culex quinquefasciatus", "Aedes aegypti"))
}


# RAREFACTION

rarefaction <- function(ps){
  p1 <- ggrare(ps,
               step = 500,
               plot = T,
               parallel = F,
               se = T)
  p1 <- p1 + 
    facet_wrap(~ Organ) + 
    geom_vline(xintercept = min(sample_sums(ps)), 
               color = "gray60") +
    xlim(0,100000) +
    ylim(0, 50) +
    labs(x = "Sample Size", y = "Species Richness")
  
  return(p1)
}


rarefaction2 <- function(ps){
  p1 <- ggrare(ps,
               step = 500,
               plot = T,
               parallel = F,
               se = T)
  p1 <- p1 + 
    facet_wrap(~ Species+Strain+Organ, labeller=label_parsed) + 
    geom_vline(xintercept = min(sample_sums(ps)), 
               color = "gray60") +
    xlim(0,100000) +
    ylim(0, 50) +
    labs(
      x = "Sample Size", y = "Species Richness")
  
  return(p1)
}

# DISTRIBUTION

distribution <- function(ps){
  readsumsdf <- data.frame(nreads = sort(taxa_sums(ps), TRUE),
                           sorted = 1:ntaxa(ps), 
                           type = "OTU")
  
  ggplot(readsumsdf, aes(x = sorted, y = nreads)) + 
    geom_bar(stat = "identity") + 
    scale_y_log10() 
  
  readsumsdf2 <- data.frame(nreads = sort(sample_sums(ps), TRUE), 
                            sorted = 1:nsamples(ps), 
                            type = "Samples")
  
  readsumsdf3 <- rbind(readsumsdf,readsumsdf2)
  
  p  <-  ggplot(readsumsdf3, 
                aes(x = sorted, y = nreads)) + 
    geom_bar(stat = "identity")+
    theme_gray()
  
  p <- p + 
    #ggtitle(paste0("Total number of reads in ", title)) + 
    scale_y_log10() + 
    facet_wrap(~type, 1, scales = "free")
  
  return(p)
}


# ALPHA DIVERSITY

alpha_div <- function(ps, measures){
  data <- ps
  
  p1 <- plot_richness(data, 
                      x="Sample", 
                      color="Location", 
                      measures=measures, 
                      nrow = 1)
  
  df <- p1$data
  new_species <- c("AA", "CP", "CQ")
  levels(df$Species) <- new_species
  
p2 <- ggplot(df,aes(Species,value,colour=Location)) +
    facet_wrap(~ variable, drop=T,scale="free")+
    scale_color_manual(values=c("#00BF7D", "#E76BF3", "#00B0F6", "#F8766D", "#A3A500"))+
    geom_boxplot(outlier.colour = NA, alpha=0.8, position = position_dodge(width=0.9)) +
    geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
    labs(x="Species", y = "Alpha Diversity Measure")

  return(p2)
}

alpha_div2 <- function(ps, measures){
  data <- ps
  
  p1 <- plot_richness(data, 
                      x="Sample", 
                      color="Location", 
                      measures=measures, 
                      nrow = 1)
  
  df <- p1$data
  new_species <- c("AA", "CP", "CQ")
  levels(df$Species) <- new_species
  

  
  p2 <- ggplot(df,aes(Species,value,colour=Location)) +
    facet_grid(variable~Location, drop=T,scale="free")+
    geom_boxplot(outlier.colour = NA, alpha=0.8, position = position_dodge(width=0.9)) +
    geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
    labs(x="Species", y = "Alpha Diversity Measure")
  
  return(p2)
}

alpha_div3 <- function(ps, measures){
  data <- ps
  
  p1 <- plot_richness(data, 
                      x="Sample", 
                      color="Strain", 
                      measures=measures, 
                      nrow = 1)
  
  df <- p1$data
  new_species <- c("AA", "CP", "CQ")
  levels(df$Species) <- new_species
  
  new_order <- c("AA", "CQ", "CP")
  df$Species <- factor(df$Species, levels = new_order)
  
  #print(levels(df$Strain))
  
  #new_location <- c("Guadeloupe", "Slab TC (Wolbachia-)", "Lavar (lab)", "Bosc", "Camping Europe")
  #new_location <- c("Guadeloupe", "Slab TC (Wolbachia-)", "Bosc", "Camping Europe", "Lavar (lab)")
  #new_location <- c("Field", "Guadeloupe", "Bosc", "Camping Europe", "Laboratory", "Slab TC (Wolbachia-)", "Lavar (lab)")
  new_location <- c("Field - Guadeloupe", 
                    "Laboratory - Slab TC (Wolbachia -)",
                    "Field - Bosc", 
                    "Field - Camping Europe", 
                    "Laboratory - Lavar")
  df$Strain <- factor(df$Strain, levels = new_location)
  
  #print(levels(df$Strain))
  
  #col <- c("#00BF7D", "#E76BF3", "#00B0F6", "#F8766D", "#A3A500")
  #col <- c(NA, "#00BF7D", "#F8766D", "#A3A500", NA, "#E76BF3", "#00B0F6")
  
  p2 <- ggplot(df,aes(Species,value,colour=Strain)) +
    facet_wrap(~ variable, drop=T,scale="free")+
    scale_color_manual("Strain", 
                       values=c("#00BF7D", "#5B6BF7", "#A3A500", "#F8766D", "#00B0F6"),
                       labels = c("Field - Guadeloupe", expression(paste("Laboratory - Slab TC (", italic("Wolbachia"), "-)")), 
                                  "Field - Bosc", "Field - Camping Europe", "Laboratory - Lavar"))+
    #geom_boxplot()+
    geom_boxplot(outlier.colour = NA, alpha=0.8, position = position_dodge2(width=1, preserve="single")) +
    #geom_point(size=2,position=position_jitterdodge(dodge.width=0.7)) +
    geom_point(size=2, position=position_dodge(width=0.7, preserve='total'))+
    #theme(legend.key = element_rect(fill = NA))
    labs(x="Species", y = "Alpha Diversity Measure")+
    #scale_x_discrete(expand=c(0.8,0))+
  theme(legend.text.align = 0)
  
  return(p2)
}

# MAKE ANOVA / TUKEY

make_anova_tukey <- function(df, var, value){
  ## ANOVA
  lm1 <- lm(value~var, data=df)
  anova(lm1) -> result_anova 
  
  ## Tukey
  mc_tukey <- glht(lm1, linfct=mcp(var="Tukey"))
  summary(mc_tukey) -> result_tukey
  
  return(list("anova"=result_anova, "tukey"=result_tukey))
}

make_anova <- function(df, var, value, effect){
  lm1 <- lm(value~var*effect, data=df)
  anova(lm1) -> result_anova 
}


# TRANSFORM TUKEY RESULT INTO DATAFRAME

glht.table <- function(x) {
  # I took this from somewehre, but cant remember the source (probably SO))
  pq <- summary(x)$test
  mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
  error <- attr(pq$pvalues, "error")
  pname <- switch(x$alternativ, less = paste("Pr(<", ifelse(x$df ==0, "z", "t"), ")", sep = ""), 
                  greater = paste("Pr(>", ifelse(x$df == 0, "z", "t"), ")", sep = ""), two.sided = paste("Pr(>|",ifelse(x$df == 0, "z", "t"), "|)", sep = ""))
  colnames(mtests) <- c("Estimate", "Std. Error", ifelse(x$df ==0, "z value", "t value"), pname)
  return(mtests)
}


# ADD SIGNIFICANCE TO ANOVA / TUKEY RESULTS

fill_significance <- function(df, var){
  df$Siginificance <- c()
  for(i in 1:nrow(df)){
    x <- df[i, var]
    if(is.na(x)==FALSE){
      if(0 <= x & x <= 0.001){
        a="***"
      }
      else if(0.001 <= x & x <= 0.01){
        a="**"
      }
      else if(0.01 <= x & x <= 0.05){
        a="*"
      }
      else if(0.05 <= x & x <= 0.1){
        a="."
      }
      else if(0.1 <= x & x <= 1){
        a="ns"
      }
    }
    
    else if(is.na(x)==TRUE){
      a=NA
    }
    
    df[i, "Significance"]=a
  }
  return(df)
}

`%ni%` <- Negate(`%in%`)


# TAXONOMIC PLOTS

taxo_data <- function(ps, top){
  top=top

  # convert ps into percent
  ps_global <- transform_sample_counts(ps, function(x) x / sum(x))

  # extract the names of the most abundant taxa depending on Top parameter
  most_abundant_taxa <-  ps_global %>% taxa_sums() %>% sort(TRUE)
  most_abundant_taxa <- most_abundant_taxa[1:top]
  names <- most_abundant_taxa %>% names()

  # melt phyloseq object
  data_for_plot <- psmelt(ps_global)

  # replace name of taxa not including in top by "Other"
  data_for_plot[!data_for_plot$OTU %in% names, "OTU"] <- "Other"

  # add a column with custom name for oligotypes
  data_for_plot$Name <- paste(data_for_plot$OTU,data_for_plot$Genus, sep=".")

  # replace incorrect names by "Other" in the new Name column
  data_for_plot[data_for_plot$OTU=="Other", "Name"] <- "Other"

  # convert Name column to factor
  data_for_plot$Name <- as.factor(data_for_plot$Name)

  # replace Abundance column name to "Relative_Abundance"
  colnames(data_for_plot)[colnames(data_for_plot) %in% "Abundance"] <- "Relative_Abundance"

  # replace name of Species column by "Species.x" to avoid duplication
  colnames(data_for_plot)[colnames(data_for_plot) %in% "sample_Species"] <- "Species.x"

  # return dataframe ready to plot
  return(data_for_plot)
}

# taxo_data <- function(ps, method, top, other){
# 
#   guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 10, face = "italic", colour = "Black", angle = 0)))
# 
#   getPalette = colorRampPalette(brewer.pal(12, "Paired"))
# 
#   top=top
# 
#   if(top>=ntaxa(ps)){
# 
#     if(method=="abundance"){
#       ps_global <- microbiome::transform(ps, "compositional")
# 
#     } else if(method=="count"){
#       ps_global <- ps
#     }
# 
#     df <- as(sample_data(ps_global), "data.frame")
# 
#     most_abundant_taxa = sort(taxa_sums(ps_global), TRUE)[1:top]
# 
#     all_taxa <- taxa_names(ps_global)
#     other_taxa <- all_taxa[!(all_taxa %in% names(most_abundant_taxa))]
# 
#     if(other==FALSE){
# 
#       ps_asv_top = prune_taxa(names(most_abundant_taxa), ps_global)
# 
#       otu_ps_asv_top <- ps_asv_top@otu_table %>% as.data.frame()
#       otu_ps_asv_top <- ps_asv_top@otu_table %>% data.frame()
#       rownames(otu_ps_asv_top) <- paste0("node_", rownames(otu_ps_asv_top))
# 
#       tax_ps_asv_top <- data.frame(ps_asv_top@tax_table)
#       rownames(tax_ps_asv_top) <- paste0("node_", rownames(tax_ps_asv_top))
#       rownames(tax_ps_asv_top) <- rownames(tax_ps_asv_top) %>% as.factor
#     } else if(other==TRUE){
#       ps_other <- prune_taxa(names(other_taxa), ps_global)
# 
#       otu_ps_asv_top <- ps_other@otu_table %>% as.data.frame()
#       otu_ps_asv_top <- ps_other@otu_table %>% data.frame()
#       rownames(otu_ps_asv_top) <- paste0("node_", rownames(otu_ps_asv_top))
# 
#       tax_ps_asv_top <- data.frame(ps_other@tax_table)
#       rownames(tax_ps_asv_top) <- paste0("node_", rownames(tax_ps_asv_top))
#       rownames(tax_ps_asv_top) <- rownames(tax_ps_asv_top) %>% as.factor
#     }
# 
# 
#   } else if(top<ntaxa(ps)){
# 
#     if(method=="abundance"){
#       ps_global <- microbiome::transform(ps, "compositional")
# 
#     } else if(method=="count"){
#       ps_global <- ps
#     }
# 
#     df <- as(sample_data(ps_global), "data.frame")
# 
#     most_abundant_taxa = sort(taxa_sums(ps_global), TRUE)[1:top]
# 
#     all_taxa <- taxa_names(ps_global)
#     other_taxa <- all_taxa[!(all_taxa %in% names(most_abundant_taxa))]
# 
#     if(other==FALSE){
#       ps_asv_top = prune_taxa(names(most_abundant_taxa), ps_global)
# 
#       otu_ps_asv_top <- ps_asv_top@otu_table %>% as.data.frame()
#       tot <- c(1-colSums(otu_ps_asv_top))
# 
#       otu_ps_asv_top <- otu_ps_asv_top %>% rbind(tot)
#       rownames(otu_ps_asv_top)[(top+1)] <- "Other"
# 
#       tax_ps_asv_top <- data.frame(ps_asv_top@tax_table)
#       other <- c("Other", seq(1,ncol(tax_ps_asv_top)))
#       tax_ps_asv_top <- tax_ps_asv_top %>% rbind(other)
#       rownames(tax_ps_asv_top)[(top+1)] <- "Other"
#     } else if(other==TRUE){
#       ps_asv_top = prune_taxa(other_taxa, ps_global)
# 
#       otu_ps_asv_top <- ps_asv_top@otu_table %>% as.data.frame()
# 
#       tax_ps_asv_top <- data.frame(ps_asv_top@tax_table)
#     }
#   }
# 
#   sample_test <- data.frame(ps_asv_top@sam_data)
# 
#   OTU = otu_table(otu_ps_asv_top, taxa_are_rows =TRUE)
#   TAX = tax_table(as.matrix(tax_ps_asv_top))
#   SAM = data.frame(ps_asv_top@sam_data) %>% sample_data
# 
#   ps_asv_top <- phyloseq(OTU, TAX, SAM)
# 
#   taxic_top <- as.data.frame(ps_asv_top@tax_table)
#   taxic_top$ASV <- rownames(taxic_top)
#   taxic_top$ASV <- as.factor(taxic_top$ASV)
# 
#   p <- microbiome::plot_composition(ps_asv_top) +
#     scale_fill_manual(values = getPalette((top)))+
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90)) +
#     ggtitle("Counts") +
#     guide_italics+
#     theme(legend.title = element_text(size = 18), legend.position="none")
# 
#   p
# 
#   data <- p$data
# 
#   colnames(data)[1] <- "ASV"
#   colnames(data)[3] <- "Relative_Abundance"
#   data2 <- data %>% merge(df, by="Sample")
# 
#   data2 <- data2 %>% left_join(taxic_top, by="ASV")
# 
#   if(top>=ntaxa(ps)){
#     data2$ASV <- data2$ASV %>% gsub(pattern = "node_", replacement ="" ) %>% as.factor()
#   }
# 
#   data2$Name <- paste(data2$ASV,data2$Genus, sep=".")
#   data2$Name <- as.factor(data2$Name)
#   levels(data2$Name)[(top+1)] <- "Other"
# 
#   return(data2)
# }

# taxo_data_fast <- function(ps, method){
#   guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 10, face = "italic", colour = "Black", angle = 0)))
#   
#   getPalette = colorRampPalette(brewer.pal(12, "Paired"))
#   
#   ps_global <- microbiome::transform(ps, "compositional")
#   df <- as(sample_data(ps_global), "data.frame")
#   
#   all_taxa <- taxa_names(ps_global)
#   x <- grep(pattern = "oligotype", x=all_taxa) %>% length()
#   most_abundant_taxa <- all_taxa[c(1:x)]
#   other_taxa <- all_taxa[!(all_taxa %in% most_abundant_taxa)]
#   
#   
#   ps_asv_top = prune_taxa(names(most_abundant_taxa), ps_global)
#   
#   otu_ps_asv_top <- ps_asv_top@otu_table %>% as.data.frame()
#   otu_ps_asv_top <- otu_ps_asv_top[c(1:x),]
#   tot <- c(1-colSums(otu_ps_asv_top))
#   
#   otu_ps_asv_top <- otu_ps_asv_top %>% rbind(tot)
#   rownames(otu_ps_asv_top)[(x+1)] <- "Other"
#   
#   tax_ps_asv_top <- ps_asv_top@tax_table %>% as.data.frame()
#   tax_ps_asv_top <- tax_ps_asv_top[c(1:x),] %>% as.data.frame()
#   colnames(tax_ps_asv_top) <- "oligotype"
#   rownames(tax_ps_asv_top) <- tax_ps_asv_top$oligotype
#   
#   other <- c("Other", seq(1,ncol(tax_ps_asv_top)))
#   tax_ps_asv_top <- tax_ps_asv_top %>% rbind(other)
#   rownames(tax_ps_asv_top)[(x+1)] <- "Other"
#   
#   sample_test <- data.frame(ps_asv_top@sam_data)
#   
#   OTU = otu_table(otu_ps_asv_top, taxa_are_rows =TRUE)
#   TAX = tax_table(as.matrix(tax_ps_asv_top))
#   SAM = data.frame(ps_asv_top@sam_data) %>% sample_data
#   
#   ps_asv_top <- phyloseq(OTU, TAX, SAM)
#   
#   taxic_top <- as.data.frame(ps_asv_top@tax_table)
#   taxic_top$ASV <- rownames(taxic_top)
#   taxic_top$ASV <- as.factor(taxic_top$ASV)
#   
#   p <- microbiome::plot_composition(ps_asv_top) + 
#     scale_fill_manual(values = getPalette((x+1)))+
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 90)) +
#     ggtitle("Counts") + 
#     guide_italics+
#     theme(legend.title = element_text(size = 18), legend.position="none")
#   
#   p
#   
#   data <- p$data
#   
#   colnames(data)[1] <- "ASV"
#   colnames(data)[3] <- "Relative_Abundance"
#   data2 <- data %>% merge(df, by="Sample")
#   
#   data2 <- data2 %>% left_join(taxic_top, by="ASV")
#   
#   data2$ASV <- data2$ASV %>% gsub(pattern = "node_", replacement ="" ) %>% as.factor()
#   
#   data2$Name <- paste0(data2$ASV,".")
#   data2$Name <- as.factor(data2$Name)
#   #levels(data2$Name)[(x+1)] <- "Other."
#   
#   return(data2)
# }
# 
# 
# # TAXONOMIC PLOTS (%)
# percent_taxonomic_plot <- function(dataframe, variable, selection, group, organ, new_names, col, mylabels){
# 
#   if(deparse(substitute(variable))=="Species.x"){
#     if(group=="Genus"){
#       test <- dataframe[dataframe$Species.x==selection,] %>%
#         group_by({{variable}}, Name) %>%
#         summarise(rel_ab = mean(Relative_Abundance))
# 
#       test$Genus <- sub(".*\\.", "", test$Name)
#       #test$percent <- test$rel_ab*100
#       test$percent <- (test$rel_ab/sum(test$rel_ab))*100
# 
#       test <- test %>% dplyr::select(c(Species.x, Genus, percent))
# 
#       df <- aggregate(.~Genus+Species.x, test, sum)
#       df <- df %>% arrange(Species.x, percent) %>% group_by({{variable}})
#       df <- df[with(df, order(Species.x,-percent)),]
# 
#       df$Genus <- factor(df$Genus, levels = unique(df$Genus))
#       title="All locations"
# 
#       df$Genus <- factor(df$Genus, levels = new_names)
#       df <- droplevels(df)
# 
#       x=df$Genus
# 
#     } else if(group=="Name"){
#       test <- dataframe[dataframe$Species.x==selection,] %>%
#         group_by({{variable}}, Name) %>%
#         summarise(rel_ab = mean(Relative_Abundance))
# 
#       test$Genus <- sub(".*\\.", "", test$Name)
#       #test$percent <- test$rel_ab*100
#       test$percent <- (test$rel_ab/sum(test$rel_ab))*100
# 
#       test <- test %>% dplyr::select(c(Species.x, Name, percent))
# 
#       df <- aggregate(.~Name+Species.x, test, sum)
#       df <- df %>% arrange(Species.x, percent) %>% group_by({{variable}})
#       df <- df[with(df, order(Species.x,-percent)),]
# 
#       df$Name <- factor(df$Name, levels = unique(df$Name))
#       title="All locations"
# 
#         df$Name <- factor(df$Name, levels=new_names)
#         df <- droplevels(df)
# 
#       x=df$Name
#     }
# 
#   }  else if(deparse(substitute(variable))=="Location"){
#     if(group=="Genus"){
#       test <- dataframe[dataframe$Location==selection,] %>%
#         group_by({{variable}}, Name) %>%
#         summarise(rel_ab = mean(Relative_Abundance))
# 
#       test$Genus <- sub(".*\\.", "", test$Name)
#       #test$percent <- test$rel_ab*100
#       test$percent <- (test$rel_ab/sum(test$rel_ab))*100
# 
#       test <- test %>% dplyr::select(c(Location, Genus, percent))
# 
#       df <- aggregate(.~Genus+Location, test, sum)
#       df <- df %>% arrange(Location, percent) %>% group_by({{variable}})
#       df <- df[with(df, order(Location,-percent)),]
# 
#       df$Genus <- factor(df$Genus, levels = unique(df$Genus))
#       title = paste(levels(df$Species.x)[1], selection, levels(df$Organ)[1])
# 
#       df$Genus <- factor(df$Genus, levels = new_names)
#       df <- droplevels(df)
# 
#       x=df$Genus
# 
#     } else if(group=="Name"){
#       test <- dataframe[dataframe$Location==selection,] %>%
#         group_by({{variable}}, Name) %>%
#         summarise(rel_ab = mean(Relative_Abundance))
# 
#       #test$percent <- test$rel_ab*100
#       test$percent <- (test$rel_ab/sum(test$rel_ab))*100
# 
#       test <- test %>% dplyr::select(c(Location, Name, percent))
# 
#       df <- aggregate(.~Name+Location, test, sum)
#       df <- df %>% arrange(Location, percent) %>% group_by({{variable}})
#       df <- df[with(df, order(Location,-percent)),]
# 
#       df$Name <- factor(df$Name, levels = unique(df$Name))
#       title = paste(levels(df$Species.x)[1], selection, levels(df$Organ)[1])
# 
#       df$Name <- factor(df$Name, levels=new_names)
#         df <- droplevels(df)
# 
#       x=df$Name
#     }
#   }
# 
#   p1 <- ggplot(df, aes(x=x, y=percent, fill = x))+
#     geom_bar(position = "dodge", stat = "identity")+
#     scale_fill_manual(values = col)+
#     theme_bw() +
#     #theme(axis.text.x = element_text(angle = 90)) +
#     theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size=12))+
#     ggtitle("") +
#     guide_italics+
#     theme(legend.title = element_text(size = 18), legend.position="bottom")+
#     theme(panel.spacing=unit(0,"lines"),
#           strip.background=element_rect(color="grey30", fill="grey90"),
#           panel.border=element_rect(color="grey90"),
#           plot.title=element_text(size=10),
#           axis.ticks.x=element_blank()) +
#     geom_text(aes(label=percent %>% round(1)),position=position_dodge(width=0.9), vjust=-0.25, size=4)+
#     scale_y_continuous(breaks=seq(0,100, by=10))+
#     scale_x_discrete(labels = mylabels)+
#     ylim(0,100)+
#     labs(title=title, x=group, y="Relative abundance (%)")
#   p1
# 
#   return(p1)
# }

percent_taxonomic_plot_test <- function(dataframe, variable, selection, group, organ, new_names, col, mylabels){
  
  if(deparse(substitute(variable))=="sample_Species"){
      test <- dataframe[dataframe$sample_Species==selection,] %>% 
        group_by({{variable}}, Genus) %>% 
        summarise(rel_ab = sum(Abundance))
      
      df <- test 
      df$Genus <- factor(df$Genus, levels = unique(df$Genus))
      title="All locations"

      df$Genus <- factor(df$Genus, levels = new_names)
      df <- droplevels(df)

      x=df$Genus
    
  }  else if(deparse(substitute(variable))=="Strain"){
      test <- dataframe[dataframe$Strain==selection,] %>%
        group_by({{variable}}, Genus) %>%
        summarise(rel_ab = sum(Abundance))
      
      df <- test 
      df$Genus <- factor(df$Genus, levels = unique(df$Genus))
      title = paste(levels(df$sample_Species)[1], selection, levels(df$Organ)[1])
      
      df$Genus <- factor(df$Genus, levels = new_names)
      df <- droplevels(df)
      
      x=df$Genus
  }

  p1 <- ggplot(df, aes(x=x, y=rel_ab, fill = x))+
    geom_bar(position = "dodge", stat = "identity")+
    scale_fill_manual(values = col)+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size=12)) +
    #theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1.1, size=15)) +
    ggtitle("") +
    guide_italics+
    theme(legend.title = element_text(size = 18), legend.position="bottom")+
    theme(panel.spacing=unit(0,"lines"),
          strip.background=element_rect(color="grey30", fill="grey90"),
          panel.border=element_rect(color="grey90"),
          plot.title=element_text(size=10),
          axis.ticks.x=element_blank()) +
    geom_text(aes(label=rel_ab),position=position_dodge(width=0.9), size=4, angle=90, hjust=-0.1, vjust=0.25)+
    #geom_text(aes(label=rel_ab),position=position_dodge(width=0.9), size=4, angle=45, hjust=0, vjust=-0.1)+
    ylim(0,2200000)+
    scale_x_discrete(labels = mylabels)+
    labs(title=title, x=group, y="Number of reads")
  p1
  
  return(p1)
}

# MAKE ALIGNMENT
make_alignment <- function(path_input, path_output, genus, oligo_zoom){
  
  # load ref table
  setwd(path_input)
  df <- read.table(paste0("2B_REF_info_", genus, ".tsv"), sep="\t", header = TRUE)
  
  # remove few 0 in MED node names
  df$MED_node <- gsub(df$MED_node, pattern="00000", replacement="")
  
  # create name for alignment plot
  df$ALIGN_name <- paste(df$MED_node, df$oligotype, df$MED_abundance, df$MED_size)
  df$ALIGN_name2 <- paste(df$oligotype, df$MED_node, df$MED_abundance, df$OLIGO_size)
  
  # order df by size of MED nodes
  df <- df[order(df$MED_size, decreasing = TRUE),]
  
  
  # Alignment with all MED nodes and oligotypes
  sequences <- Biostrings::DNAStringSet(df$seq)
  names(sequences) <- df$ALIGN_name
  sequences
  
  setwd(path_output)
  MED_alignment_complete <- msa(sequences, order="input")
  MED_alignment_zoom <- msa(sequences, order="input")
  
  msaPrettyPrint(MED_alignment_complete, showNames="left", showLogo="none", showLogoScale="right", showNumbering="none", askForOverwrite=FALSE,
                 verbose=FALSE, furtherCode=c("\\showruler{1}{top}"), file=paste0("2C_ALIGNMENT_MED_complete_", genus, ".pdf"), 
                 output="pdf")
  
  
  # Alignment with oligotypes only
  df2 <- df[!is.na(df$OLIGO_oligotype_frequency_size),]
  sequences2 <- Biostrings::DNAStringSet(df2$seq)
  names(sequences2) <- df2$ALIGN_name2
  sequences2
  
  setwd(path_output)
  OLIGO_alignment_complete <- msa(sequences2, order="input")
  OLIGO_alignment_zoom <- msa(sequences2, order="input")
  
  msaPrettyPrint(OLIGO_alignment_complete, showNames="left", showLogo="none", showLogoScale="right", showNumbering="none", askForOverwrite=FALSE,
                 verbose=FALSE, furtherCode=c("\\showruler{1}{top}"), file=paste0("2C_ALIGNMENT_oligotype_complete_", genus, ".pdf"), 
                 output="pdf")
  
  msaPrettyPrint(OLIGO_alignment_zoom, showNames="left", showLogo="none", showLogoScale="right", showNumbering="none", askForOverwrite=FALSE,
                 verbose=FALSE, y=oligo_zoom, furtherCode=c("\\showruler{1}{top}"), file=paste0("2C_ALIGNMENT_oligotype_zoom_", genus, ".pdf"),
                 output="pdf")
}

make_reference_table <- function(path_MED, path_OLIGO, path_RDATA, path_PLOT, genus){
  # Reference table between MED nodes and size
  
  # import fasta from MED results
  setwd(path_MED)
  fasta_med <- readDNAStringSet("NODE-REPRESENTATIVES.fasta")
  
  # arrange fasta from MED
  fasta_med <- fasta_med %>% as.data.frame()
  colnames(fasta_med) <- "seq"
  rownames(fasta_med) <- gsub(x=rownames(fasta_med), pattern="00000", replacement="")
  rownames(fasta_med) <- paste0("N", rownames(fasta_med))
  
  # remove gaps
  fasta_med$seq[1]
  fasta_med$seq <- fasta_med$seq %>% gsub(pattern = "-", replacement = "")
  fasta_med$seq[1]
  
  # change name of MED nodes_size
  fasta_med$MED_node_size <- rownames(fasta_med)
  
  #create name and size columns
  fasta_med$MED_node <- fasta_med$MED_node_size %>% gsub(pattern="\\|.*", replace="")
  fasta_med$MED_size <- fasta_med$MED_node_size %>% gsub(pattern=".*\\|size:", replace="")
  
  # replace rownames by name of MED node
  rownames(fasta_med) <- fasta_med$MED_node
  
  # create main_table from fasta_med
  main_table <- fasta_med
  
  
  # Add frequency of MED nodes to the reference table
  
  # import matrix count from MED results
  setwd(path_MED)
  MED_count <- read.table("MATRIX-COUNT.txt", header = TRUE) %>% t()
  
  # arrange MED_count
  colnames(MED_count) <- MED_count[1,]
  MED_count <- MED_count[-1,]
  MED_count <- MED_count %>% as.data.frame()
  
  # change name of MED nodes and remove characters in sample names
  rownames(MED_count) <- gsub(x=rownames(MED_count), pattern="X", replacement="N")
  rownames(MED_count) <- gsub(x=rownames(MED_count), pattern="00000", replacement="")
  colnames(MED_count) <- gsub(x=colnames(MED_count), pattern="Sample-", replacement="")
  
  # extract frequency of MED nodes
  MED_count <- data.frame(sapply(MED_count, function(x) as.numeric(as.character(x))), row.names=rownames(MED_count))
  MED_count[MED_count>0] <- 1 # convert count matrix to presence/absence matrix
  
  # add frequency of MED nodes in MED_count
  MED_count$MED_frequency <- rowSums(MED_count)
  
  # replace rownames by name of MED node
  MED_count$MED_node <- rownames(MED_count)
  
  # add frequency in main table
  main_table <- main_table %>% merge(MED_count %>% dplyr::select(c(MED_node, MED_frequency)), by="MED_node", all.x=TRUE)
  
  # create a summary columun for MED node (name, size, frequency)
  main_table$MED_node_frequency_size <- paste0(main_table$MED_node, " (", main_table$MED_frequency, ") | size:", main_table$MED_size)
  
  
  # Add Genus to the reference table
  
  # import tax table from assignment
  setwd(path_tsv)
  taxo <- read.table(paste0("1B_MED_nodes_taxo.tsv"))
  
  # arrange tax table
  taxo$seq <- rownames(taxo)
  taxo <- taxo %>% dplyr::select(c(seq, everything()))
  colnames(taxo)[7] <- "MED_Genus"
  
  # add taxonomy information in main table
  main_table <- main_table %>% merge(taxo %>% dplyr::select(c(seq, MED_Genus)))
  
  # create a column with name et Genus of MED node
  main_table$MED_node_genus <- paste0(main_table$MED_node, ".", main_table$MED_Genus)
  
  
  # Add oligotypes to the reference table
  
  # import fasta from oligotyping results
  setwd(path_OLIGO)
  fasta_oligo <- readDNAStringSet("OLIGO-REPRESENTATIVES.fasta")
  
  # arrange fasta from oligotyping
  fasta_oligo <- fasta_oligo %>% as.data.frame()
  colnames(fasta_oligo) <- "seq"
  #fasta_oligo$seq
  
  # remove gaps
  fasta_oligo$seq <- fasta_oligo$seq %>% gsub(pattern = "-", replacement = "")
  #fasta_oligo$seq
  
  # create column with name of oligotype
  fasta_oligo$oligotype <- rownames(fasta_oligo)
  
  # add oligotype in main table
  main_table <- main_table %>% merge(fasta_oligo, by="seq", all.x=TRUE)
  
  
  # Add frequency of oligotypes to the reference table
  
  # import matrix count from oligotyping results
  setwd(path_OLIGO)
  matrix_count <- read.table("MATRIX-COUNT.txt", header = TRUE) %>% t()
  
  # arrange matrix count
  colnames(matrix_count) <- matrix_count[1,]
  matrix_count <- matrix_count[-1,]
  matrix_count <- matrix_count %>% as.data.frame()
  
  # extract size and frequency of MED nodes
  matrix_count <- data.frame(sapply(matrix_count, function(x) as.numeric(as.character(x))), row.names=rownames(matrix_count))
  size <- rowSums(matrix_count)
  matrix_count[matrix_count>0] <- 1 # convert count matrix to presence/absence matrix
  frequency <- rowSums(matrix_count)
  
  # add size and frequency in matrix count
  matrix_count$OLIGO_frequency <- frequency
  matrix_count$OLIGO_size <- size
  
  # add column with name of oligotype
  matrix_count$oligotype <- rownames(matrix_count)
  
  # add size and frequency of oligotype in main table
  main_table <- main_table %>% merge(matrix_count %>% dplyr::select(c(oligotype, OLIGO_frequency, OLIGO_size)), by="oligotype", all.x=TRUE, sort=FALSE)
  
  # create column with name and frequency of oligotype
  main_table$OLIGO_oligotype_frequency <- paste0(main_table$oligotype, " (", main_table$OLIGO_frequency,")")
  main_table[main_table$OLIGO_oligotype_frequency=="NA (NA)", "OLIGO_oligotype_frequency"] <- NA
  
  # create column with name, frequency and size of oligotype
  main_table$OLIGO_oligotype_frequency_size <- paste0(main_table$oligotype, " (", main_table$OLIGO_frequency, ") | size:", main_table$OLIGO_size)
  main_table[main_table$OLIGO_oligotype_frequency_size=="NA (NA) | size:NA", "OLIGO_oligotype_frequency_size"] <- NA
  
  # change order of columns in main table
  main_table <- main_table %>% dplyr::select(c(seq, MED_Genus,MED_node_frequency_size, OLIGO_oligotype_frequency_size, MED_node, oligotype, MED_frequency, OLIGO_frequency, MED_size, OLIGO_size, MED_node_genus))
  
  # change rownames and order main table by name of MED node
  rownames(main_table) <- main_table$MED_node
  main_table <- main_table[order(main_table$MED_node),]
  
  
  # Add information about abundance of MED nodes (among the most abundant genus or among Other)
  
  setwd(path_RDATA)
  ps.filter <- readRDS(paste0("1D_MED_phyloseq_decontam.rds"))
  
  # select Culex from phyloseq object
  ps.filter.culex <- subset_samples(ps.filter, Species!="Aedes aegypti" & Organ!="Salivary gland")
  ps.filter.culex <- prune_taxa(taxa_sums(ps.filter.culex) >= 1, ps.filter.culex)
  ps.filter.culex <- prune_samples(sample_sums(ps.filter.culex) >= 1, ps.filter.culex)
  ps.filter.culex
  
  # create a long format dataset with abundance of all MED nodes
  data_for_plot_all <- taxo_data(ps.filter.culex, method = "abundance", top=ntaxa(ps.filter.culex), other=FALSE)
  nodes_all <- levels(data_for_plot_all$Name %>% as.factor())
  nodes_all <- data.frame(nodes_all)
  nodes_all$nodes <- nodes_all$nodes_all
  
  # create a long format dataset with abundance of the 15 most abundant MED nodes
  data_for_plot_15 <- taxo_data(ps.filter.culex, method = "abundance", top=15, other=FALSE)
  nodes_15 <- levels(data_for_plot_15$Name %>% as.factor())
  nodes_15 <- data.frame(nodes_15)
  nodes_15$nodes <- nodes_15$nodes_15
  
  # create a long format dataset with abundance of the Other MED nodes (that not corresponding to the 15 most abundant nodes)
  data_for_plot_15_other <- taxo_data(ps.filter.culex, method = "abundance", top=15, other=TRUE)
  nodes_15_other <- levels(data_for_plot_15_other$Name %>% as.factor())
  nodes_15_other <- data.frame(nodes_15_other)
  nodes_15_other$nodes <- nodes_15_other$nodes_15_other
  
  # merge all, the 15 most abundant and Other MED nodes
  df <- nodes_all %>% merge(nodes_15, by="nodes", all.x=TRUE)
  df <- df %>% merge(nodes_15_other, by="nodes", all.x=TRUE)
  df <- df %>% dplyr::select(-c(nodes))
  colnames(df)[1] <- "MED_node_genus"
  df <- df[-nrow(df),]
  
  # add these informations in main table
  main_table <- main_table %>% merge(df, by="MED_node_genus", all.x=TRUE)
  rownames(main_table) <- main_table$node
  colnames(main_table)[12:13] <- c("MED_node_15", "MED_node_other")
  
  main_table <- main_table %>% dplyr::select(c(seq, MED_Genus,MED_node_frequency_size, OLIGO_oligotype_frequency_size, MED_node, oligotype, MED_frequency, OLIGO_frequency, MED_size, OLIGO_size, MED_node_genus, MED_node_15, MED_node_other))
  
  main_table$MED_abundance <- main_table$oligotype
  
  for(i in 1:nrow(main_table)){
    if(!is.na(main_table[i, "MED_node_15"])){
      main_table[i, "MED_abundance"] <- "a"
    } else if(!is.na(main_table[i, "MED_node_other"])){
      main_table[i, "MED_abundance"] <- "o"
    }
  }
  
  
  # # Distribution
  # 
  # p1 <- ggplot(main_table, aes(x=as.factor(MED_node), y=as.numeric(as.character(MED_size))))+
  #   geom_bar(stat="identity")+
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #   labs(x="Genus", y="Count")
  # 
  # 
  # p2 <- ggplot(main_table, aes(x=as.factor(oligotype), y=as.numeric(as.character(OLIGO_size))))+
  #   geom_bar(stat="identity")+
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #   labs(x="Genus", y="Count")+
  #   geom_text(aes(label=OLIGO_size), position=position_dodge(width=0.9), vjust=-0.25)
  # 
  # 
  # setwd(path_PLOT)
  # 
  # tiff("2B_MED_distribution_plot.tiff", units="in", width=16, height=12, res=300)
  # p1
  # dev.off()
  # 
  # tiff("2B_OLIGO_distribution_plot.tiff", units="in", width=16, height=12, res=300)
  # p2
  # dev.off()
  # 
  # 
  # png("2B_MED_distribution_plot.png", units="in", width=16, height=12, res=300)
  # p1
  # dev.off()
  # 
  # png("2B_OLIGO_distribution_plot.png", units="in", width=16, height=12, res=300)
  # p2
  # dev.off()
  
  
  
  # Save node and oligotype infos
  
  setwd(path_tsv)
  write.table(main_table, paste0("2B_REF_info_", genus, ".tsv"), sep="\t", quote = FALSE, row.names = FALSE)
  #write.table(main_table_oligotype, paste0("REF_oligotype_info_", genus, extension, ".tsv"), sep="\t", quote = FALSE, row.names = FALSE)
}