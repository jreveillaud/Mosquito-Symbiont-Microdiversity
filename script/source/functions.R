# MAKE ITALIC
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))

# PREPARE DATA FOR MED 

load_sequence_file <- function(path_input, filename){
  # COMBINED_SEQUENCES.FASTA FILE
  
  # import fasta file
  fasta <- readDNAStringSet(paste0(path_input, filename))
  
  # sequence headers
  seq_name = names(fasta)
  
  # keep only sample names and number of reads
  seq_name = sub("\\|.*", "", seq_name) 
  
  # sequences
  sequence = paste(fasta)
  
  # create df with new header and sequence
  fasta_df <- data.frame(seq_name, sequence)
  colnames(fasta_df) <- c("Sample", "Seq")
  
  # remove a few r objects 
  rm(seq_name, sequence, fasta)
  
  # create column with only sample name
  fasta_df$Only_sample <- fasta_df$Sample %>% gsub(pattern="_.*", replacement="")
  
  # check the levels of only sample
  samples <- fasta_df$Only_sample %>% as.factor() %>% levels()
  
  # nread by sample
  nread <- c()
  for(i in 1:length(samples)){
    #paste0(samples[i], " : ", fasta_df[fasta_df$Only_sample==samples[i],] %>% nrow(), "\n") %>% cat()
    n <- fasta_df[fasta_df$Only_sample==samples[i],] %>% nrow()
    nread <- c(nread, n)
  }
  
  # count table
  count_table <- data.frame(samples, nread)
  
  # results
  results <- list(count_table, fasta_df)
  return(results)
}


# PREPARE DATA FOR OLIGOTYPING 

prepare_data_for_oligotyping <- function(path_RDATA, path_MED, path_NODES, path_metadata, metadata, path_output, genus){
  
  # Tables preparation
  
  ## Taxonomy
  taxo <- read.table(paste0(path_RDATA, "/1B_MED_nodes_taxo.tsv"))
  
  ## Metadata
  metadata <- read.csv(paste0(path_metadata, metadata), sep=";")

  ## Seqtab

  # files import
  matrix_count <- read.table(paste0(path_MED, "/MATRIX-COUNT.txt"), header = TRUE) %>% t()
  colnames(matrix_count) <- matrix_count[1,]
  matrix_count <- matrix_count[-1,]
  matrix_count <- matrix_count %>% as.data.frame()

  # sequences of representative nodes from MED
  fasta <- readDNAStringSet(paste0(path_MED, "/NODE-REPRESENTATIVES.fasta"))

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

  df <- data.frame()
  for(i in 1:length(nodes_genus)){
    fastaFile <- readDNAStringSet(paste0(path_NODES, nodes_genus[i], ".fa"))
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

  ## Export Wolbachia sequences
  write.fasta(sequences=as.list(df$Seq), names=df$Sample, file.out=paste0(path_output, "2A_oligotyping_", genus, "_sequences.fasta"))
  write.table(count_nodes, paste0(path_output, "2A_MED_nodes_", genus, "_counts.tsv"), row.names=FALSE, quote=FALSE, sep="\t")
}


# COMPUTE STATS ON PHYLOSEQ OBJECTS

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


# RAREFACTION

rarefaction <- function(ps){
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
                      color="Strain", 
                      measures=measures, 
                      nrow = 1)
  
  df <- p1$data
  new_species <- c("AA", "CP", "CQ")
  levels(df$Species) <- new_species
  
  new_order <- c("AA", "CQ", "CP")
  df$Species <- factor(df$Species, levels = new_order)
  
  
  new_location <- c("Field - Guadeloupe", 
                    "Laboratory - Slab TC (Wolbachia -)",
                    "Field - Bosc", 
                    "Field - Camping Europe", 
                    "Laboratory - Lavar")
  df$Strain <- factor(df$Strain, levels = new_location)
  
  p2 <- ggplot(df,aes(Species,value,colour=Strain)) +
    facet_wrap(~ variable, drop=T,scale="free")+
    scale_color_manual("Strain", 
                       values=c("#00BF7D", "#5B6BF7", "#A3A500", "#F8766D", "#00B0F6"),
                       labels = c("Field - Guadeloupe", expression(paste("Laboratory - Slab TC (", italic("Wolbachia"), "-)")), 
                                  "Field - Bosc", "Field - Camping Europe", "Laboratory - Lavar"))+
    geom_boxplot(outlier.colour = NA, alpha=0.8, position = position_dodge2(width=1, preserve="single")) +
    geom_point(size=2, position=position_dodge(width=0.7, preserve='total'))+
    labs(x="Species", y = "Alpha Diversity Measure")+
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


# TRANSFORM TUKEY RESULT INTO DATAFRAME

glht.table <- function(x) {
  # I took this from somewhere, but cant remember the source (probably stackoverflow))
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


# CONTRARY OF %in%
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

percent_taxonomic_plot <- function(dataframe, variable, selection, group, organ, new_names, col, mylabels){
  
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
  df <- read.table(paste0(path_input, "/2B_REF_info_", genus, ".tsv"), sep="\t", header = TRUE)
  
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
  
  MED_alignment_complete <- msa(sequences, order="input")
  MED_alignment_zoom <- msa(sequences, order="input")
  
  setwd(path_output)
  msaPrettyPrint(MED_alignment_complete, showNames="left", showLogo="none", showLogoScale="right", showNumbering="none", askForOverwrite=FALSE,
                 verbose=FALSE, furtherCode=c("\\showruler{1}{top}"), file=paste0("2C_ALIGNMENT_MED_complete_", genus, ".pdf"), 
                 output="pdf")
  
  
  # Alignment with oligotypes only
  df2 <- df[!is.na(df$OLIGO_oligotype_frequency_size),]
  sequences2 <- Biostrings::DNAStringSet(df2$seq)
  names(sequences2) <- df2$ALIGN_name2
  sequences2
  
  
  OLIGO_alignment_complete <- msa(sequences2, order="input")
  OLIGO_alignment_zoom <- msa(sequences2, order="input")
  
  # Export
  msaPrettyPrint(OLIGO_alignment_complete, showNames="left", showLogo="none", showLogoScale="right", showNumbering="none", askForOverwrite=FALSE,
                 verbose=FALSE, furtherCode=c("\\showruler{1}{top}"), file=paste0("2C_ALIGNMENT_oligotype_complete_", genus, ".pdf"), 
                 output="pdf")
  
  msaPrettyPrint(OLIGO_alignment_zoom, showNames="left", showLogo="none", showLogoScale="right", showNumbering="none", askForOverwrite=FALSE,
                 verbose=FALSE, y=oligo_zoom, furtherCode=c("\\showruler{1}{top}"), file=paste0("2C_ALIGNMENT_oligotype_zoom_", genus, ".pdf"),
                 output="pdf")
}

make_reference_table <- function(path_MED, path_OLIGO, path_RDATA, path_TSV, path_OUTPUT, genus){
  # Reference table between MED nodes and size
  
  # import fasta from MED results
  fasta_med <- readDNAStringSet(paste0(path_MED, "/NODE-REPRESENTATIVES.fasta"))
  
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
  MED_count <- read.table(paste0(path_MED, "/MATRIX-COUNT.txt"), header = TRUE) %>% t()
  
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
  taxo <- read.table(paste0(path_TSV, "/1B_MED_nodes_taxo.tsv"))
  
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
  fasta_oligo <- readDNAStringSet(paste0(path_OLIGO, "/OLIGO-REPRESENTATIVES.fasta"))
  
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
  matrix_count <- read.table(paste0(path_OLIGO, "/MATRIX-COUNT.txt"), header = TRUE) %>% t()
  
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
  
  # import ps object
  ps.filter <- readRDS(paste0(path_RDATA, "/1D_MED_phyloseq_decontam.rds"))
  
  # create a long format dataset with abundance of all MED nodes
  data_for_plot_all <- taxo_data(ps.filter, top=ntaxa(ps.filter))
  nodes_all <- levels(data_for_plot_all$Name %>% as.factor())
  nodes_all <- data.frame(nodes_all)
  nodes_all$nodes <- nodes_all$nodes_all
  
  # create a long format dataset with abundance of the 15 most abundant MED nodes
  data_for_plot_15 <- taxo_data(ps.filter, top=15)
  nodes_15 <- levels(data_for_plot_15$Name %>% as.factor())
  nodes_15 <- data.frame(nodes_15)
  nodes_15$nodes <- nodes_15$nodes_15
  
  # create a dataset with abundance of the Other MED nodes (that not corresponding to the 15 most abundant nodes)
  other <- setdiff(nodes_all$nodes_all, nodes_15$nodes_15)
  nodes_other <- data.frame(nodes=1, nodes_other=other)
  nodes_other$nodes <- nodes_other$nodes_other
  
  # merge all the nodes, the 15 most abundant nodes and the other MED nodes
  df <- nodes_all %>% merge(nodes_15, by="nodes", all.x=TRUE)
  df <- df %>% merge(nodes_other, by="nodes", all.x=TRUE)
  df <- df %>% dplyr::select(-c(nodes))
  colnames(df)[1] <- "MED_node_genus"
  
  # add these informations in main table
  main_table <- main_table %>% merge(df, by="MED_node_genus", all.x=TRUE)
  rownames(main_table) <- main_table$MED_node
  colnames(main_table)[12:13] <- c("MED_node_15", "MED_node_other")
  
  # select variables of interest
  main_table <- main_table %>% dplyr::select(c(seq, MED_Genus,MED_node_frequency_size, OLIGO_oligotype_frequency_size, MED_node, oligotype, MED_frequency, OLIGO_frequency, MED_size, OLIGO_size, MED_node_genus, MED_node_15, MED_node_other))
  
  # add a letter system to know if a nodes is among the 15 abudant nodes or not
  main_table$MED_abundance <- main_table$oligotype
  for(i in 1:nrow(main_table)){
    if(!is.na(main_table[i, "MED_node_15"])){
      main_table[i, "MED_abundance"] <- "a"
    } else if(!is.na(main_table[i, "MED_node_other"])){
      main_table[i, "MED_abundance"] <- "o"
    }
  }
  
  # save node and oligotype infos
  write.table(main_table, paste0(path_OUTPUT, "/2B_REF_info_", genus, ".tsv"), sep="\t", quote = FALSE, row.names = FALSE)
}
