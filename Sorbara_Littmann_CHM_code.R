## Data Analysis and figures generated in R for Sorbara and Littmann et al; Cell Host and Microbe 2020


# SETUP: Load Required Packages -------------------------------------------


library(RPostgreSQL)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape)
library(seqinr)
library(gridExtra)
library(muscle)
library(msa)
library(ape)
library(ggtree)
library(gridExtra)
library(rPython)
library(umap)
library(phangorn)
library(ggplot2)
library(dplyr)
library(reshape)
library(treeio)
library(viridis)
library(ComplexHeatmap)
library(circlize)

# SUMMARY: Final Figure Plots -----------------------------------------------------
#Main text figures
  
  #16S rRNA Phylogenetic trees:
  Figure_2
  Figure_2_inset

  #Core Genome
  Figure_3A_plot
  #Figure 3B: values given below
  Figure_3D
  
  #Core Genome Phylogenetic Tree
  Figure_4_concat 

  #UMAP Analysis
  Figure_5A
  Figure_5B
  Figure_5C
  grid.arrange(Figure_5B,Figure_5C,ncol=2)
  
  #16S rRNA distance verse UMAP distance.
  Figure_6_A_Boot
  Figure_6B <- do.call(grid.arrange,c(individual_species_distances[c(17,7,19)],ncol=3))
  
  #Intraspecies differences
  Figure_7_A
  Figure_7A <- do.call(grid.arrange,c(next_lachno,ncol=4))
  Figure_7B_table
  Figure_7C_table
  Figure_7D_Top
  Figure_7D_Bottom
  

#Supplemental Figures
  #Figure S1
  Figure_S1_b
  Figure_s1_b_inset
  Figure_S1_c
  grid.arrange(Figure_S1_b,Figure_S1_c,ncol=1)
  
  #Figure S2
  do.call(grid.arrange,c(Figure_s2_left,ncol=1))
  do.call(grid.arrange,c(Figure_s2_right,ncol=1))
  
  #Figure S3
  figure_S3A_plot
  Figure_S3_b_plot
  figure_s3_c
  figure_s3_d_plot
  Figure_S3_e
  Figure_S3_f
  grid.arrange(Figure_S3_e,Figure_S3_f,ncol=2)
  
  #Figure S4
  Figure_S4_A
  Figure_S4_B
  Figure_S4_acetate
  Figure_S4_acetoacetate
  Figure_S4_butyrate_kinase
  
  #Figure S5
  Figure_S5A <- do.call(grid.arrange,c(lachno_ind_annotate,ncol=4))
  Figure_S5b <- do.call(grid.arrange,c(lachno_ind_hypo,ncol=4))
  #Figure S6
  Figure_S6
  

#Reviewer Figures
  Figure_1_reviewers
  Figure_1b_reviewers
  grid.arrange(Figure_1_reviewers,Figure_1b_reviewers,ncol=1)
  
  Figure_r_3
  Figure_r_3b
  
# SETUP: Define functions and color palettes --------------------------------------------------------


'%!in%' <- function(x,y)!('%in%'(x,y))

colors_species <- read.csv("chm_species_colors.csv",header=TRUE,sep=",")
colors_species$species_name <- as.character(colors_species$species_name)
colors_species$color <- as.character(colors_species$color)
species.colors <- colors_species$color
names(species.colors) <- colors_species$species_name
species.colors
#https://gist.github.com/shaunpwilkinson/2c4ded3c99a3fe8a08d01c8352bac012
read.aa <- function(file = file.choose(), bin = FALSE){
  x <- readLines(file)
  namelines <- grepl("^>", x)
  f <- cumsum(namelines)
  res <- split(x, f)
  resnames <- sapply(res, function(s) s[1])
  resnames <- gsub("^>", "", resnames)
  resnames <- gsub("\\|.+", "", resnames) 
  res <- lapply(res, function(s) paste0(s[-1], collapse = ""))
  names(res) <- resnames
  if(bin){
    res <- lapply(res, charToRaw)
    class(res) <- "AAbin"
  }
  return(res)
}
# SETUP: Load RAW data --------------------------------------------------------------------

  #1. Annotations of Full Biobank
  load("kraken_library_972_int.RData")
  
  #2. List of Lachnospiraceae isolates included in analysis
  lachnospiraceae_tax <- read.csv("lachnospiraceae_taxonomy.csv",header=TRUE,sep=",")
  
  #3.  List of Lachnospiraceae isolates that have been sequenced more than once
  exclude_dup_seq <- read.csv("seq_id_exclude_duplicates.csv",header=TRUE,sep=",")

  #4 Lachnospiraceae pH
  lachno_pH <- read.csv("lachnospiraceae_pH_final.csv",header=TRUE,sep=",") 
  
  #5 type_strain info.
  lachno_type_strains <- read.csv("16S_type_strains.csv",header=TRUE,sep=",")
  lachno_type_strains$type <- "type_strain"

  
  #6 Butyrate pathway annotations
  
  butyrate_genes <- read.csv("butyrate_pathway_annotations.csv",stringsAsFactors = F) %>%
    mutate(prokka_target=gsub("^s","S",prokka_target))
  butyrate_genes$prokka_target <- iconv(butyrate_genes$prokka_target, from = 'UTF-8', to = 'ASCII//TRANSLIT')
  
  #7 Kegg lookup table
  
  con <- dbConnect(dbDriver("PostgreSQL"),host="10.151.15.23",dbname="dfi_commensal_library",user="dfi_lab",password="dfilab")
  kegg_uniprot_lookup <- tbl(con,"kegg_uniprot_lookup")
  kegg_uniprot_lookup <- as.data.frame(kegg_uniprot_lookup)
  
  #8 Centroid of hypothetical clusters identification
  matt_clusters_refseq <- tbl(con,"matt_clusters_refseq_lookup")
  matt_clusters_refseq <- as.data.frame(matt_clusters_refseq)
  
  
  
  #9 top hits for Blast Identified Clusters.
  tophits <- read.csv(file="hypothetical_protein_clusters_refseq.csv",header=TRUE,sep=",")
  tophits$qseqid <- as.character(tophits$qseqid)
  tophits$sseqid <- as.character(tophits$sseqid)  
  tophits$name <- as.character(tophits$name)
  
  #10 KEGG info.
  load("~/Dropbox/Matt's Cloudtop/Postdoc Cloudtop/Klebsiella/MS DATA/R/data_to_upload/prokka_kegg_lookup.RData")
  load("~/Dropbox/Matt's Cloudtop/Postdoc Cloudtop/Klebsiella/MS DATA/R/data_to_upload/locus_uniprot_lookup.RData")
  
  #11 R gnavus glucorhamnan loci coverage
  polysaccharide_groups <- read.csv(file="coverage_groups.csv",header=TRUE,sep=",")
  polysaccharide_coverage <- read.csv(file="polysaccharide_coverage.csv",header=TRUE,sep=",")
# SETUP: Working_lists_of_isolates -----------------------------------------------


exclude_dup_seq$seq_id_exlcude <- as.character(exclude_dup_seq$seq_id_exlcude)

lachnospiraceae_working <- lachnospiraceae_tax %>%
  inner_join(lookup,by="msk_id") %>%
  filter(is_lachno==1,
         seq_id %!in% exclude_dup_seq$seq_id_exlcude)

lachnospiraceae_ruminococcaceae <- lachnospiraceae_tax %>%
  inner_join(lookup,by="msk_id") %>%
  filter(is_lachno ==1 | is_lachno ==2,
         seq_id %!in% exclude_dup_seq$seq_id_exlcude)

lachnospiraceae_clostridiaceae <- lachnospiraceae_tax %>%
  inner_join(lookup,by="msk_id") %>%
  filter(is_lachno ==1 | is_lachno == 3,
         seq_id %!in% exclude_dup_seq$seq_id_exlcude)


# Figure 2, Figure S1: 16S_rRNA_Phylogenetic_trees ---------------------------------------------

#generate list of 16S locus tags, using longest 16S rRNA sequence for each isolate
rRNA_genes_lachno_rumino <- prokka_genes %>% 
  filter (seq_id %in% lachnospiraceae_ruminococcaceae$seq_id,prokka_genes$ftype == 'rRNA')

rRNA_genes_lachno_clostrid <- prokka_genes %>% 
  filter (seq_id %in% lachnospiraceae_clostridiaceae$seq_id,prokka_genes$ftype == 'rRNA')

rRNA_16S_lachno_rumino <- rRNA_genes_lachno_rumino[grep("16S",rRNA_genes_lachno_rumino$product),] %>%
  group_by(seq_id)  %>%
  arrange(-length_bp) %>%
  dplyr::slice(1)
rRNA_16S_lachno_clostrid <- rRNA_genes_lachno_clostrid[grep("16S",rRNA_genes_lachno_clostrid$product),] %>%
  group_by(seq_id)  %>%
  arrange(-length_bp) %>%
  dplyr::slice(1)

# Retrieve nucleotide sequences and export fasta files (provided)
lachno_rumino_to_align <- prokka_seq %>% 
  filter(locus_tag %in% rRNA_16S_lachno_rumino$locus_tag) %>%
  inner_join(rRNA_16S_lachno_rumino,by="locus_tag") %>%
  inner_join(lachnospiraceae_ruminococcaceae,by="seq_id")

lachno_clostrid_to_align <- prokka_seq %>% 
  filter(locus_tag %in% rRNA_16S_lachno_clostrid$locus_tag) %>%
  inner_join(rRNA_16S_lachno_clostrid,by="locus_tag") %>%
  inner_join(lachnospiraceae_clostridiaceae,by="seq_id")

#for (i in 1:(nrow(lachno_rumino_to_align))){
#  write.fasta(lachno_rumino_to_align[i,]$nuc_sequence,lachno_rumino_to_align[i,]$msk_id,file.out="lachno_rumino_16S_rRNA.fasta",open="a")
#}

#for (i in 1:(nrow(lachno_clostrid_to_align))){
#  write.fasta(lachno_clostrid_to_align[i,]$nuc_sequence,lachno_clostrid_to_align[i,]$msk_id,file.out="lachno_clostrid_16S_rRNA.fasta",open="a")
#}


#Generate MUSCLE Alignments:
  # run muscle: muscle -in lachno_rumino_16S_rRNA.fasta -out lachno_rumino_16S.afa
  # run muscle: muscle -in lachno_clostrid_16S_rRNA.fasta -out lachno_clostrid_16S.afa

  # convert muscle .afa file to phylip.
  
  # generate tree using PhyML 3.0: hosted at: http://www.atgc-montpellier.fr/phyml/
    # automatic model selection using SMS
    # starting tree with BioNJ
    # NNI nearest neighbor interchange
    # aLRT-SH Fast likelihood method.


#FIGURE 2 TREE Visualization:
  lachno_rumino_tree <- read.tree(file="lachno_rumino_tree.nwk") %>%
     root("MSK.23.74")
 
  

  Figure_2 <- ggtree(lachno_rumino_tree,layout='circular')+#geom_label(aes(label=node), size=3)+
    geom_hilight(node=397,alpha=0.4,fill="#8AD188",extendto=0.37)+
    geom_hilight(node=280,alpha=0.4,fill="#8C855E",extendto=0.37)+
    geom_hilight(node=302,alpha=0.4,fill="#638E5E",extendto=0.37)+
    geom_hilight(node=402,alpha=0.4,fill="#CEF948",extendto=0.37)+
    geom_hilight(node=471,alpha=0.4,fill="#F6A6FF",extendto=0.37)+
    geom_hilight(node=443,alpha=0.4,fill="#ED3580",extendto=0.37)+ #LUTI
    geom_hilight(node=447,alpha=0.4,fill="#ED3580",extendto=0.37)+
    geom_hilight(node=451,alpha=0.4,fill="#ED3580",extendto=0.37)+
    geom_hilight(node=442,alpha=0.4,fill="#ED3580",extendto=0.37)+
    geom_hilight(node=174,alpha=0.4,fill="#ED3580",extendto=0.37)+
    geom_hilight(node=154,alpha=0.4,fill="#ED3580",extendto=0.37)+
    geom_hilight(node=469,alpha=0.4,fill="#A79AFF",extendto=0.37)+ #glucersea
    geom_hilight(node=194,alpha=0.4,fill="#A79AFF",extendto=0.37)+ #glucersea
    geom_hilight(node=324,alpha=0.4,fill="#DAAD86",extendto=0.37)+ #Selmin.
    geom_hilight(node=354,alpha=0.4,fill="#FFF5BA",extendto=0.37)+ #Ty nexil.
    geom_hilight(node=534,alpha=0.4,fill="#254A9E",extendto=0.37)+ #C. clostridif
    geom_hilight(node=530,alpha=0.4,fill="#2E9DB7",extendto=0.37)+ #C. symbosium
    geom_hilight(node=268,alpha=0.4,fill="#242851",extendto=0.37)+ #C. celer
    geom_hilight(node=391,alpha=0.4,fill="#A44DF7",extendto=0.37)+ #F. faecal
    geom_hilight(node=533,alpha=0.4,fill="#AFF8DB",extendto=0.37)+ #C. aldenense
    geom_hilight(node=427,alpha=0.4,fill="#B5B5B5",extendto=0.37)+ #B. producta
    geom_hilight(node=512,alpha=0.4,fill="#DCD3FF",extendto=0.37)+ #F. saccharivorans
    geom_hilight(node=195,alpha=0.4,fill="#4C4C4C",extendto=0.37)+ #B. caecimuris
    geom_hilight(node=150,alpha=0.4,fill="#8221F4",extendto=0.37)+ #B. hansenii
    geom_hilight(node=125,alpha=0.4,fill="#F2B15C",extendto=0.37)+ #R. intestinalis
    geom_hilight(node=330,alpha=0.4,fill="#383838",extendto=0.37)+ #Unclassified
    geom_hilight(node=466,alpha=0.4,fill="#9E2556",extendto=0.37)+#B. faecis
    geom_hilight(node=187,alpha=0.4,fill="#9E2556",extendto=0.37)+#B. faecis
    geom_hilight(node=436,alpha=0.4,fill="#FFB5E8",extendto=0.37)+ #B. schinkii
    geom_hilight(node=338,alpha=0.4,fill="#6EB5FF",extendto=0.37)+
    geom_hilight(node=331,alpha=0.4,fill="#ACE7FF",extendto=0.37)+ #D. formicigenerans
    geom_hilight(node=355,alpha=0.4,fill="#F4B571",extendto=0.37)+
    geom_hilight(node=336,alpha=0.4,fill="#877B5C",extendto=0.37)+ #C. scindens
    geom_hilight(node=455,alpha=0.4,fill="#FFCCF9",extendto=0.37)+
    geom_treescale()
  
  Figure_2 <- Figure_2 %>% rotate(335) %>% rotate(510)%>% rotate(425) %>% rotate(465) 
  
  for(i in 1:nrow(Figure_2$data)){
    ifelse(Figure_2$data$branch.length[i]>1,Figure_2$data$branch.length[i]<-.35,Figure_2$data$branch.length[i])
    ifelse(Figure_2$data$x[i]>1,Figure_2$data$x[i]<-.35,Figure_2$data$x[i])
  }


  figure_2_tip_labels <- lachno_rumino_to_align %>% 
      select(msk_id,Donor,species)
  figure_2_tip_labels$seq <- figure_2_tip_labels$msk_id
  figure_2_tip_labels <- figure_2_tip_labels %>%
      select(seq,species)
  figure_2_tip_labels$species <- as.character(figure_2_tip_labels$species)

  Figure_2 <- Figure_2 %<+% figure_2_tip_labels + 
    geom_tippoint(aes(fill=species),size=2,pch=21,color="black",alpha=0.7)+
    scale_fill_manual(values=species.colors)+
    theme(legend.position = "none")
  Figure_2
  
#Figure 2 Inset:
  Figure_2_inset <- ggtree(lachno_rumino_tree)
  Figure_2_inset <- Figure_2_inset %>% rotate(335) %>% rotate(510)%>% rotate(425) %>% rotate(465) 
  Figure_2_inset

#Figure 1a Reviewers - Unrooted Tree
  Figure_1_reviewers <- ggtree(lachno_rumino_tree,layout="equal_angle")
  Figure_1_reviewers <- Figure_1_reviewers %<+% figure_2_tip_labels + 
    geom_tippoint(aes(fill=species),size=2,pch=21,color="black",alpha=0.7)+
    scale_fill_manual(values=species.colors)+
    theme(legend.position = "none")

  #Figure S1b

  lachno_clostrid_tree <- read.tree(file="lachno_clostrid_tree.nwk") %>%
    root("MSK.15.2")

  Figure_S1_b <- ggtree(lachno_clostrid_tree,layout="circular")+
    geom_hilight(node=453,alpha=0.4,fill="#F6A6FF",extendto=0.37)+ #Wexlerae
    geom_hilight(node=448,alpha=0.4,fill="#9E2556",extendto=0.37)+ #faecis
    geom_hilight(node=169,alpha=0.4,fill="#9E2556",extendto=0.37)+ #faecis
    geom_hilight(node=177,alpha=0.4,fill="#4C4C4C",extendto=0.37)+ #caecimuris
    geom_hilight(node=451,alpha=0.4,fill="#A79AFF",extendto=0.37)+ #glucerase
    geom_hilight(node=176,alpha=0.4,fill="#A79AFF",extendto=0.37)+ #glucerase
    geom_hilight(node=437,alpha=0.4,fill="#FFCCF9",extendto=0.37)+ #obeum
    geom_hilight(node=420,alpha=0.4,fill="#FFB5E8",extendto=0.37)+ #schinkii
    geom_hilight(node=426,alpha=0.4,fill="#ED3580",extendto=0.37)+ #luti
    geom_hilight(node=427,alpha=0.4,fill="#ED3580",extendto=0.37)+#luti
    geom_hilight(node=433,alpha=0.4,fill="#ED3580",extendto=0.37)+#luti
    geom_hilight(node=413,alpha=0.4,fill="#ED3580",extendto=0.37)+#luti
    geom_hilight(node=409,alpha=0.4,fill="#B5B5B5",extendto=0.37)+#producta
    geom_hilight(node=132,alpha=0.4,fill="#8221F4",extendto=0.37)+#hansenii
    geom_hilight(node=494,alpha=0.4,fill="#DCD3FF",extendto=0.37)+#hansenii
    geom_hilight(node=515,alpha=0.4,fill="#AFF8DB",extendto=0.37)+#aldense
    geom_hilight(node=516,alpha=0.4,fill="#254A9E",extendto=0.37)+#clostridoforme
    geom_hilight(node=512,alpha=0.4,fill="#2E9DB7",extendto=0.37)+#symbo
    geom_hilight(node=250,alpha=0.4,fill="#242851",extendto=0.37)+#celercres
    geom_hilight(node=525,alpha=0.4,fill="#CEF948",extendto=0.37)+#rectale
    geom_hilight(node=400,alpha=0.4,fill="#A44DF7",extendto=0.37)+#faecalicatena
    geom_hilight(node=276,alpha=0.4,fill="#8C855E",extendto=0.37)+#hadrus
    geom_hilight(node=304,alpha=0.4,fill="#8AD188",extendto=0.37)+#eutactus
    geom_hilight(node=349,alpha=0.4,fill="#6EB5FF",extendto=0.37)+#longica
    geom_hilight(node=347,alpha=0.4,fill="#877B5C",extendto=0.37)+#scin
    geom_hilight(node=342,alpha=0.4,fill="#ACE7FF",extendto=0.37)+#formici
    geom_hilight(node=341,alpha=0.4,fill="#383838",extendto=0.37)+#unclass
    geom_hilight(node=335,alpha=0.4,fill="#DAAD86",extendto=0.37)+#selimon
    geom_hilight(node=313,alpha=0.4,fill="#638E5E",extendto=0.37)+#comes
    geom_hilight(node=311,alpha=0.4,fill="#FFF5BA",extendto=0.37)+#nexilis
    geom_hilight(node=365,alpha=0.4,fill="#F4B571",extendto=0.37)+#gnavus
    geom_hilight(node=251,alpha=0.4,fill="#F2B15C",extendto=0.37)+#rose
    geom_treescale()
  #Shorten C. cochlearium branch
  
  for(i in 1:nrow(Figure_S1_b$data)){
    ifelse(Figure_S1_b$data$branch.length[i]>1,Figure_S1_b$data$branch.length[i]<-.35,Figure_S1_b$data$branch.length[i])
    ifelse(Figure_S1_b$data$x[i]>1,Figure_S1_b$data$x[i]<-.35,Figure_S1_b$data$x[i])
  }
  Figure_S1_b<- Figure_S1_b %>% rotate(346)
  
  figure_s1_b_tip_labels <- lachno_clostrid_to_align %>% 
    select(msk_id,Donor,species)
  figure_s1_b_tip_labels$seq <-  figure_s1_b_tip_labels$msk_id
  figure_s1_b_tip_labels <-  figure_s1_b_tip_labels %>%
    select(seq,species)
  figure_s1_b_tip_labels$species <- as.character( figure_s1_b_tip_labels$species)
  
  Figure_S1_b <- Figure_S1_b %<+% figure_s1_b_tip_labels + 
    geom_tippoint(aes(fill=species),size=2,pch=21,color="black",alpha=0.7)+
    scale_fill_manual(values=species.colors)+
    theme(legend.position = "none")
  Figure_S1_b
#Figure S1_b inset:
  Figure_s1_b_inset <- ggtree(lachno_clostrid_tree) %>%
    rotate (346)
  Figure_s1_b_inset  
  
#Figure 1b Reviewers
  
  Figure_1b_reviewers <- ggtree(lachno_clostrid_tree,layout="equal_angle")
  Figure_1b_reviewers <- Figure_1b_reviewers %<+% figure_s1_b_tip_labels + 
    geom_tippoint(aes(fill=species),size=2,pch=21,color="black",alpha=0.7)+
    scale_fill_manual(values=species.colors)+
    theme(legend.position = "none")
  
#Figure S1 B: 16S tree with type strains
 #for (i in 1:(nrow(lachno_rumino_to_align))){
#    write.fasta(lachno_rumino_to_align[i,]$nuc_sequence,lachno_rumino_to_align[i,]$msk_id,file.out="figure_s1.fasta",open="a")
#  }

  # combine with fasta file of 16S_type_strains.fasta
  # align with muscle: muscle -in figure_s1_fasta -out figure_s1_b.afa
  # convert align fasta file to phylip format
  # analyze with PhyML 3.0
  # output tree in newick format: 16S_type_strains_tree.tre
  
  type_strain_16S_tree <- ape::read.tree (file="16S_type_strains_tree.tre",) %>%
    root("MSK.23.74")
  
  
  Figure_S1_c <- ggtree(type_strain_16S_tree,layout="circular")+
    geom_treescale()
  
  Figure_S1_c <- Figure_S1_c %>% 
    rotate(435) %>% rotate(369) %>% rotate(490) %>% rotate(528)
  
  #shorten the F. plautii branch
  for(i in 1:nrow(Figure_S1_c$data)){
    ifelse(Figure_S1_c$data$branch.length[i]>1,Figure_S1_c$data$branch.length[i]<-.35,Figure_S1_c$data$branch.length[i])
    ifelse(Figure_S1_c$data$x[i]>1,Figure_S1_c$data$x[i]<-.35,Figure_S1_c$data$x[i])
  }
  
  #Create dataframe for tip_labels
  type_strains_tip_labels <- lachno_type_strains %>%
    select(fasta_code,species_name,type)
  colnames(type_strains_tip_labels) <- c("seq","species","type")
  
  Figure_S1_c_tip_labels <- lachno_rumino_to_align %>%
    select(msk_id,species) %>%
    dplyr::mutate(seq=msk_id)%>%
    select(seq,species) %>%
    dplyr::mutate(type ="Isolate") %>%
    dplyr::mutate(species=as.character(species)) %>%
    rbind(type_strains_tip_labels) %>%
    dplyr::mutate(species=as.character(species))
  
  
  Figure_S1_c <- Figure_S1_c %<+% Figure_S1_c_tip_labels +
    geom_tippoint(aes(fill=species,shape=type,size=type),color="black",alpha=0.7)+
    scale_size_manual(values=c(1,3))+
    scale_shape_manual(values=c(21,22))+
    scale_fill_manual(values=species.colors)+
    theme(legend.position = "none")
  Figure_S1_c


# Figure S2 16S rRNA intra-species diversity, percent identities of hits ------------------------------------------------

  lachno_16S_distances <- cophenetic.phylo(lachno_rumino_tree)
  lachno_16S_distances[lower.tri(lachno_16S_distances)]<- NA
  lachno_16S_distances <- as.data.frame(lachno_16S_distances)
  lachno_16S_distances <- lachno_16S_distances[row.names(lachno_16S_distances)!='MSK.23.74',colnames(lachno_16S_distances)!='MSK.23.74']
  
  lachno_16S_distances$id <- rownames(lachno_16S_distances)
  
  melted_16S_distance <- melt(lachno_16S_distances)
  melted_16S_distance$variable <- as.character(melted_16S_distance$variable)
  melted_16S_distance <- melted_16S_distance %>% 
    filter(value!='NA') %>%
    filter(variable != id) %>%
    merge(lachnospiraceae_working,by.x='id',by.y='msk_id') %>%
    select(id,variable,value,species) %>%
    merge(lachnospiraceae_working,by.x="variable",by.y="msk_id") %>%
    select(id,variable,value,species.x,species.y)
  
  melted_16S_distance$species.x <- as.character(melted_16S_distance$species.x)
  melted_16S_distance$species.y <- as.character(melted_16S_distance$species.y)
  
  number_per_species <- lachnospiraceae_working %>%
    group_by(species) %>%
    tally() %>%
    arrange(desc(n))
  
  
  num_species <- length(unique(lachnospiraceae_working$species))
  next_species <- unique(number_per_species$species)
  next_species
  average_distances_by_species <- data.frame(species=character(),average_distance=numeric(),number_in_species=integer())
  plot_distances <- NULL
  x <- 0
  
  for(i in 1:num_species){
    print(i)
    
    intra_species_16S <- melted_16S_distance %>% 
      filter(species.x== next_species[i] | species.y== next_species[i] ) %>%
      filter(species.x==species.y)
    
    data_intra_species_16S <- data.frame(species=next_species[i],average_distance=mean(intra_species_16S$value),number_in_species=nrow(intra_species_16S))
    average_distances_by_species <- rbind(average_distances_by_species,data_intra_species_16S)
    
    
    x <- ifelse(nrow(intra_species_16S)>1,x+1,x)
    
   ifelse(nrow(intra_species_16S)>1,
           plot_distances[[x]] <- ggplot(intra_species_16S,aes(x=value))+
               geom_density(aes(y=..scaled..,fill=species.x,color=species.x),alpha=.7)+
              scale_fill_manual(values=species.colors)+
              scale_color_manual(values=species.colors)+
               scale_y_continuous(breaks=c(0,1),limits=c(0,1))+
            theme(legend.position = "none",axis.text.x = element_text(size=5),axis.text.y = element_text(size=5), axis.title.y = element_blank(),axis.title.x = element_blank(),panel.grid.major = element_blank(),panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),panel.grid.minor = element_blank() )+
            xlim(0,0.08),0)
    
  }
  
  Figure_s2_left <- plot_distances
 
  Figure_s2_right <- NULL
  for(i in 1:num_species){
    print(i)
    intra_species_df <- lachnospiraceae_working %>% 
      filter(species==next_species[i])
    Figure_s2_right[[i]] <- ggplot(intra_species_df,aes(x=species,y=pident,fill=species))+
      geom_bar(stat="summary",fun.y="mean",alpha=0.7)+
      geom_jitter(position = position_jitter(height = 0, width = .3),size=0.5,alpha=0.7)+
      scale_fill_manual(values=species.colors)+
      theme(legend.position = "none", 
            axis.text.y = element_blank(),
            axis.title=element_blank(),
            axis.text.x = element_text(size=5),
            panel.grid.major = element_blank(),
            panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
            panel.grid.minor = element_blank())+
      coord_flip(ylim=c(93.5,100))
  }
# Figure S3 Generation/Identification of Hypothetical Protein Clusters, Types of CDS Annotation per isolate----------------------------------------------------------------
 
   #select all loci from Lachnospiraceae
  lachno_genes <- prokka_genes %>%
    filter(seq_id %in% lachnospiraceae_working$seq_id) %>%
    filter(ftype=='CDS')
  glimpse(lachno_genes)
  lachno_genes$count <- 1
  isolate_summary <- lachno_genes %>%
    group_by(seq_id) %>%
    dplyr::summarise(sum(count))
  
  colnames(isolate_summary) <- c("seq_id","num_genes")
  median(isolate_summary$num_genes)
  max(isolate_summary$num_genes)
  min(isolate_summary$num_genes)
  
  #retrieve sequences of all hypothetical proteins in lachnospiraceae
  lachno_genes_hypothetical_seqs <- prokka_seq %>% 
    filter(locus_tag %in% lachno_genes$locus_tag) %>%
    filter(prot_sequence != '') %>%
    filter(product=='hypothetical protein') %>%
    inner_join(lachno_genes,by='locus_tag')
  
  # determine protein length, and order in descending order
  lachno_genes_hypothetical_seqs$length <- nchar(lachno_genes_hypothetical_seqs$prot_sequence)
  lachno_genes_hypothetical_seqs <- lachno_genes_hypothetical_seqs[order(-lachno_genes_hypothetical_seqs$length),]
  
  #output a FASTA file
  for (i in 1:nrow(lachno_genes_hypothetical_seqs)){
    write.fasta(sequences=lachno_genes_hypothetical_seqs[i,]$prot_sequence,names=lachno_genes_hypothetical_seqs[i,]$locus_tag,file.out = "hypothetical_proteins.fasta",open="a")
  }
  
  # run UCLUST
  # ./usearch -cluster_fast hypothetical_proteins.fasta -id 0.5 -centroids hypothetical_centroids.fasta -uc hypothetical_clusters.uc
  # convert .uc file to .csv
  
  # read in cluster file, link with isolates
  clusters <- read.csv("hypothetical_clusters.csv",header=TRUE, sep=",")
  clusters <- clusters %>%
    filter(record_type!='C') %>%
    select(locus_tag,cluster_num)
  
  
  lachno_genes_hypothetical_clusters <- inner_join(lachno_genes_hypothetical_seqs,
                                                   clusters,
                                                   by='locus_tag') %>%
    select(seq_id,cluster_num,count)
  
  #Generate list of genes in lachno with COGs/No COGS/Annotated/
  lachno_genes$count <- 1
  
  lachno_genes_annotated <- lachno_genes %>%
    select(seq_id,product,count) %>%
    filter(product!='hypothetical protein')
  
  lachno_genes_COG <- lachno_genes %>%
    select(seq_id,product,COG,count) %>%
    filter(product!='hypothetical protein') %>%
    filter(COG != '')
 
  Figure_S3_a_data_COGs <- aggregate(lachno_genes_COG$count,
                                     by=list(seq_id=lachno_genes_COG$seq_id),
                                     FUN=sum)
  
  colnames(Figure_S3_a_data_COGs) <- c("seq_id","COG")
  #lachno_genes_COG <- unique(lachno_genes_COG$product)
  
  
  lachno_genes_no_COG <- lachno_genes %>%
    select(seq_id,product,COG,count) %>%
    filter(product!='hypothetical protein') %>%
    filter(COG == '')
  #lachno_genes_no_COG$count <- 1 
  
  Figure_S3_a_data_products <- aggregate(lachno_genes_no_COG$count,
                                         by=list(seq_id=lachno_genes_no_COG$seq_id),
                                         FUN=sum)
  
  colnames(Figure_S3_a_data_products) <- c("seq_id","product")
  #lachno_genes_no_COG <- unique(lachno_genes_no_COG$product)
  
  lachno_genes_hypothetical <- lachno_genes %>%
    select(seq_id,product,COG,count) %>%
    filter(product=='hypothetical protein')
  #lachno_genes_hypothetical$count <- 1 
  
  Figure_S3_a_data_hypothetical <- aggregate(lachno_genes_hypothetical$count,
                                             by=list(seq_id=lachno_genes_hypothetical$seq_id),
                                             FUN=sum)
  
  colnames(Figure_S3_a_data_hypothetical) <- c("seq_id","hypothetical")
  
  Figure_S3_a_data <- Figure_S3_a_data_COGs %>%
    inner_join(Figure_S3_a_data_products,by="seq_id") %>%
    inner_join(Figure_S3_a_data_hypothetical,by="seq_id") %>%
    melt(id="seq_id")
  
  Figure_S3_a_data$variable <- factor(Figure_S3_a_data$variable,levels= c("hypothetical", "product","COG")) 
  
  reorder_figure_S3A <- read.csv(file="Figure_3A_reorder.csv", header=TRUE,sep=",")
  reorder_figure_S3A <- merge(reorder_figure_S3A,lachnospiraceae_working,by.x="toReorder",by.y="msk_id")
  Figure_S3_a_data <- merge(Figure_S3_a_data,reorder_figure_S3A,by="seq_id")
  
  figure_S3A_plot <- ggplot(Figure_S3_a_data, aes(x=properorder, y=value, fill = variable))+
    geom_bar(stat="identity",color="black",size=0.1)+theme_minimal()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.text.x = element_text(angle=90,size=10),axis.line.y = element_line(color="black",size=.2))+
    scale_fill_manual(values=c("#000000","#FBA27D","#6FC7CF" ))
  figure_S3A_plot
  
  
  Figure_S3_b_data <- Figure_S3_a_data %>%
    select(species,variable,value) %>%
    group_by(species,variable) %>%
    summarise_each(funs(mean, sd))
  Figure_S3_b_data$rsd <- 100*Figure_S3_b_data$sd/Figure_S3_b_data$mean
  
  Figure_S3_b_data <- cast(Figure_S3_b_data,species~variable)
  
  Figure_S3_b_plot <- ggplot(Figure_S3_b_data)+geom_density(aes(x=hypothetical),alpha=.6,fill="#000000")+geom_density(aes(x=product),alpha=.6,fill="#FBA27D")+geom_density(aes(x=COG),alpha=.6,fill="#6FC7CF")+
    theme_minimal()+
    theme(panel.grid.major = element_line(color="#000000",size=0.2), panel.grid.minor = element_blank(),strip.text.x = element_text(angle=90,size=10),axis.line.y = element_line(color="black",size=.2))
  
  Figure_S3_b_plot
  
  ## FIGURE S3C Plot
  
 
  figure_s3c_data_products <- lachno_genes %>%
    filter(product != 'hypothetical protein') %>%
    filter(COG=="")%>%
    select(length_bp)
    
  figure_s3c_data_products$type <- "product"
  figure_s3c_data_COGS <- lachno_genes %>%
    filter(product != 'hypothetical protein') %>%
    filter(COG!="")%>%
    select(length_bp) 
  figure_s3c_data_COGS$type <- "COG"
     
  Figure_s3c_data_hypothetical <- lachno_genes %>%
    filter(product =='hypothetical protein') %>%
    select(length_bp)
  
  Figure_s3c_data_hypothetical$type <- "hypothetical"
  
  Figure_s3c_master_lengths <- rbind(Figure_s3c_data_hypothetical,figure_s3c_data_COGS,figure_s3c_data_products)
  colnames(Figure_s3c_master_lengths) <- c("length","type")
  
  Figure_s3c_master_lengths$prot_length <- Figure_s3c_master_lengths$length/3
  Figure_s3c_master_lengths$logval <- log(Figure_s3c_master_lengths$prot_length,2)
  figure_s3_c <- ggplot(Figure_s3c_master_lengths,aes(x=logval,group=type,fill=type)) + geom_density(alpha=.6)+
    scale_fill_manual(values=c("#6FC7CF","#000000","#FBA27D"))+
    theme_minimal()+
    theme(panel.grid.major = element_line(color="#000000",size=0.2), panel.grid.minor = element_blank(),strip.text.x = element_text(angle=90,size=10),axis.line.y = element_line(color="black",size=.2),legend.position = "none")
  figure_s3_c
  
  
  anova_master_lengths <- aov(prot_length ~type,data=Figure_s3c_master_lengths)
  summary(anova_master_lengths)
  TukeyHSD(anova_master_lengths)
  
  #merge lists of annotated genes in each isolate with list of hypothetical protein clusters in each isolate
 colnames(lachno_genes_annotated) <- c("seq_id","cluster_num","count")
  lachno_genes_combo <- rbind(lachno_genes_annotated,lachno_genes_hypothetical_clusters)
  
  #Identifying Clusters by Blast
  
  tophits <- tophits %>%
    group_by(qseqid) %>%
    arrange(desc(bitscore)) %>%
    dplyr::slice(1) #
  
  top_hit_hypo <- tophits %>%
    filter(name %>% startsWith(" hypothetical protein"))
  top_hit_hypo_multi <- tophits %>% 
    filter(name %>% startsWith(" MULTISPECIES: hypothetical protein"))
  top_hit_hypo <- rbind(top_hit_hypo,top_hit_hypo_multi)
  tophits$include <- tophits$name %!in% top_hit_hypo$name
  tophits <- as.data.frame(tophits) %>%
    filter(include==TRUE)
  tophits$correct_names <- str_replace(tophits$name," MULTISPECIES:","")
  clusters_id <- read.csv("hypothetical_clusters.csv",header=TRUE, sep=",")
  colnames(tophits)[colnames(tophits)=="qseqid"] <- "locus_tag"
  identified_clusters <- clusters_id %>% 
    filter(record_type =='C') %>%
    merge(tophits,by='locus_tag',all.x=TRUE)
  
  identified_clusters$combo <- ifelse(is.na(identified_clusters$correct_names),identified_clusters$cluster_num,paste(identified_clusters$correct_name,identified_clusters$cluster_num))
  identified_clusters$type <- ifelse(is.na(identified_clusters$correct_names),"hypothetical protein","identified cluster")
  
  
  #colnames(lachno_genes_annotated)[colnames(lachno_genes_annotated)=="product"] <- "cluster_num"
  
  identified_cluster_list <- identified_clusters %>%
    filter(type=="identified cluster")
  identified_cluster_list <- unique(identified_cluster_list$cluster_num)
  lachno_id_cluster <- lachno_genes_combo %>%
    filter(cluster_num %in% identified_cluster_list)
  
  Figure_S3_d_data_identified<- aggregate(lachno_id_cluster$count,
                                             by=list(seq_id=lachno_id_cluster$seq_id),
                                             FUN=sum)
  colnames(Figure_S3_d_data_identified) <- c("seq_id","identified_cluster")

  Figure_S3_d_data <- Figure_S3_a_data_COGs %>%
    inner_join(Figure_S3_a_data_products,by="seq_id") %>%
    inner_join(Figure_S3_a_data_hypothetical,by="seq_id") %>%
    inner_join(Figure_S3_d_data_identified,by="seq_id")
    
  Figure_S3_d_data$hypothetical <- Figure_S3_d_data$hypothetical- Figure_S3_d_data$identified_cluster
  Figure_S3_d_data <- Figure_S3_d_data %>%
    melt(id="seq_id")
  
  Figure_S3_d_data$variable <- factor(Figure_S3_d_data$variable,levels= c("hypothetical","identified_cluster", "product","COG")) 
  
  Figure_S3_d_data <- merge(Figure_S3_d_data,reorder_figure_S3A,by="seq_id")
  
  
  figure_s3_d_plot <- ggplot(Figure_S3_d_data, aes(x=properorder, y=value, fill = variable))+
    geom_bar(stat="identity",color="black",size=0.1)+theme_minimal()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.text.x = element_text(angle=90,size=10),axis.line.y = element_line(color="black",size=.2))+
    scale_fill_manual(values=c("#000000","#89EF9A","#FBA27D","#6FC7CF" ))
  figure_s3_d_plot


# Figure S3 Average Blautia wexlerae Isolate Annotations  -----------------------------
  
  figure_s3_e_data <- Figure_S3_d_data %>%
    filter(species=="Blautia wexlerae") %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(mean(value)) %>%
    dplyr::ungroup(variable)
  colnames(figure_s3_e_data) <- c("type","average")
  figure_s3_e_data$percentage <- 100*(figure_s3_e_data$average/sum(figure_s3_e_data$average))
  
  Figure_S3_e <- ggplot(figure_s3_e_data,aes(x=1,y=percentage,fill=type))+
    geom_bar(stat="identity",color="black")+
    theme_minimal()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x = element_text(angle=90,size=10),
          axis.line.y = element_line(color="black",size=.2),
          legend.position = "none")+
    scale_fill_manual(values=c("#000000","#89EF9A","#FBA27D","#6FC7CF" ))
  Figure_S3_e
  
  
  
# Figure 3 Core Genome Analysis per Isolate -------------------------------

  
  #generate table of genes, with presence or absence in each isolate (isolates are rows, genes are columns))
  cast_lachno_genes <- cast(lachno_genes_combo,seq_id~cluster_num) %>%
    mutate_at(vars(-one_of("seq_id")),list(~ ifelse( . > 1, 1, .)))
  
  #determine core genome of entire Lachnospiraceae family
  
  
  totalcounts <- colSums(cast_lachno_genes[,-1])
  family_level_counts <- as.data.frame(totalcounts)
  
  family_level_counts$include <-family_level_counts$totalcounts == (nrow(cast_lachno_genes))
  family_level_counts$product <- row.names(family_level_counts)
  list_of_core_genes <- family_level_counts %>%
    filter(include==TRUE)
  list_of_core_genes$product <- as.character(list_of_core_genes$product)
  list_of_core_genes$hypothetical <- !is.na(as.numeric(list_of_core_genes$product))
  list_of_core_genes$cluster_num <- as.numeric(list_of_core_genes$product)
  
  family_hypothetical_count <- list_of_core_genes %>% 
    filter(hypothetical==TRUE) %>%
    nrow()
  
  lachno_genes_prokka_annotated <- lachno_genes_no_COG %>%
    filter(product %!in% lachno_genes_COG$product)
  family_with_cog_count <- list_of_core_genes %>%
    filter(product %in% lachno_genes_COG$product) 
  family_with_cog_count <- length(unique(family_with_cog_count$product))
  family_no_cog_count <- list_of_core_genes %>% 
    filter(product %in% lachno_genes_prokka_annotated$product) 
  family_no_cog_count <- length(unique(family_no_cog_count$product))
  

  family_identified_cluster_count <- merge(list_of_core_genes,identified_clusters,by='cluster_num',all.x=TRUE) #%>%
  family_identified_cluster_count <-  sum(family_identified_cluster_count$type=='identified cluster',na.rm = TRUE)
  family_hypothetical_count <- family_hypothetical_count - family_identified_cluster_count
  
  
  #Core genome of individual genera
  
  #list of genera
  lachnospiraceae_genera<- data.frame(table(lachnospiraceae_working$genus)) %>%
    filter(Freq >= 1)
  lachnospiraceae_genera_list <- paste(lachnospiraceae_genera$Var1)
  
  
  df_name_genus <- paste("lachno_",lachnospiraceae_genera$Var1,sep="")
  lachno_genus_ <- list()
  percent_genus_output <- list()
  matched_genus <- list()
  counts_genus_genes <- data.frame(genus=character(),genus_genes=integer(),number_in_genus=integer())
  filter_core_genes_by_master_genus <- data.frame("totalcountsgenus"=numeric(),"include"=logical(),"gene"=character(),"genus"=character(),number_of_species_in_genus =integer()) 

  for(i in 1:nrow(lachnospiraceae_genera)){

    print(i)

    
    lachno_genus_working <- lachnospiraceae_working %>% 
      filter(genus==lachnospiraceae_genera_list[i])
    print(paste("Filtering by",lachnospiraceae_genera_list[i]))
   
    
    cast_genus_genes <- cast_lachno_genes %>% 
      filter(seq_id%in% lachno_genus_working$seq_id)
    
    totalcountsgenus <- colSums(cast_genus_genes[,-1])
    
    countedgenus <- as.data.frame(totalcountsgenus)
    countedgenus$include <-countedgenus$totalcountsgenus == (nrow(cast_genus_genes))
    countedgenus$gene <- row.names(countedgenus)
    countedgenus$genus <-lachnospiraceae_genera_list[i]
    countedgenus$number_of_species_in_genus <- length(unique(lachno_genus_working$species))
    #shared is total for genus, includes those shared at family level as well
    
    core_genus_gene_names <- countedgenus %>% 
      filter(include==TRUE)
    
    filter_core_genes_by_master_genus <- rbind(filter_core_genes_by_master_genus,core_genus_gene_names)
    
   

    
  }
  
  filter_core_genes_by_master_genus <- filter_core_genes_by_master_genus %>%
    filter(number_of_species_in_genus > 1)
  
  
  filter_core_genes_by_master_genus$count <- 1
  Genus_core_plot <- aggregate(filter_core_genes_by_master_genus$count,by=list(genus=filter_core_genes_by_master_genus$genus),FUN=sum)
  colnames(Genus_core_plot) <- c("genus","genus_plot_not_corrected")
  
  filter_core_genes_by_master_genus <- filter_core_genes_by_master_genus %>%
    filter(gene %!in% list_of_core_genes$product)
  
  
  
  
  filter_core_genes_by_master_genus$hypothetical <- !is.na(as.numeric(filter_core_genes_by_master_genus$gene))
  genus_hypothetical <- filter_core_genes_by_master_genus %>%
    filter(hypothetical==TRUE)
  genus_hypothetical_count <- length(unique(genus_hypothetical$gene))
  
  genus_with_cog <- filter_core_genes_by_master_genus %>% 
    filter(gene %in% lachno_genes_COG$product)
  genus_with_cog_count <- length(unique(genus_with_cog$gene))
  
  genus_no_cog <- filter_core_genes_by_master_genus %>% 
    filter(gene %in% lachno_genes_prokka_annotated$product)
  genus_no_cog_count <- length(unique(genus_no_cog$gene))
  
  filter_core_genes_by_master_genus$cluster_num <- as.numeric(filter_core_genes_by_master_genus$gene)
  filter_core_genes_by_master_genus <- merge(filter_core_genes_by_master_genus,identified_clusters,by.x='cluster_num',by.y='cluster_num',all.x=TRUE)
  genus_identified_clusters <- filter_core_genes_by_master_genus %>% 
    filter(filter_core_genes_by_master_genus$type=='identified cluster')
  
  genus_identified_cluster_count <- length(unique(genus_identified_clusters$combo))
  genus_hypothetical_count <- genus_hypothetical_count - genus_identified_cluster_count

  #Core genome of individual species 
  
  #list of species
  lachnospiraceae_species<- data.frame(table(lachnospiraceae_working$species)) %>%
    filter(Freq >= 1)
  lachnospiraceae_species_list <- paste(lachnospiraceae_species$Var1)
  
  
  lachno_species_ <- list()
  percent_species_output <- list()
  matched_species <- list()
  counts_species_genes <- data.frame(species=character(),species_genes=integer(),number_in_species=integer())
  filter_core_genes_by_master_species <- data.frame("totalcountsgenus"=numeric(),"include"=logical(),"gene"=character(),"species"=character,"genus"=character(),number_of_isolates_in_species=integer()) 
  
  for(i in 1:nrow(lachnospiraceae_species)){
    
    print(i)
    
    
    lachno_species_working <- lachnospiraceae_working %>% 
      filter(species==lachnospiraceae_species_list[i])
    print(paste("Filtering by",lachnospiraceae_species_list[i]))
    
    
    cast_species_genes <- cast_lachno_genes %>% 
      filter(seq_id%in% lachno_species_working$seq_id)
    
    totalcountsspecies <- colSums(cast_species_genes[,-1])
    
    countedspecies <- as.data.frame(totalcountsspecies)
    countedspecies$include <-countedspecies$totalcountsspecies == (nrow(cast_species_genes))
    countedspecies$gene <- row.names(countedspecies)
    countedspecies$species <-lachnospiraceae_species_list[i]
    countedspecies$number_of_isolates_in_species <- length(unique(lachno_species_working$seq_id))
    countedspecies$genus <- unique(lachno_species_working$genus)
    
    
    core_species_gene_names <- countedspecies %>% 
      filter(include==TRUE)
    
    filter_core_genes_by_master_species <- rbind(filter_core_genes_by_master_species,core_species_gene_names)
    
    
    
  }
  filter_core_genes_by_master_species$count <- 1
 
  species_core_plot <- aggregate(filter_core_genes_by_master_species$count,by=list(species=filter_core_genes_by_master_species$species),FUN=sum)
  colnames(species_core_plot) <- c("species","species_plot_not_corrected_for_genus_family")
 
   filter_core_genes_by_master_species <- filter_core_genes_by_master_species %>%
    filter(gene %!in% list_of_core_genes$product)
  
  
  
  #filter_core_genes_by_master_species= list of genes that are core in one of the species group cores
  
  filter_core_genes_by_master_species <- filter_core_genes_by_master_species %>%
    filter(gene %!in% list_of_core_genes$product) %>%
    filter(gene %!in% filter_core_genes_by_master_genus$gene) %>%
    filter(number_of_isolates_in_species > 1)
  

  
  filter_core_genes_by_master_species$hypothetical <- !is.na(as.numeric(filter_core_genes_by_master_species$gene))
  
  species_hypothetical <- filter_core_genes_by_master_species %>%
    filter(hypothetical==TRUE)
  species_hypothetical_count <- length(unique(species_hypothetical$gene))
  
  species_with_cog <- filter_core_genes_by_master_species %>% 
    filter(gene %in% lachno_genes_COG$product)
  species_with_cog_count <- length(unique(species_with_cog$gene))
  
  species_no_cog <- filter_core_genes_by_master_species %>% 
    filter(gene %in% lachno_genes_prokka_annotated$product)
  species_no_cog_count <- length(unique(species_no_cog$gene))
  
  filter_core_genes_by_master_species$cluster_num <- as.numeric(filter_core_genes_by_master_species$gene)
  filter_core_genes_by_master_species <- merge(filter_core_genes_by_master_species,identified_clusters,by.x='cluster_num',by.y='cluster_num',all.x=TRUE)
  species_identified_clusters <- filter_core_genes_by_master_species %>% 
    filter(type=='identified cluster')
  
  species_identified_cluster_count <- length(unique(species_identified_clusters$combo))
  species_hypothetical_count <- species_hypothetical_count - species_identified_cluster_count

#DONOR level
  lachnospiraceae_working$donor <- lachnospiraceae_working$Donor
  
  #prep for loop
  
  tablename <- paste("lachno_",lachnospiraceae_species$Var1,sep="")
  lachno_ <- list()
  percent_species_output <- list()
  lachno_donor_ <- list()
  percent_donor_output <- list()
  master_donor_table <- data.frame(ID=character(),donor=character(),species=character(),identity=integer(),numberdonor=integer(),shared_donor_genes = integer())
  master_unique_donor <- data.frame(value=factor(),genes_unique_within_donor=integer(),seq_id=integer())
  filter_core_genes_by_master_donor <- data.frame("totalcounts_donor"=numeric(),"include"=logical(),"gene"=character(),"donor"=character(),"species"=character()) 
  
  for(i in 1:nrow(lachnospiraceae_species)){
    print(paste("i =",i))
    lachnospiraceae_species_list[i]
    lachno_species_working <- lachnospiraceae_working %>% 
      filter(species==lachnospiraceae_species_list[i])
    #now we have a table with just one species
    # now we count the donors in that table
    donorcount<- data.frame(table(lachno_species_working$donor)) %>%
      filter(Freq >=1)
    
    #prep for loop in the loop
    donor_filter <- paste(donorcount$Var1)
    
    for(x in 1:nrow(donorcount)){
      print(paste("x =",x))
      donor_filter[x]
      
      
      donor_table_name <- paste("lachno_",lachnospiraceae_species_list[i],donor_filter[x],sep="")
      #That gives us a unique name for the donor and species
      
      donortable <- lachno_species_working  %>% 
        filter(donor==donor_filter[x])
      #now we should have a table with 1 donor
      
      cast_donor_genes <- cast_lachno_genes %>% 
        filter(seq_id %in% donortable$seq_id)
      
      totalcounts_donor <- colSums(cast_donor_genes[,-1]!=0)
      counted_donor <- as.data.frame(totalcounts_donor)
      shared <- sum(counted_donor$totalcounts_donor == nrow(cast_donor_genes))
      
      counted_donor$include <-counted_donor$totalcounts_donor == (nrow(cast_donor_genes)) 
      counted_donor$gene <- row.names(counted_donor)
      counted_donor$species <-lachnospiraceae_species_list[i]
      counted_donor$donor <- donor_filter[x]
      counted_donor$number_of_isolates_in_donor <- length(unique(cast_donor_genes$seq_id))
      counted_donor$number_of_donors_for_species <- length(unique(lachno_species_working$Donor))
      core_donor_gene_names <- counted_donor %>% 
        filter(counted_donor$include==TRUE)
      filter_core_genes_by_master_donor <- rbind(filter_core_genes_by_master_donor,core_donor_gene_names)
      
    }
    
  }
  filter_core_genes_by_master_donor <- filter_core_genes_by_master_donor %>% 
    filter(number_of_isolates_in_donor >1) %>%
    filter(number_of_donors_for_species >1)
  filter_core_genes_by_master_donor$count <- 1
  donor_core_plot <- aggregate(filter_core_genes_by_master_donor$count,by=list(species=filter_core_genes_by_master_donor$species,donor=filter_core_genes_by_master_donor$donor),FUN=sum)
  colnames(donor_core_plot) <- c("species","donor","un_corrected_donor_plot")
  
  
  filter_core_genes_by_master_donor <- filter_core_genes_by_master_donor %>% 
    filter(gene %!in% list_of_core_genes$product) %>%
    filter(gene %!in% filter_core_genes_by_master_genus$gene) %>%
    filter(gene %!in% filter_core_genes_by_master_species$gene)
  

  
  filter_core_genes_by_master_donor$hypothetical <- !is.na(as.numeric(filter_core_genes_by_master_donor$gene))
  donor_total_unique_core <- length(unique(filter_core_genes_by_master_donor$gene))
  donor_hypothetical <- filter_core_genes_by_master_donor %>%
    filter(hypothetical==TRUE)
  donor_hypothetical_count <- length(unique(donor_hypothetical$gene))
  
  donor_with_cog <- filter_core_genes_by_master_donor %>% 
    filter(gene %in% lachno_genes_COG$product)
  donor_with_cog_count <- length(unique(donor_with_cog$gene))
  
  donor_no_cog <- filter_core_genes_by_master_donor %>% 
    filter(gene %in% lachno_genes_prokka_annotated$product)
  donor_no_cog_count <- length(unique(donor_no_cog$gene))
  
  filter_core_genes_by_master_donor$cluster_num <- as.numeric(filter_core_genes_by_master_donor$gene)
  filter_core_genes_by_master_donor <- merge(filter_core_genes_by_master_donor,identified_clusters,by.x='cluster_num',by.y='cluster_num',all.x=TRUE)
  donor_identified_clusters <- filter_core_genes_by_master_donor %>% 
    filter(type=='identified cluster')
  
  donor_identified_cluster_count <- length(unique(donor_identified_clusters$combo))
  donor_hypothetical_count <- donor_hypothetical_count - donor_identified_cluster_count

  
  filter_non_core <- lachno_genes_combo %>% 
    filter(cluster_num %!in% list_of_core_genes$product) %>%
    filter(cluster_num %!in% filter_core_genes_by_master_genus$gene) %>%
    filter(cluster_num %!in% filter_core_genes_by_master_species$gene) %>%
    filter(cluster_num%!in% filter_core_genes_by_master_donor$gene)
  
  colnames(filter_non_core) <- c("seq_id","gene","count")
  
non_core_per_isolate <- aggregate(filter_non_core$count,by=list(seq_id=filter_non_core$seq_id),FUN=sum)
colnames(non_core_per_isolate) <- c("seq_id","non_core_plot")
#colnames(filter_non_core) <- c("gene")

filter_non_core$hypothetical <- !is.na(as.numeric(filter_non_core$gene))
non_core_unique <- length(unique(filter_non_core$gene))
non_core_hypothetical <- filter_non_core %>%
  filter(hypothetical==TRUE)
non_core_hypothetical_count <- length(unique(non_core_hypothetical$gene))

non_core_with_cog <- filter_non_core %>% 
  filter(gene %in% lachno_genes_COG$product)
non_core_with_cog_count <- length(unique(non_core_with_cog$gene))

non_core_no_cog <- filter_non_core %>% 
  filter(gene %in% lachno_genes_prokka_annotated$product)
non_core_no_cog_count <- length(unique(non_core_no_cog$gene))

filter_non_core$cluster_num <- as.numeric(filter_non_core$gene)
filter_non_core <- merge(filter_non_core,identified_clusters,by.x='cluster_num',by.y='cluster_num',all.x=TRUE)
non_core_identified_clusters <- filter_non_core %>% 
  filter(type=='identified cluster')

non_core_identified_cluster_count <- length(unique(non_core_identified_clusters$combo))
non_core_hypothetical_count <- non_core_hypothetical_count - non_core_identified_cluster_count
#Figure 3b: plot:
  #Family
  family_hypothetical_count
  family_identified_cluster_count
  family_no_cog_count
  family_with_cog_count
  #genus level
  genus_hypothetical_count
  genus_identified_cluster_count
  genus_with_cog_count
  genus_no_cog_count
  #Species level
  species_hypothetical_count
  species_identified_cluster_count
  species_no_cog_count
  species_with_cog_count
  #donor level
  donor_hypothetical_count
  donor_identified_cluster_count
  donor_no_cog_count
  donor_with_cog_count
  #non_core
  non_core_hypothetical_count
  non_core_identified_cluster_count
  non_core_no_cog_count
  non_core_with_cog_count
  

total_genes_per_isolate <- lachno_genes_combo %>%
  group_by(seq_id) %>%
  summarise(n_distinct(cluster_num))

colnames(total_genes_per_isolate) <- c("seq_id","total_genes")
lachnospiraceae_working$donor <- as.character(lachnospiraceae_working$donor)
Figure_3A_data <- lachnospiraceae_working %>%
  select(seq_id,genus,species,donor) %>%
  inner_join(total_genes_per_isolate,by="seq_id")

Figure_3A_data$family_plot <- nrow(list_of_core_genes)

Figure_3A_data <- merge(Figure_3A_data,Genus_core_plot,by="genus",all.x=TRUE)
Figure_3A_data$genus_plot_corrected <- Figure_3A_data$genus_plot_not_corrected-Figure_3A_data$family_plot
Figure_3A_data <- merge(Figure_3A_data,species_core_plot,by="species",all.x=TRUE)

Figure_3A_data[is.na(Figure_3A_data)] <- 0

Figure_3A_data$species_plot_corrected <- Figure_3A_data$species_plot_not_corrected_for_genus_family - Figure_3A_data$genus_plot_corrected - Figure_3A_data$family_plot
Figure_3A_data <- merge(Figure_3A_data,donor_core_plot,by=c("species","donor"),all.x=TRUE)

Figure_3A_data$donor_plot_corrected <- Figure_3A_data$un_corrected_donor_plot-Figure_3A_data$species_plot_corrected - Figure_3A_data$genus_plot_corrected - Figure_3A_data$family_plot

Figure_3A_data[is.na(Figure_3A_data)] <- 0

Figure_3A_data$non_core_plot <- Figure_3A_data$total_genes - Figure_3A_data$donor_plot_corrected -Figure_3A_data$species_plot_corrected - Figure_3A_data$genus_plot_corrected - Figure_3A_data$family_plot

Figure_3A_data <- Figure_3A_data %>%
  select(seq_id,genus,species,donor,total_genes,family_plot,genus_plot_corrected,species_plot_corrected,donor_plot_corrected,non_core_plot)

Figure_3A_data_melt <- melt(Figure_3A_data, id=c("seq_id","genus","species","donor","total_genes"))
Figure_3A_data_melt$variable <- factor(Figure_3A_data_melt$variable,levels= c("non_core_plot", "donor_plot_corrected","species_plot_corrected","genus_plot_corrected","family_plot"))

reorder_figure_3A <- read.csv(file="Figure_3A_reorder.csv", header=TRUE,sep=",")
reorder_figure_3A <- merge(reorder_figure_3A,lachnospiraceae_working,by.x="toReorder",by.y="msk_id")
Figure_3A_data_melt <- merge(Figure_3A_data_melt, reorder_figure_3A, by.x = "seq_id", by.y = "seq_id", all.x = TRUE,all.y=TRUE)

Figure_3A_plot <- ggplot(Figure_3A_data_melt, aes(x=properorder, y=value, fill = variable))+
  geom_bar(stat="identity",color="black",size=0.1)+theme_minimal()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.text.x = element_text(angle=90,size=10),axis.line.y = element_line(color="black",size=.2))+
  scale_fill_manual(values=c("#778899","#FFFFFF","#ADADAD","#636363","#000000"  ))
  theme(legend.position = "none")
Figure_3A_plot

Figure_3A_stats <- Figure_3A_data

Figure_3A_stats$percent_species_shared = 100*(Figure_3A_stats$family_plot+Figure_3A_stats$genus_plot_corrected+Figure_3A_stats$species_plot_corrected)/Figure_3A_stats$total_genes
mean(Figure_3A_stats$percent_species_shared)
# Figure S3F Average B. wexlerae Core Genome------------------------------------------

figure_s3_f_data <- Figure_3A_data_melt %>%
  filter(species.x=="Blautia wexlerae") %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(mean(value))

colnames(figure_s3_f_data) <- c("type","average")

figure_s3_f_data$percentage <- 100*(figure_s3_f_data$average/sum(figure_s3_f_data$average))

Figure_S3_f <- ggplot(figure_s3_f_data,aes(x=1,y=percentage,fill=type))+
  geom_bar(stat="identity",color="black")+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(angle=90,size=10),
        axis.line.y = element_line(color="black",size=.2),
        legend.position = "none")+
  scale_fill_manual(values=c("#778899","#FFFFFF","#ADADAD","#636363","#000000" ))
Figure_S3_f

# Figure 3 Metabolic Pathways ---------------------------------------------

#generate lists of kegg pathways associated with genes that become part of the core genome at different taxonomic ranks
  #family
  family_kegg_lookup <- lachno_genes %>% 
    filter(product %in% list_of_core_genes$product) %>%
    select(seq_id,product,locus_tag) %>%
    inner_join(prokka_locus_uniprot_lookup,by='locus_tag') %>%
    select(seq_id,locus_tag,uniprot)
  family_kegg_lookup$uniprot <- paste("up.",family_kegg_lookup$uniprot,sep="")  

  family_kegg_pathways <- family_kegg_lookup %>%
    inner_join(kegg_uniprot_lookup,by="uniprot")
  family_kegg_pathways <- unique(family_kegg_pathways$kegg_id)
  family_kegg_table <- kegg_pathway_ko %>%
    filter(kegg_id %in% family_kegg_pathways)
 
  kegg_summary <- family_kegg_table %>%
    dplyr::count(level3, sort=TRUE)
  kegg_summary$family_rank <- seq.int(nrow(kegg_summary))
  colnames(kegg_summary) <- c("level3","family_count","family_rank")
   
  #genus"
  genus_kegg_lookup <- lachno_genes %>% 
    filter(product %in% filter_core_genes_by_master_genus$gene) %>%
    select(seq_id,product,locus_tag) %>%
    inner_join(prokka_locus_uniprot_lookup,by='locus_tag') %>%
    select(seq_id,locus_tag,uniprot)
  genus_kegg_lookup$uniprot <- paste("up.",genus_kegg_lookup$uniprot,sep="")  

  genus_kegg_pathways <- genus_kegg_lookup %>%
    inner_join(kegg_uniprot_lookup,by="uniprot")
  genus_kegg_pathways <- unique(genus_kegg_pathways$kegg_id)
  genus_kegg_table <- kegg_pathway_ko %>%
    filter(kegg_id %in% genus_kegg_pathways)
  genus_kegg_summary <- genus_kegg_table %>%
    dplyr::count(level3,sort=TRUE)
  genus_kegg_summary$genus_rank <- seq.int(nrow(genus_kegg_summary))
  colnames(genus_kegg_summary) <- c("level3","genus_count","genus_rank")
  
  master_kegg_summary <- merge(kegg_summary,genus_kegg_summary,by="level3",all.x=TRUE,all.y=TRUE)
  #species
  species_kegg_lookup <- lachno_genes %>% 
    filter(product %in% filter_core_genes_by_master_species$gene) %>%
    select(seq_id,product,locus_tag) %>%
    inner_join(prokka_locus_uniprot_lookup,by='locus_tag') %>%
    select(seq_id,locus_tag,uniprot)
  species_kegg_lookup$uniprot <- paste("up.",species_kegg_lookup$uniprot,sep="")  

  species_kegg_pathways <- species_kegg_lookup %>%
    inner_join(kegg_uniprot_lookup,by="uniprot")
  species_kegg_pathways <- unique(species_kegg_pathways$kegg_id)
  species_kegg_table <- kegg_pathway_ko %>%
    filter(kegg_id %in% species_kegg_pathways)
  
  species_kegg_summary <- species_kegg_table %>%
    dplyr::count(level3,sort=TRUE)
  species_kegg_summary$species_rank <- seq.int(nrow(species_kegg_summary))
  colnames(species_kegg_summary) <- c("level3","species_count","species_rank")
  master_kegg_summary <- merge(master_kegg_summary,species_kegg_summary,by="level3",all.x=TRUE,all.y=TRUE)
  
  #donor
  donor_kegg_lookup <- lachno_genes %>% 
    filter(product %in% filter_core_genes_by_master_donor$gene) %>%
    select(seq_id,product,locus_tag) %>%
    inner_join(prokka_locus_uniprot_lookup,by='locus_tag') %>%
    select(seq_id,locus_tag,uniprot)
  donor_kegg_lookup$uniprot <- paste("up.",donor_kegg_lookup$uniprot,sep="")  

  donor_kegg_pathways <- donor_kegg_lookup %>%
    inner_join(kegg_uniprot_lookup,by="uniprot")
  donor_kegg_pathways <- unique(donor_kegg_pathways$kegg_id)
  donor_kegg_table <- kegg_pathway_ko %>%
    filter(kegg_id %in% donor_kegg_pathways)
  donor_kegg_table <- donor_kegg_table %>%
    dplyr::count(level3,sort=TRUE)
  
  donor_kegg_pathways <- donor_kegg_lookup %>%
    inner_join(kegg_uniprot_lookup,by="uniprot")
  donor_kegg_pathways <- unique(donor_kegg_pathways$kegg_id)
  donor_kegg_table <- kegg_pathway_ko %>%
    filter(kegg_id %in% donor_kegg_pathways)
  
  donor_kegg_summary <- donor_kegg_table %>%
    dplyr::count(level3,sort=TRUE)
  donor_kegg_summary$donor_rank <- seq.int(nrow(donor_kegg_summary))
  colnames(donor_kegg_summary) <- c("level3","donor_count","donor_rank")
  
  master_kegg_summary <- merge(master_kegg_summary,donor_kegg_summary,by="level3",all.x=TRUE,all.y=TRUE)
  
  #noncore
  non_core_kegg_lookup <- lachno_genes %>% 
    filter(product %in% filter_non_core$gene) %>%
    select(seq_id,product,locus_tag) %>%
    inner_join(prokka_locus_uniprot_lookup,by='locus_tag') %>%
    select(seq_id,locus_tag,uniprot)
  non_core_kegg_lookup$uniprot <- paste("up.",non_core_kegg_lookup$uniprot,sep="")  
  
  non_core_kegg_pathways <- non_core_kegg_lookup %>%
    inner_join(kegg_uniprot_lookup,by="uniprot")
  non_core_kegg_pathways <- unique(non_core_kegg_pathways$kegg_id)
  non_core_kegg_table <- kegg_pathway_ko %>%
    filter(kegg_id %in% non_core_kegg_pathways)
  non_core_kegg_summary <- non_core_kegg_table %>%
    dplyr::count(level3,sort=TRUE)
  non_core_kegg_summary$donor_rank <- seq.int(nrow(non_core_kegg_summary))
  colnames(non_core_kegg_summary) <- c("level3","non_core_count","non_core_rank")
  
  master_kegg_summary <- merge(master_kegg_summary,non_core_kegg_summary,by="level3",all.x=TRUE,all.y=TRUE)
  
  master_kegg_summary$total <- master_kegg_summary$family_count+master_kegg_summary$species_count+master_kegg_summary$genus_count
  
  #list of pathways with 20 or more genes
  selected_pathways_10 <- master_kegg_summary %>%
    filter(total >19)
  
  master_kegg_summary$family_percent <- master_kegg_summary$family_count/(sum(master_kegg_summary$family_count,na.rm=TRUE))
  master_kegg_summary$genus_percent <- master_kegg_summary$genus_count/(sum(master_kegg_summary$genus_count,na.rm=TRUE))
  master_kegg_summary$species_percent <- master_kegg_summary$species_count/(sum(master_kegg_summary$species_count,na.rm=TRUE))
  master_kegg_summary$donor_percent <- master_kegg_summary$donor_count/(sum(master_kegg_summary$donor_count,na.rm=TRUE))
  master_kegg_summary$non_core_percent <- master_kegg_summary$non_core_count/(sum(master_kegg_summary$non_core_count,na.rm=TRUE))
  master_kegg_summary[is.na(master_kegg_summary)] <- 0
  
  level1_summary <- kegg_pathway_ko %>%
    group_by(level3) %>%
    dplyr::count(level1) %>%
    arrange(desc(n)) %>%
    top_n(1)
  
  figure_3_c_data <- master_kegg_summary %>%
    select(level3,family_percent,genus_percent,species_percent,donor_percent,non_core_percent) %>%
    inner_join(level1_summary,by="level3") %>%
    select(level3,family_percent,genus_percent,species_percent,donor_percent,non_core_percent,level1)
  
  figure_3_c_data <- melt(figure_3_c_data,id=c("level1","level3"))
  figure_3_c_data <- figure_3_c_data %>%
    dplyr::group_by(variable,level1) %>%
    dplyr::summarise(sum(value))
  figure_3_c_data
  kegg_heatmap <- master_kegg_summary %>%
    select(level3,family_percent,genus_percent,species_percent) %>%
    filter(level3 %in% selected_pathways_10$level3) %>%
    inner_join(level1_summary,by="level3")
  
 
 
    
  
  
  kegg_heatmap$level3 <- substr(kegg_heatmap$level3,9,nchar(kegg_heatmap$level3))
  row.names(kegg_heatmap) <- kegg_heatmap[,1]
  kegg_heatmap_annotate <- as.data.frame(kegg_heatmap[,5])
  colnames(kegg_heatmap_annotate) <- c("Pathway type")
  vector_heat <- as.character(kegg_heatmap_annotate$`Pathway type`)
  level1_annotation = HeatmapAnnotation(Path_type=vector_heat,which="row",width = unit(.10, "cm"),col=list(Path_type=c(" 09100 Metabolism"="#4A385B"," 09120 Genetic Information Processing"="#51b3e5"," 09130 Environmental Information Processing"="#EA60C0"," 09140 Cellular Processes"="#4FAF5D"," 09160 Human Diseases"="#000000")))
  kegg_heatmap <- kegg_heatmap[,c(2,3,4)]
  row.names(kegg_heatmap) <- str_replace(row.names(kegg_heatmap),"metabolism","meta.")
  row.names(kegg_heatmap) <- str_replace(row.names(kegg_heatmap),"biosynthesis","biosyn.")
  row.names(kegg_heatmap) <- str_replace(row.names(kegg_heatmap),"interconversions","intercon.")
  row.names(kegg_heatmap) <- str_replace(row.names(kegg_heatmap),"system","sys.")
  
  kegg_heatmap_matrix <- as.matrix(kegg_heatmap)
  Figure_3D <- ComplexHeatmap::Heatmap(kegg_heatmap_matrix,
                                  cluster_columns = FALSE,
                                  col=colorRamp2(c(0,0.035,.15),
                                                 c("#FFFFFF","#383838","#000000")),
                                  column_names_gp=gpar(fontsize=5),
                                  row_names_gp = gpar(fontsize=5),
                                  rect_gp = gpar(col = "black", lwd = 0.5),
                                  row_dend_gp = gpar(col="black",lwd=0.5),
                                  width = unit(1.5,"cm"),
                                  height=unit((0.15*nrow(kegg_heatmap_matrix)),"cm"),
                                  left_annotation = level1_annotation)
  
  Figure_3D
  
# Figure 4 Core_genome_trees on annotated products-------------------------------------------------------

#get list of core genes across lachnospiraceae that are annotated by prokka
annotated_core_genes <- list_of_core_genes %>%
    filter(hypothetical==FALSE)
  list_of_core_genes_align <- annotated_core_genes$product

#Save list
#write.csv(file="3_4_20_core_gene_list.csv",list_of_core_genes_align)

## Amino acid alignments:
core_gene_sequences <- list()
core_trees <- list()
core_alignments <-list()
core_alignments_seqinr <- list()
core_distances <- list()
core_njs <- list()

for(i in 1:nrow(annotated_core_genes)){
  print(i) 
  #take the longest copy in the case of multiple copies of a gene
  core_gene_loci_list <- prokka_genes %>%
     filter (seq_id %in% lachnospiraceae_working$seq_id,prokka_genes$product == list_of_core_genes_align[i]) %>%
     group_by(seq_id) %>%
     arrange(-length_bp) %>%
     dplyr::slice(1)

  core_gene_sequences <- prokka_seq %>%
    filter(locus_tag %in% core_gene_loci_list$locus_tag) %>%
    inner_join(core_gene_loci_list,by="locus_tag") %>%
    inner_join(lachnospiraceae_working,by="seq_id")
  
  #get protein sequence
  AA_xstring <- AAStringSet(core_gene_sequences$prot_sequence)
  names(AA_xstring) <- core_gene_sequences$msk_id
  
  #Align with muscle
  print("starting alignment")
  myAln <- msa(AA_xstring,"Muscle")
  
  #store alignment in two formats
  core_alignments[[i]] <- myAln
  print("finishing alignment")
  
  myAln2 <- msaConvert(myAln, type="seqinr::alignment")
  core_alignments_seqinr[[i]] <- myAln2
  
  #generate an NJ tree on each gene --> used in ASTRAL-III approach
  d <- dist.alignment(myAln2, "identity")
  core_distances[[i]] <- d
  myTree <- njs(d)
  core_njs[[i]] <- myTree
  
  myTree2 <- ggtree(myTree,layout="circular",branch.length = "none")+geom_treescale()
  
  tip_labels <- core_gene_sequences  %>% 
    select(msk_id,species,gene)
  tip_labels$seq <- tip_labels$msk_id
  tip_labels <- tip_labels %>%
    select(seq,species,gene)
  tip_labels$species <- as.character(tip_labels$species)
  tip_labels$gene <- as.character(tip_labels$gene)
  core_trees[[i]] <- myTree2 %<+% tip_labels + 
    geom_tippoint(aes(fill=species),size=2,pch=21,color="black",alpha=0.7)+
    scale_fill_manual(values=species.colors)+
    theme(legend.position = "none")+
    ggtitle(label=list_of_core_genes_align[i])
  
}

#Output all the Fasta files of the alignments
for(i in 1:nrow(annotated_core_genes)){
  file_name <- paste("gene",i,".fasta",sep="")
  gene_1_aa <- msaConvert(core_alignments[[i]],type="ape::AAbin")
  write.dna(gene_1_aa,file=file_name,format="fasta",nbcol=1,colw=60)
}

# read in fasta files and make data frame - protein

msk_id_list <- unique(lachnospiraceae_working$msk_id)
master_core_genome_alignments <- as.data.frame(msk_id_list)
colnames(master_core_genome_alignments) <- c("msk_id")

for(i in 1:nrow(annotated_core_genes)){
  print(i)
  file_name <- paste("gene",i,".fasta",sep="")
  individual_gene_alignment <- read.aa(file=file_name,bin=FALSE)
  data <- t(as.data.frame(individual_gene_alignment))
  data <- as.data.frame(data)
  data$V1 <- as.character(data$V1)
  colnames(data) <- list_of_core_genes[i]
  data$msk_id <- rownames(data)
  data <- data %>%
    select(msk_id,1)
  master_core_genome_alignments <- master_core_genome_alignments %>%
    inner_join(data,by="msk_id")
}

seq <- master_core_genome_alignments %>%
  reshape2::melt(id.vars="msk_id") %>%
  arrange(variable) %>%
  group_by(msk_id) %>%
  summarize(big_seq=paste(value,collapse="")) %>%
  mutate(big_seq=str_replace_all(big_seq,pattern="X",replacement = "-"))

for (i in 1:(nrow(seq))){
  write.fasta(seq[i,]$big_seq,seq[i,]$msk_id,file.out="big_seq3.fasta",open="a")
}

#Align concatenated sequences using MUSCLE.  Generate neighbor joining tree

concat_seq_tree <- read.tree(file="concatenate_aa_muscle.phy") %>%
  root("MSK.17.84")

concat_seq_tree_plot <- ggtree(concat_seq_tree,layout="circular",size=0.2)+
  ggtree::geom_hilight(node=352,alpha=0.4,fill="#CEF948",extendto=.35)+ #Rectale
  ggtree::geom_hilight(node=321,alpha=0.4,fill="#8AD188",extendto=.35)+#eutactus
  ggtree::geom_hilight(node=325,alpha=0.4,fill="#8C855E",extendto=.35)+#hadrus
  ggtree::geom_hilight(node=476,alpha=0.4,fill="#2E9DB7",extendto=.35)+#symbo
  ggtree::geom_hilight(node=475,alpha=0.4,fill="#AFF8DB",extendto=.35)+#aldense
  ggtree::geom_hilight(node=467,alpha=0.4,fill="#254A9E",extendto=.35)+#clostridoforme
  ggtree::geom_hilight(node=411,alpha=0.4,fill="#DAAD86",extendto=.35)+#selimon
  ggtree::geom_hilight(node=377,alpha=0.4,fill="#F4B571",extendto=.35)+#gnavus
  ggtree::geom_hilight(node=464,alpha=0.4,fill="#FFF5BA",extendto=.35)+#nexilis
  ggtree::geom_hilight(node=443,alpha=0.4,fill="#638E5E",extendto=.35)+#comes
  ggtree::geom_hilight(node=439,alpha=0.4,fill="#ACE7FF",extendto=.35)+#formici
  ggtree::geom_hilight(node=422,alpha=0.4,fill="#383838",extendto=.35)+#unclass
  ggtree::geom_hilight(node=420,alpha=0.4,fill="#877B5C",extendto=.35)+#scin
  ggtree::geom_hilight(node=423,alpha=0.4,fill="#6EB5FF",extendto=.35)+#longica
  ggtree::geom_hilight(node=478,alpha=0.4,fill="#A44DF7",extendto=.35)+#faecalicatena
  ggtree::geom_hilight(node=200,alpha=0.4,fill="#8221F4",extendto=.35)+#hansenii
  ggtree::geom_hilight(node=484,alpha=0.4,fill="#B5B5B5",extendto=.35)+#producta
  ggtree::geom_hilight(node=219,alpha=0.4,fill="#242851",extendto=.35)+#celercres
  ggtree::geom_hilight(node=487,alpha=0.4,fill="#DCD3FF",extendto=.35)+#f. saccharo
  ggtree::geom_hilight(node=522,alpha=0.4,fill="#FFB5E8",extendto=.35)+ #schinkii
  ggtree::geom_hilight(node=244,alpha=0.4,fill="#ED3580",extendto=.35)+ #luti
  ggtree::geom_hilight(node=506,alpha=0.4,fill="#ED3580",extendto=.35)+ #luti
  ggtree::geom_hilight(node=253,alpha=0.4,fill="#4C4C4C",extendto=.35)+ #caecimuris
  ggtree::geom_hilight(node=529,alpha=0.4,fill="#FFCCF9",extendto=.35)+ #obeum
  ggtree::geom_hilight(node=261,alpha=0.4,fill="#A79AFF",extendto=.35)+ #glucerase
  ggtree::geom_hilight(node=260,alpha=0.4,fill="#A79AFF",extendto=.35)+ #glucerase
  ggtree::geom_hilight(node=256,alpha=0.4,fill="#A79AFF",extendto=.35)+ #glucerase
  ggtree::geom_hilight(node=257,alpha=0.4,fill="#9E2556",extendto=.35)+ #faecis
  ggtree::geom_hilight(node=541,alpha=0.4,fill="#9E2556",extendto=.35)+ #faecis
  ggtree::geom_hilight(node=258,alpha=0.4,fill="#9E2556",extendto=.35)+ #faecis
  ggtree::geom_hilight(node=259,alpha=0.4,fill="#9E2556",extendto=.35)+ #faecis
  ggtree::geom_hilight(node=311,alpha=0.4,fill="#F6A6FF",extendto=.35)+ #Wexlerae
  ggtree::geom_hilight(node=85,alpha=0.4,fill="#F2B15C",extendto=.35)#rose
  #geom_text(aes(label=node), size=1,color="red")
concat_seq_tree_plot

concat_seq_tree_labels <- lachnospiraceae_working %>%
  select(msk_id,species)

Figure_4_concat <- concat_seq_tree_plot %<+% concat_seq_tree_labels + 
  geom_tippoint(aes(fill=species),size=1.5,pch=21,color="black",alpha=0.7)+
  scale_fill_manual(values=species.colors)+
  theme(legend.position = "none")
Figure_4_concat  


# Figure 5A & B UMAP Analysis --------------------------------------------------
#RUN UMAP for Figure 5A:
set.seed(1)
stop<- ncol(cast_lachno_genes)
custom_config <- umap.defaults
custom_config$n_neighbors = 100
custom_config$metric ='manhattan'
custom_config
comparison_map <- umap(cast_lachno_genes[,2:stop],custom_config)


coord <- as.data.frame(comparison_map$layout)
coord$seq_id <- cast_lachno_genes$seq_id

coord <- merge(coord,lachnospiraceae_working,by='seq_id')

Figure_5A <- ggplot(coord,aes(x=V1,y=V2))+
  geom_point(aes(fill=species),pch=21,color="black",size=4,alpha=0.7)+
  scale_fill_manual(values=species.colors) +
  theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"))+
  theme(legend.background=element_blank())

Figure_5A
glimpse(lachno_genes_annotated)



# Overlay pH values

coord_pH <- merge(coord,lachno_pH,by="msk_id",all.x=TRUE)
coord_pH_pos <- merge(coord,lachno_pH,by="msk_id")

Figure_5B <- ggplot(coord_pH,aes(x=V1,y=V2))+
  geom_point(aes(fill=mean_pH),pch=21,color="black",size=4)+
  scale_fill_gradientn(colors=c("#160993","#AFE2FA","#ED3780"),na.value = "#EAEAEA")+
  theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"))+
  theme(legend.background=element_blank())+
  theme(legend.position = "none",axis.title = element_blank())

Figure_5B

grid.arrange(Figure_5B,Figure_5B,ncol=2)


# Figure 5C and Figure S4B ------------------------------------------------



butyrate_pathway <- lachnospiraceae_working %>%
  select(seq_id) %>%
  left_join(lookup) %>%
  left_join(prokka_genes %>%
              #            filter(grepl(paste(butyrate_genes$prokka_target,collapse="|"),product)) %>%
              mutate(product=gsub("Crotonyl-CoA hydratase","Short-chain-enoyl-CoA hydratase",product)) %>%
              mutate(product=gsub("Phosphate acetyltransferase","Ethanolamine utilization protein EutD",product)) %>%
              mutate(product=gsub("Butyrate--acetoacetate CoA-transferase subunit A","Acetate CoA-transferase subunit alpha",product)) %>%
              mutate(product=gsub("Butyrate--acetoacetate CoA-transferase subunit B","Acetate CoA-transferase subunit beta",product)) %>%
              filter(product %in% butyrate_genes$prokka_target) %>%
              mutate(present=1)) %>%
  group_by(seq_id,product) %>%
  summarize(num_genes=sum(present)) %>%
  left_join(butyrate_genes %>%
              select(product=prokka_target,class)) %>%
  replace_na(list(present=0)) %>%
  reshape2::dcast(class+product ~ seq_id,fill=0,value.var="num_genes") %>%
  reshape2::melt(id.vars=c("product","class")) %>%
  dplyr::rename(seq_id=variable,num_genes=value) %>%
  mutate(num_genes=ifelse(product=="Phosphate acetyltransferase",num_genes-1,num_genes),
         present=ifelse(num_genes>1,1,num_genes)) %>%
  filter(!is.na(class)) %>%
  group_by(class,seq_id) %>%
  summarize(total=sum(present)) %>%
  ungroup() %>%
  left_join(butyrate_genes %>%
              group_by(class) %>%
              summarize(ptotal=n()) %>%
              ungroup()) %>%
  mutate(butyrate=ifelse(total==ptotal,"complete",
                         ifelse(total==0,"missing",
                                "incomplete"))) %>%
  left_join(lookup) %>%
  left_join(lachnospiraceae_working %>%
              select(msk_id,seq_id,species)) %>%
  filter(seq_id %in% lachnospiraceae_working$seq_id)


pathway_butyrate_acetate <- butyrate_pathway %>%
  #filter(butyrate=='complete') %>%
  filter(class=='B')
pathway_butyrate_acetoacetate <- butyrate_pathway %>%
  #filter(butyrate=='complete') %>%
  filter(class=='A')
pathway_butyrate_butyrate_kinase <- butyrate_pathway %>%
  #filter(butyrate=='complete') %>%
  filter(class=='F')



butyrate_isolates <- butyrate_pathway %>%
  select(seq_id,butyrate) %>%
  filter(butyrate=="complete") %>%
  group_by(seq_id) %>%
  summarise()

butyrate_isolates$include <- 1
butyrate_isolates <- butyrate_isolates %>%
  select(seq_id,include)

coord_butyrate_genes <- coord %>%
  select(seq_id,V1,V2) %>%
  merge(butyrate_isolates,by="seq_id",all.x=TRUE)

Figure_5C  <- ggplot(coord_butyrate_genes,aes(x=V1,y=V2))+
  geom_point(aes(fill=as.factor(include)),pch=21,color="black",size=4,alpha=0.7)+
  scale_fill_manual(values=c("red"))+
  theme(legend.position="none",panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"))+
  theme(legend.background=element_blank(),axis.title = element_blank())
Figure_5C

coord_butyrate_genes <- coord %>% 
  select(seq_id,V1,V2) %>%
  inner_join(pathway_butyrate_acetate,by="seq_id")

Figure_S4_acetate <- ggplot(coord_butyrate_genes,aes(x=V1,y=V2))+
  geom_point(aes(color=total),size=4,alpha=0.7)+
  scale_color_gradient(low="#FFFFFF",high="#383838",na.value = "#EAEAEA")+
  geom_point(aes(fill=butyrate,alpha=butyrate),pch=21,color="black",size=4)+
  scale_fill_manual(values=c("blue",NA,NA)) +
  scale_alpha_manual(values=c(.7,0.7,0.7))+
  theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"))+
  theme(legend.background=element_blank())
Figure_S4_acetate

coord_butyrate_genes <- coord %>% 
  select(seq_id,V1,V2) %>%
  inner_join(pathway_butyrate_butyrate_kinase,by="seq_id")

Figure_S4_butyrate_kinase <- ggplot(coord_butyrate_genes,aes(x=V1,y=V2))+
  geom_point(aes(color=total),size=4,alpha=0.7)+
  scale_color_gradient(low="#FFFFFF",high="#383838",na.value = "#EAEAEA")+
  geom_point(aes(fill=butyrate),pch=21,color="black",size=4,alpha=0.7)+
  scale_fill_manual(values=c("green",NA,NA)) +
  theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"))+
  theme(legend.background=element_blank())
Figure_S4_butyrate_kinase

coord_butyrate_genes <- coord %>% 
  select(seq_id,V1,V2) %>%
  inner_join(pathway_butyrate_acetoacetate,by="seq_id")
Figure_S4_acetoacetate <- ggplot(coord_butyrate_genes,aes(x=V1,y=V2))+
  geom_point(aes(color=total),size=4,alpha=0.7)+
  scale_color_gradient(low="#FFFFFF",high="#383838",na.value = "#EAEAEA")+
  geom_point(aes(fill=butyrate),pch=21,color="black",size=4,alpha=0.7)+
  scale_fill_manual(values=c("darkred",NA,NA)) +
  theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"))+
  theme(legend.background=element_blank())
Figure_S4_acetoacetate

# Figure S4A --------------------------------------------------------------

#generate table of genes, with presence or absence in each isolate for either hypothetical protein clusters or annotated proteins (isolates are rows, genes are columns))
cast_lachno_hypothetical_genes <- cast(lachno_genes_hypothetical_clusters,seq_id~cluster_num) %>%
  mutate_at(vars(-one_of("seq_id")),list(~ ifelse( . > 1, 1, .)))
cast_lachno_annotated_genes <- cast(lachno_genes_annotated,seq_id~cluster_num) %>%
  mutate_at(vars(-one_of("seq_id")),list(~ ifelse( . > 1, 1, .)))

#FIGURE S5 A - Hypothetical protein clusters alone
set.seed(1)
stop<- ncol(cast_lachno_hypothetical_genes)
custom_config <- umap.defaults
custom_config$n_neighbors = 100
custom_config$metric ='manhattan'
custom_config
comparison_map <- umap(cast_lachno_hypothetical_genes[,2:stop],custom_config)


coord_hypo <- as.data.frame(comparison_map$layout)
coord_hypo$seq_id <- cast_lachno_hypothetical_genes$seq_id

coord_hypo <- merge(coord_hypo,lachnospiraceae_working,by='seq_id')

Figure_S4_A <- ggplot(coord_hypo,aes(x=V1,y=V2))+
  geom_point(aes(fill=species),pch=21,color="black",size=4,alpha=0.7)+
  scale_fill_manual(values=species.colors) +
  theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"))+
  theme(legend.background=element_blank())

Figure_S4_A

#FIGURE S5B Annotated proteins only

set.seed(1)
stop<- ncol(cast_lachno_annotated_genes)
custom_config <- umap.defaults
custom_config$n_neighbors = 100
custom_config$metric ='manhattan'
custom_config
comparison_map <- umap(cast_lachno_annotated_genes[,2:stop],custom_config)


coord_annotated <- as.data.frame(comparison_map$layout)
coord_annotated$seq_id <- cast_lachno_annotated_genes$seq_id

coord_annotated <- merge(coord_annotated,lachnospiraceae_working,by='seq_id')

Figure_S4_B <- ggplot(coord_annotated,aes(x=V1,y=V2))+
  geom_point(aes(fill=species),pch=21,color="black",size=4,alpha=0.7)+
  scale_fill_manual(values=species.colors) +
  theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"))+
  theme(legend.background=element_blank())

Figure_S4_B


# Figure 6 Multiple UMAP vs 16S ----------------------------------------------------------

#Get 16S intra-isolate distances (calculated for Figure S3)
fig_6_melted_16S_distance <- melted_16S_distance
fig_6_melted_16S_distance$id_A <- ifelse(melted_16S_distance$id < melted_16S_distance$variable,
                                         melted_16S_distance$id,
                                         melted_16S_distance$variable)
fig_6_melted_16S_distance$id_B <- ifelse(melted_16S_distance$id > melted_16S_distance$variable,
                                         melted_16S_distance$id,
                                         melted_16S_distance$variable)

fig_6_melted_16S_distance <- fig_6_melted_16S_distance %>%
  select(id_A,id_B,species.y,species.x,value)

#generate intra-isolate distances from single seed UMAP


for(i in 1:1){
  print(i)
  set.seed(i)
  stop<- ncol(cast_lachno_genes)
  custom_config <- umap.defaults
  custom_config$n_neighbors = 100
  custom_config$metric ='manhattan'
  custom_config
  print("Running UMAP")
  comparison_map <- umap(cast_lachno_genes[,2:stop],custom_config)
  print("Finished UMAP")
  
  coord_boot <- as.data.frame(comparison_map$layout)
  coord_boot$seq_id <- cast_lachno_genes$seq_id
  
  coord_boot <- merge(coord_boot,lachnospiraceae_working,by='seq_id')
  
  coord_boot <- coord_boot %>% 
    select(msk_id,V1,V2)
  
  distances_boot <- dist(coord_boot)
  distances_boot <- as.matrix(distances_boot)
  distances_boot[lower.tri(distances_boot)] <- NA 
  row.names(distances_boot) <- coord_boot$msk_id
  colnames(distances_boot) <- c(row.names(distances_boot))
  distances_boot <- as.data.frame(distances_boot)
  distances_boot$id <- row.names(distances_boot)
  
  melted_boot_distance <- melt(distances_boot)
  
  melted_boot_distance$variable <- as.character(melted_boot_distance$variable)
  melted_boot_distance <- melted_boot_distance %>% 
    filter(value!='NA') %>%
    filter(variable != id) %>%
    merge(lachnospiraceae_working,by.x='id',by.y='msk_id') %>%
    select(id,variable,value,species) %>%
    merge(lachnospiraceae_working,by.x="variable",by.y="msk_id") %>%
    select(id,variable,value,species.x,species.y)
  
  melted_boot_distance$species.x <- as.character(melted_boot_distance$species.x)
  melted_boot_distance$species.y <- as.character(melted_boot_distance$species.y)
  
  melted_boot_distance$id_A <- ifelse(melted_boot_distance$id < melted_boot_distance$variable,
                                      melted_boot_distance$id,
                                      melted_boot_distance$variable)
  melted_boot_distance$id_B <- ifelse(melted_boot_distance$id > melted_boot_distance$variable,
                                      melted_boot_distance$id,
                                      melted_boot_distance$variable)
  
  melted_boot_distance <- melted_boot_distance %>%
    select(id_A,id_B,species.y,species.x,value)
  
 master_boot_distance <-melted_boot_distance
  fig_6_melted_single_distance <- melted_boot_distance
}

for(i in 2:100){
  print(i)
  set.seed(i)
  stop<- ncol(cast_lachno_genes)
  custom_config <- umap.defaults
  custom_config$n_neighbors = 100
  custom_config$metric ='manhattan'
  custom_config
  print("Running UMAP")
  comparison_map <- umap(cast_lachno_genes[,2:stop],custom_config)
  print("Finished UMAP")
  
  coord_boot <- as.data.frame(comparison_map$layout)
  coord_boot$seq_id <- cast_lachno_genes$seq_id
  
  coord_boot <- merge(coord_boot,lachnospiraceae_working,by='seq_id')
  
  coord_boot <- coord_boot %>% 
    select(msk_id,V1,V2)
  
  distances_boot <- dist(coord_boot)
  distances_boot <- as.matrix(distances_boot)
  distances_boot[lower.tri(distances_boot)] <- NA 
  row.names(distances_boot) <- coord_boot$msk_id
  colnames(distances_boot) <- c(row.names(distances_boot))
  distances_boot <- as.data.frame(distances_boot)
  distances_boot$id <- row.names(distances_boot)
  
  melted_boot_distance <- melt(distances_boot)
  
  melted_boot_distance$variable <- as.character(melted_boot_distance$variable)
  melted_boot_distance <- melted_boot_distance %>% 
    filter(value!='NA') %>%
    filter(variable != id) %>%
    merge(lachnospiraceae_working,by.x='id',by.y='msk_id') %>%
    select(id,variable,value,species) %>%
    merge(lachnospiraceae_working,by.x="variable",by.y="msk_id") %>%
    select(id,variable,value,species.x,species.y)
  
  melted_boot_distance$species.x <- as.character(melted_boot_distance$species.x)
  melted_boot_distance$species.y <- as.character(melted_boot_distance$species.y)
  
  melted_boot_distance$id_A <- ifelse(melted_boot_distance$id < melted_boot_distance$variable,
                                      melted_boot_distance$id,
                                      melted_boot_distance$variable)
  melted_boot_distance$id_B <- ifelse(melted_boot_distance$id > melted_boot_distance$variable,
                                      melted_boot_distance$id,
                                      melted_boot_distance$variable)
  
  melted_boot_distance <- melted_boot_distance %>%
    select(id_A,id_B,species.y,species.x,value)
  print("Joining")
  master_boot_distance <- inner_join(master_boot_distance,melted_boot_distance,by=c("id_A"="id_A","id_B"="id_B","species.y"="species.y","species.x"="species.x"))
  
}

#write.csv(file="master_boot_distance.csv",master_boot_distance)
master_boot_distance <- read.csv(file="master_boot_distance.csv",header=TRUE,sep=",")
master_boot_distance <- master_boot_distance[,-1]
master_boot_distance$average <- rowMeans(master_boot_distance[,5:104],na.rm=TRUE)

plot_boot_distance <- master_boot_distance %>%
  select(id_A,id_B,species.y,species.x,average)

fig_r3_umap_vs_boot <- inner_join(plot_boot_distance,fig_6_melted_single_distance,by=c("id_A"="id_A","id_B"="id_B"))

Figure_r_3 <- ggplot(fig_r3_umap_vs_boot,aes(x=value,y=average))+
  geom_point(pch=21,fill="#2E9DB7",color="#2E9DB7",size=2,alpha=.1)+
  #geom_smooth(method="lm",se=TRUE,color='black')+
  stat_density_2d(color="black")+
  theme(axis.title = element_blank(), plot.title=element_text(size=6), panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"))+
  theme(legend.position="none")
Figure_r_3

fig_6A_16S_vs_boot <- inner_join(fig_6_melted_16S_distance,plot_boot_distance,by=c("id_A"="id_A","id_B"="id_B"))

Figure_6_A_Boot <- ggplot(fig_6A_16S_vs_boot,aes(x=value,y=average))+
  geom_point(pch=21,fill="#2E9DB7",color="#2E9DB7",size=2,alpha=.1)+
  geom_smooth(method="lm",se=TRUE,color='black')+
  stat_density_2d(color="black")+
  theme(axis.title = element_blank(), plot.title=element_text(size=6), panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"))+
  theme(legend.position="none")

Figure_6_A_Boot

Fig_6A_correlate <- fig_6A_16S_vs_boot %>%
  select(value,average)
cor(Fig_6A_correlate,use="all.obs",method="pearson")


test_species <- unique(as.character(lachnospiraceae_working$species))
test_species[1]
individual_species_distances <- NULL
test_species
for(i in 1:27){
  print(i)
  
  fig_6A_16S_vs_boot$plot_species_test <- ifelse(as.character(fig_6A_16S_vs_boot$species.x.x==test_species[i]),fig_6A_16S_vs_boot$species.y.x,'NA')
  fig_6A_16S_vs_boot$plot_species_test <- ifelse(as.character(fig_6A_16S_vs_boot$species.y.x==test_species[i]),fig_6A_16S_vs_boot$species.x.x,fig_6A_16S_vs_boot$plot_species_test)
  plot_selected <- fig_6A_16S_vs_boot %>% 
    filter(species.x.x==test_species[i] | species.y.x==test_species[i])
 
  individual_species_distances[[i]] <- ggplot(fig_6A_16S_vs_boot,aes(x=value,y=average))+
    geom_point(pch=21,fill="#383838",color="#383838",size=2,alpha=.05)+
    geom_point(data=plot_selected,aes(x=value,y=average,fill=plot_species_test),pch=21,color="black",size=2,alpha=.7)+
    scale_fill_manual(values=species.colors) +
    #stat_density_2d(aes(fill=..level..),geom="polygon",color='white')+
    geom_smooth(method="lm",se=TRUE,color='black')+
    ggtitle(paste(test_species[i]," vs."))+
    theme(axis.title = element_blank(), plot.title=element_text(size=6), panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"),axis.text = element_blank())+
    theme(legend.position="none")
}

do.call(grid.arrange,c(individual_species_distances[c(19,17,7)],ncol=3))


# Figure R3 n_neighbors=15 UMAP  -------------------------------------------


for(i in 1:1){
  print(i)
  set.seed(i)
  stop<- ncol(cast_lachno_genes)
  custom_config <- umap.defaults
  custom_config$n_neighbors = 15
  custom_config$metric ='manhattan'
  custom_config
  print("Running UMAP")
  comparison_map <- umap(cast_lachno_genes[,2:stop],custom_config)
  print("Finished UMAP")
  
  coord_boot <- as.data.frame(comparison_map$layout)
  coord_boot$seq_id <- cast_lachno_genes$seq_id
  
  coord_boot <- merge(coord_boot,lachnospiraceae_working,by='seq_id')
  
  coord_boot <- coord_boot %>% 
    select(msk_id,V1,V2)
  
  distances_boot <- dist(coord_boot)
  distances_boot <- as.matrix(distances_boot)
  distances_boot[lower.tri(distances_boot)] <- NA 
  row.names(distances_boot) <- coord_boot$msk_id
  colnames(distances_boot) <- c(row.names(distances_boot))
  distances_boot <- as.data.frame(distances_boot)
  distances_boot$id <- row.names(distances_boot)
  
  melted_boot_distance <- melt(distances_boot)
  
  melted_boot_distance$variable <- as.character(melted_boot_distance$variable)
  melted_boot_distance <- melted_boot_distance %>% 
    filter(value!='NA') %>%
    filter(variable != id) %>%
    merge(lachnospiraceae_working,by.x='id',by.y='msk_id') %>%
    select(id,variable,value,species) %>%
    merge(lachnospiraceae_working,by.x="variable",by.y="msk_id") %>%
    select(id,variable,value,species.x,species.y)
  
  melted_boot_distance$species.x <- as.character(melted_boot_distance$species.x)
  melted_boot_distance$species.y <- as.character(melted_boot_distance$species.y)
  
  melted_boot_distance$id_A <- ifelse(melted_boot_distance$id < melted_boot_distance$variable,
                                      melted_boot_distance$id,
                                      melted_boot_distance$variable)
  melted_boot_distance$id_B <- ifelse(melted_boot_distance$id > melted_boot_distance$variable,
                                      melted_boot_distance$id,
                                      melted_boot_distance$variable)
  
  melted_boot_distance <- melted_boot_distance %>%
    select(id_A,id_B,species.y,species.x,value)
  
  master_boot_n10_distance <-melted_boot_distance
  single_seed_n10_distance <- melted_boot_distance
}

for(i in 2:100){
  print(i)
  set.seed(i)
  stop<- ncol(cast_lachno_genes)
  custom_config <- umap.defaults
  custom_config$n_neighbors = 15
  custom_config$metric ='manhattan'
  custom_config
  print("Running UMAP")
  comparison_map <- umap(cast_lachno_genes[,2:stop],custom_config)
  print("Finished UMAP")
  
  coord_boot <- as.data.frame(comparison_map$layout)
  coord_boot$seq_id <- cast_lachno_genes$seq_id
  
  coord_boot <- merge(coord_boot,lachnospiraceae_working,by='seq_id')
  
  coord_boot <- coord_boot %>% 
    select(msk_id,V1,V2)
  
  distances_boot <- dist(coord_boot)
  distances_boot <- as.matrix(distances_boot)
  distances_boot[lower.tri(distances_boot)] <- NA 
  row.names(distances_boot) <- coord_boot$msk_id
  colnames(distances_boot) <- c(row.names(distances_boot))
  distances_boot <- as.data.frame(distances_boot)
  distances_boot$id <- row.names(distances_boot)
  
  melted_boot_distance <- melt(distances_boot)
  
  melted_boot_distance$variable <- as.character(melted_boot_distance$variable)
  melted_boot_distance <- melted_boot_distance %>% 
    filter(value!='NA') %>%
    filter(variable != id) %>%
    merge(lachnospiraceae_working,by.x='id',by.y='msk_id') %>%
    select(id,variable,value,species) %>%
    merge(lachnospiraceae_working,by.x="variable",by.y="msk_id") %>%
    select(id,variable,value,species.x,species.y)
  
  melted_boot_distance$species.x <- as.character(melted_boot_distance$species.x)
  melted_boot_distance$species.y <- as.character(melted_boot_distance$species.y)
  
  melted_boot_distance$id_A <- ifelse(melted_boot_distance$id < melted_boot_distance$variable,
                                      melted_boot_distance$id,
                                      melted_boot_distance$variable)
  melted_boot_distance$id_B <- ifelse(melted_boot_distance$id > melted_boot_distance$variable,
                                      melted_boot_distance$id,
                                      melted_boot_distance$variable)
  
  melted_boot_distance <- melted_boot_distance %>%
    select(id_A,id_B,species.y,species.x,value)
  print("Joining")
  master_boot_n10_distance <- inner_join(master_boot_n10_distance,melted_boot_distance,by=c("id_A"="id_A","id_B"="id_B","species.y"="species.y","species.x"="species.x"))
  
}

master_boot_n10_distance$average <- rowMeans(master_boot_n10_distance[,5:104],na.rm=TRUE)

plot_n10_boot_distance <- master_boot_n10_distance %>%
  select(id_A,id_B,species.y,species.x,average)

fig_r3b_umap_vs_boot_n10 <- inner_join(plot_n10_boot_distance,single_seed_n10_distance,by=c("id_A"="id_A","id_B"="id_B"))

Figure_r_3b <- ggplot(fig_r3b_umap_vs_boot_n10,aes(x=value,y=average))+
  geom_point(pch=21,fill="#2E9DB7",color="#2E9DB7",size=2,alpha=.1)+
  #geom_smooth(method="lm",se=TRUE,color='black')+
  stat_density_2d(color="black")+
  theme(axis.title = element_blank(), plot.title=element_text(size=6), panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"))+
  theme(legend.position="none")+
  ylim(0,40)+
  xlim(0,40)
Figure_r_3b




# Figure 7A UMAP individual species ---------------------------------------------------------------

species_to_analyze <- lachnospiraceae_working %>%
  group_by(species) %>%
  dplyr::summarise(n()) %>%
  ungroup(species) %>%
  dplyr::rename("count"='n()') %>%
  filter(count >9)

next_lachno_species <- as.character(species_to_analyze$species)

next_lachno <- NULL
r_gnavus_coor <- NULL

for(i in 1:8){
  print(i)
  print("analyzing")
  print(next_lachno_species[i])
  
  lachno_umap <- lachnospiraceae_working %>% 
    filter(species %in% next_lachno_species[i])
  #castpidshypo <- cast(pidssort4,seq_id~cluster_num)
  #castpidshypo2 <- mutate_at( castpidshypo, vars(-one_of("seq_id")),
  #   list(~ ifelse( . > 1, 1, .)))
  
  castpids_ind <- cast_lachno_genes %>% 
    filter(seq_id %in% lachno_umap$seq_id)
  
  set.seed(1)
  stop<- ncol(castpids_ind)
  custom_config <- umap.defaults
  custom_config$n_neighbors = 10
  custom_config$metric ='manhattan'
  custom_config
  comparison_map_ind <- umap(castpids_ind[,2:stop],custom_config)
  
  
  coord_ind <- as.data.frame(comparison_map_ind$layout)
  coord_ind$seq_id <- castpids_ind$seq_id
  
  coord_ind <- merge(coord_ind,lachnospiraceae_working,by='seq_id')
  
  r_gnavus_coor[[i]] <- coord_ind

  next_lachno[[i]]<- ggplot(coord_ind,aes(x=V1,y=V2))+
    geom_point(aes(fill=species),pch=21, color="black", size=4,alpha=0.7)+
    scale_fill_manual(values=species.colors)+
    theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"),
          legend.background=element_blank(),legend.position = "none",axis.title = element_blank())+
    xlim(-25,25)+
    ylim(-17.5,17.5)
  
  next_lachno_species[i]=='Blautia wexlerae'
  
}

r_gnavus_coor[2]
Figure_7C <- next_lachno[2]
Figure_7B <- next_lachno[5]

Figure_7A <- do.call(grid.arrange,c(next_lachno,ncol=4))


# Figure_S5 Individual Species UMAP Hypothetical/Annotated ----------------

species_to_analyze <- lachnospiraceae_working %>%
  group_by(species) %>%
  dplyr::summarise(n()) %>%
  ungroup(species) %>%
  dplyr::rename("count"='n()') %>%
  filter(count >9)

next_lachno_species <- as.character(species_to_analyze$species)

lachno_ind_annotate <- NULL
for(i in 1:8){
  print(i)
  print("analyzing")
  print(next_lachno_species[i])
  
  lachno_umap <- lachnospiraceae_working %>% 
    filter(species %in% next_lachno_species[i])

  castpids_ind_annotate <- cast_lachno_annotated_genes %>% 
    filter(seq_id %in% lachno_umap$seq_id)
  
  set.seed(1)
  stop<- ncol(castpids_ind_annotate)
  custom_config <- umap.defaults
  custom_config$n_neighbors = 10
  custom_config$metric ='manhattan'
  custom_config
  comparison_map_ind_annotate <- umap(castpids_ind_annotate[,2:stop],custom_config)
  
  
  coord_ind_annotate <- as.data.frame(comparison_map_ind_annotate$layout)
  coord_ind_annotate$seq_id <- castpids_ind_annotate$seq_id
  
  coord_ind_annotate <- merge(coord_ind_annotate,lachnospiraceae_working,by='seq_id')
  
  
  lachno_ind_annotate[[i]]<- ggplot(coord_ind_annotate,aes(x=V1,y=V2))+
    geom_point(aes(fill=species),pch=21, color="black", size=4,alpha=0.7)+
    scale_fill_manual(values=species.colors)+
    theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"),
          legend.background=element_blank(),legend.position = "none",axis.title = element_blank())+
    xlim(-25,25)+
    ylim(-17.5,17.5)
}

Figure_S5A <- do.call(grid.arrange,c(lachno_ind_annotate,ncol=4))


lachno_ind_hypo <- NULL
for(i in 1:8){
  print(i)
  print("analyzing")
  print(next_lachno_species[i])
  
  lachno_umap <- lachnospiraceae_working %>% 
    filter(species %in% next_lachno_species[i])
  
  castpids_ind_hypo <- cast_lachno_hypothetical_genes %>% 
    filter(seq_id %in% lachno_umap$seq_id)
  
  set.seed(1)
  stop<- ncol(castpids_ind_hypo)
  custom_config <- umap.defaults
  custom_config$n_neighbors = 10
  custom_config$metric ='manhattan'
  custom_config
  comparison_map_ind_hypo <- umap(castpids_ind_hypo[,2:stop],custom_config)
  
  
  coord_ind_hypo <- as.data.frame(comparison_map_ind_hypo$layout)
  coord_ind_hypo$seq_id <- castpids_ind_hypo$seq_id
  
  coord_ind_hypo <- merge(coord_ind_hypo,lachnospiraceae_working,by='seq_id')
  
  
  lachno_ind_hypo[[i]]<- ggplot(coord_ind_hypo,aes(x=V1,y=V2))+
    geom_point(aes(fill=species),pch=21, color="black", size=4,alpha=0.7)+
    scale_fill_manual(values=species.colors)+
    theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"),
          legend.background=element_blank(),legend.position = "none",axis.title = element_blank())+
    xlim(-25,25)+
    ylim(-17.5,17.5)
}

Figure_S5b <- do.call(grid.arrange,c(lachno_ind_hypo,ncol=4))








# Figure 7B ---------------------------------------------------------------
blautia_wexlerae <- lachnospiraceae_working %>% 
  filter(species == 'Blautia wexlerae')
r_gnavus <- lachnospiraceae_working %>% 
  filter(species =='[Ruminococcus] gnavus')

#read in cluster assignments

blautia_wexlerae_cluster <- read.csv(file="blautia_wexlerae_cluster_assignment.csv",header=TRUE,sep=",")%>%
  inner_join(blautia_wexlerae,by="seq_id")
r_gnavus_cluster <- read.csv(file="r_gnavus_cluster_assignment.csv",header=TRUE,sep=",") %>%
  inner_join(r_gnavus,by="seq_id")


#assess core genome by cluster

bw_cluster_master <- data.frame("bw_totalcountscluster"=numeric(),"core_in_cluster"=logical(), "gene"=character(),"cluster"=integer()) 
rg_cluster_master <- data.frame("bw_totalcountscluster"=numeric(),"core_in_cluster"=logical(), "gene"=character(),"cluster"=integer()) 

for(i in 1:3){
  print(paste("printing by ",i))
  bw_clusterfilter <- blautia_wexlerae_cluster %>% 
    filter(cluster==i)
  rg_clusterfilter <- r_gnavus_cluster %>% 
    filter(cluster==i)
  
  bw_cluster_genes <- cast_lachno_genes %>% 
    filter(seq_id %in% bw_clusterfilter$seq_id)
  rg_cluster_genes <- cast_lachno_genes %>% 
    filter(seq_id %in% rg_clusterfilter$seq_id)
  bw_cluster_gene_names <- colnames(bw_cluster_genes[,-1])
  bw_totalcountscluster <- colSums(bw_cluster_genes[,-1])
  
  
  rg_totalcountscluster <- colSums(rg_cluster_genes[,-1])
  rg_cluster_gene_names <- colnames(rg_cluster_genes[,-1])
  
  bw_cluster <- as.data.frame(bw_totalcountscluster) %>%
    dplyr::mutate(core_in_cluster=bw_totalcountscluster == (nrow(bw_cluster_genes)))%>%
    dplyr::mutate(gene=bw_cluster_gene_names) %>%
    dplyr::mutate(cluster=i)
    
  rg_cluster <- as.data.frame(rg_totalcountscluster) %>%
    dplyr::mutate(core_in_cluster=rg_totalcountscluster == (nrow(rg_cluster_genes)))%>%
    dplyr::mutate(gene=rg_cluster_gene_names) %>%
    dplyr::mutate(cluster=i)
  
  
  bw_cluster_master <- rbind(bw_cluster_master,bw_cluster)
  rg_cluster_master <- rbind(rg_cluster_master,rg_cluster)
  
}


Figure_7B_Counts <- bw_cluster_master %>%
  filter(core_in_cluster==TRUE) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(num_clusters = n())%>%
  filter(num_clusters==1) %>%
  dplyr::mutate(hypothetical=!is.na(as.numeric(gene)))

Figure_7B_table <- as.data.frame(unique(blautia_wexlerae_cluster$cluster))
colnames(Figure_7B_table) <- "cluster"

Figure_7B_number_ID_Products <- Figure_7B_Counts %>%
  filter(hypothetical==FALSE) %>%
  dplyr::group_by(cluster)%>%
  dplyr::summarise(n())

Figure_7B_number_ID_clusters <- Figure_7B_Counts %>%
  filter(gene %in% identified_cluster_list)%>%
  dplyr::group_by(cluster)%>%
  dplyr::summarise(n())

Figure_7B_number_hypothetical <- Figure_7B_Counts %>%
  filter(hypothetical==TRUE)%>%
  filter(gene %!in% identified_cluster_list)%>%
  dplyr::group_by(cluster)%>%
  dplyr::summarise(n())

Figure_7B_table <- Figure_7B_table %>%
  merge(Figure_7B_number_ID_Products,by="cluster",all.x=TRUE) %>%
  merge(Figure_7B_number_ID_clusters,by="cluster",all.x=TRUE)%>%
  merge(Figure_7B_number_hypothetical,by="cluster",all.x=TRUE)

colnames(Figure_7B_table) <- c("cluster","ID_products","ID_Clusters","hypothetical")


Figure_7C_table <- as.data.frame(unique(r_gnavus_cluster$cluster))
colnames(Figure_7C_table) <- "cluster"
Figure_7C_Counts <- rg_cluster_master %>%
  filter(core_in_cluster==TRUE) %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(num_clusters = n())%>%
  filter(num_clusters==1) %>%
  dplyr::mutate(hypothetical=!is.na(as.numeric(gene)))


Figure_7C_number_ID_Products <- Figure_7C_Counts %>%
  filter(hypothetical==FALSE) %>%
  dplyr::group_by(cluster)%>%
  dplyr::summarise(n())

Figure_7C_number_ID_clusters <- Figure_7C_Counts %>%
  filter(gene %in% identified_cluster_list)%>%
  dplyr::group_by(cluster)%>%
  dplyr::summarise(n())

Figure_7C_number_hypothetical <- Figure_7C_Counts %>%
  filter(hypothetical==TRUE)%>%
  filter(gene %!in% identified_cluster_list)%>%
  dplyr::group_by(cluster)%>%
  dplyr::summarise(n())

Figure_7C_table <- Figure_7C_table %>%
  merge(Figure_7C_number_ID_Products,by="cluster",all.x=TRUE) %>%
  merge(Figure_7C_number_ID_clusters,by="cluster",all.x=TRUE)%>%
  merge(Figure_7C_number_hypothetical,by="cluster",all.x=TRUE)
colnames(Figure_7C_table) <- c("cluster","ID_products","ID_Clusters","hypothetical")



# Figure 7D Glucorhamnan coverage -----------------------------------------

row.names(polysaccharide_coverage) <- polysaccharide_coverage[,1]
polysaccharide_coverage <- polysaccharide_coverage[,-1]

polysaccharide_coverage <- polysaccharide_coverage %>%
  t()
polysaccharide_coverage <- as.data.frame(polysaccharide_coverage)

polysaccharide_coverage$seq_id <- row.names(polysaccharide_coverage)

polysaccharide_coverage[is.na(polysaccharide_coverage)] <- 0

#glimpse(polysaccharide_coverage)
summary_coverage <- polysaccharide_coverage %>%
  gather(key="position",value="reads",-seq_id,) %>%
  inner_join(polysaccharide_groups,by="seq_id") %>%
  select(-seq_id)%>%
  dplyr::group_by(coverage_group,position) %>%
  dplyr::summarise_all(funs(mean(reads),sd(reads)))%>%
  dplyr::mutate(position=as.numeric(position))%>%
  arrange(position) %>%
  dplyr::group_by(coverage_group) %>%
  dplyr::mutate(bin=cut(position,600)) %>%
  select(bin,coverage_group,mean,sd)%>%
  dplyr::group_by(bin,coverage_group) %>%
  dplyr::summarise_all(funs(mean))%>%
  dplyr::mutate(upper=mean+sd)%>%
  dplyr::mutate(lower=mean-sd) %>%
  dplyr::group_by(coverage_group) %>%
  dplyr::mutate(position=50.185*seq(1:n()))

colnames(summary_coverage) <- c("bin","coverage_group","group_mean","group_sd","upper","lower","position")
#head(summary_coverage)
Figure_7D_Bottom <- summary_coverage %>%
  ggplot(aes(x=position,y=group_mean))+
    facet_wrap(~coverage_group,ncol=1)+  
      geom_line(aes(color=coverage_group))+
      geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.3)+
      scale_color_manual(values=c("red","green","blue"))+
  theme_minimal()+
  theme(panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"),
        panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        legend.position="none",
        strip.text = element_blank())
  
Figure_7D_Top <-   summary_coverage %>%
  ggplot(aes(x=position,y=group_mean))+
     
    geom_line(aes(color=coverage_group))+
  
    scale_color_manual(values=c("red","green","blue"))+
    theme_minimal()+
    theme(panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"),
        panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        legend.position="none",
        strip.text = element_blank())
figure_7E_df <- r_gnavus_coor[[2]]
figure_7E_df <- merge(figure_7E_df,polysaccharide_groups,by="seq_id")
Figure_7E <- ggplot(figure_7E_df,aes(x=V1,y=V2))+
  geom_point(aes(fill= coverage_group),pch=21, color="black", size=4,alpha=0.7)+
  #scale_fill_manual(values=species.colors)+
  theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#000000",size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#a8a8a8"),
        legend.background=element_blank(),axis.title = element_blank())+
  xlim(-25,25)+
  ylim(-17.5,17.5)
Figure_7E


# Figure S6 Lachnospiraceae  ----------------------------------------------
gnavus_colors <- r_gnavus_cluster

gnavus_pH <- lachnospiraceae_working %>%
  filter(seq_id %in% gnavus_colors$seq_id) %>%
  select(msk_id,seq_id) %>%
  inner_join(lachno_pH,by="msk_id") %>%
  inner_join(gnavus_colors,by="seq_id")

r_gnavus_cluster

gnavus_pH$cluster <- as.factor(gnavus_pH$cluster)
gnavus_pH$cluster <- factor(gnavus_pH$cluster,levels=c("1","3","2"))


Figure_S6 <- ggplot(gnavus_pH,aes(x=cluster,y=pH))+
  geom_bar(stat="summary",fun.y="mean",fill="#383838",alpha=.7,color="black")+
  geom_jitter(aes(fill=species),pch=21,size=3,width=.2)+
  scale_fill_manual(values=species.colors)+
  theme_minimal()+
  theme(legend.position = "none")+
  coord_cartesian(ylim=c(5.85,6.5))+
  scale_x_discrete(labels=c("1" = "A", "2" = "C",
                            "3" = "B"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5))
Figure_S6
anova_gnavus_pH <- aov(pH ~as.factor(cluster),data=gnavus_pH)
summary(anova_gnavus_pH)
TukeyHSD(anova_gnavus_pH)
