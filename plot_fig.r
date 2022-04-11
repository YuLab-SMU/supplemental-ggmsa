library(ggmsa)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(gggenes)
library(ape)
library(phangorn)
library(Biostrings)
library(ggnewscale)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(ggplotify)
library(aplot)
library(magick)

protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")
nt_sequence <- system.file("extdata", "LeaderRepeat_All.fa", package = "ggmsa")
miRNA_sequences <- system.file("extdata", "seedSample.fa", package = "ggmsa")
tp53_sequences <-  system.file("extdata", "tp53.fa", package = "ggmsa")
tp53_genes <- system.file("extdata", "TP53_genes.xlsx", package = "ggmsa")

##Fig.2A sequence logo + sequence bundles##
p2A <- seqlogo("data/Gram-NP-merge.fa", 
               color = "Chemistry_AA", 
               font = "DroidSansMono") + 
  coord_cartesian()


negative <-  system.file("extdata", "Gram-negative_AKL.fasta", 
                         package = "ggmsa")
positive <-  system.file("extdata", "Gram-positive_AKL.fasta", 
                         package = "ggmsa")

pos <- data.frame(x= c(4, 7, 9, 24, 27, 29,
                       4, 7, 24, 27), 
                  y = c(c(21, 11, 20, 17, 12, 18) +.3,
                        c(13, 13, 13, 13) +.5
                  ), 
                  label = c("H", "S", "R", "D", "T", "E",
                            "C", "C", "C", "C"),
                  color = c(rep("#ff4700",6), 
                            rep("#0443d0",4)))

p2A2 <- ggSeqBundle(list(negative, positive), 
                    alpha = 0.1, 
                    bundle_color = c("#FC8D62","#8DA0CB"))+ #RColorBrewer: Set2:2-3
  geom_text(data = pos, 
            mapping = aes(x, y, 
                          label = label, 
                          color = I(color)),
            inherit.aes = FALSE, 
            size = 4)

##Fig.2B stacked MSA + logo&bar annotations##
p2B <- ggmsa(protein_sequences, 221, 280, seq_name = TRUE, char_width = 0.5,border = "white") +
  geom_seqlogo(color = "Chemistry_AA") +
  geom_msaBar()

##Fig.2C MSA + Arc##
p_RF03120_SS <- image_read_svg("data/RF03120_SS.svg")
q_RF03120_SS <- as.ggplot(p_RF03120_SS)

RF03120_msa<- system.file("extdata", "Rfam", "RF03120.fasta", package = "ggmsa")
RF03120_ss <- system.file("extdata", "Rfam", "RF03120_SS.txt", package = "ggmsa")
RF_arc <- readSSfile(RF03120_ss, type = "Vienna" )
p2c <- ggmsa(RF03120_msa,
             font = NULL,
             color = "Chemistry_NT",
             seq_name = F,
             show.legend = F,
             border = NA) +
  geom_helix(helix_data = RF_arc) +
  theme(axis.text.y = element_blank())
require(patchwork)
p2C <- p2c + inset_element(q_RF03120_SS, 
                           left = 0, 
                           bottom = 0.6,
                           right = 0.4, 
                           top = 1, 
                           align_to = 'full') 
##Fig.2D 2Arc##
tpp_seq <- "data/tpp_riboswitch.fasta"
arc_4NYG <- "data/riboswitch_thiamine.txt"
arc_4NYD <- "data/riboswitch_hypoxanthine.txt"
thiamine <- readSSfile(arc_4NYG, type = "Vienna" )
hypoxanthine <- readSSfile(arc_4NYD, type = "Vienna")
p_double_arc <- ggmsa(tpp_seq, 
                      color = "Chemistry_NT", 
                      seq_name = F, 
                      show.legend = F, 
                      border = NA) +
  geom_helix(helix_data = list(known = hypoxanthine, 
                               predicted = thiamine)) + 
  theme(axis.text.y = element_blank())
p_loop1 <- image_read_pdf(path = "data/bpRNA_PDB_590_ColorCodedStructures_4NYG.pdf",
                          density = 300)
p_loop2 <- image_read_pdf(path = "data/bpRNA_PDB_589_ColorCodedStructures_4NYD.pdf",
                          density = 300)
q1 <- as.ggplot(p_loop1)
q2 <- as.ggplot(p_loop2)
p_loop <- plot_list(gglist = list(q1, q2), ncol = 1, labels = c("D"," "))
p2D <- plot_list(gglist = list(p_loop, p_double_arc), ncol = 2)

#merge C&D
p2CD <- plot_list(gglist = list(p2C, p2D), 
                  ncol = 2,
                  widths = c(0.6,0.4))

##plot all##
pp <- plot_list(gglist = list(p2A, p2A2,p2B, p2CD), ncol = 1, 
                heights = c(0.3, 1, 0.6, 0.8),
                labels = c("A", " ", " B", "C"))

ggsave("Fig2.pdf", plot = pp, width = 12, height = 11)
ggsave("Fig2.png", plot = pp, width = 12, height = 11)



################################################################################
##Fig 3 sequence recombination##
fas <- c("data/HM_KP.fa","data/CK_KP.fa")
xx <- lapply(fas, seqdiff)
plts <- lapply(xx, plot, width = 100)
plts[[3]] <- simplot("data/CK_HM_KP.fa", 'KP827649',smooth = F) + theme(legend.position = "bottom")

p3 <- plot_list(gglist=plts, ncol=1, tag_levels = list(c("A",' ',"B", ' ',"C")))


ggsave("Fig3.pdf", plot = p3, width = 10, height = 14)
ggsave("Fig3.png", plot = p3, width = 10, height = 14)


################################################################################
##Fig 4 MAF plot
maf <- "data/chr1_KI270707v1_random.txt.maf"
ref = "hg38.chr1_KI270707v1_random"
seq_df <- read_maf(maf)
tidy_df <- tidy_maf_df(seq_df, ref = ref)

p4 <- ggmaf(data = tidy_df, 
            ref = ref, 
            block_start = 1, 
            block_end = 10, 
            facet_field = 5,
            facet_heights = c(0.45,0.55))

ggsave("Fig4.png",p5, width = 12, height = 10)
ggsave("Fig4.pdf",p5, width = 12, height = 10)
################################################################################

##Fig 5 graphics combination##

##Fig 5 tree + msa + genes locus
dat <- read.aa(tp53_sequences, format = "fasta") %>% phyDat(type = "AA", levels = NULL)
tree <- dist.ml(dat, model = "JTT") %>% bionj()
dd <- ggimage::phylopic_uid(tree$tip.label)

p_tp53 <- ggtree(tree, branch.length = 'none') %<+% dd +
  geom_tiplab(aes(image=uid), geom = "phylopic", offset =1.9) +
  geom_tiplab(aes(label=label)) +
  geom_treescale(x = 0,y = -1)
#msa
data_53 <- readAAMultipleAlignment(tp53_sequences) %>% tidy_msa()
#gene maps
TP53_arrow <- readxl::read_xlsx(tp53_genes)
TP53_arrow$direction <- 1
TP53_arrow[TP53_arrow$strand == "reverse","direction"] <- -1

#color
mapping = aes(xmin = start, xmax = end, fill = gene, forward = direction)
my_pal <- colorRampPalette(rev(brewer.pal(n = 10, name = "Set3")))

#tree + gene maps + msa
p5a <- p_tp53  + xlim_tree(4) +
  geom_facet(geom = geom_msa, data = data_53,
             panel = 'Multiple Sequence Alignment of the TP53 Protein', font = NULL,
             border = NA) +
  new_scale_fill() +
  scale_fill_manual(values = my_pal(10)) +
  geom_facet(geom = geom_motif,
             mapping = mapping, data = TP53_arrow,
             panel = 'Genome_Locus',  on = 'TP53',
             arrowhead_height = unit(3, "mm"),
             arrowhead_width = unit(1, "mm")) +
  theme(strip.background=element_blank(),
        strip.text = element_text(size = 13))


p5A <- facet_widths(p5a, c(Tree = 0.35, Genome_Locus = 0.3))

# ggsave("Fig4A.pdf", width = 13, height = 4)
# ggsave("Fig4A.png", width = 13, height = 4)

##Fig 5B tree + msa + 2boxplot
seq <-readDNAStringSet("data//btuR.fa")
aln <- tidy_msa(seq)
btuR_tree <- read.tree("data/btuR.nwk")
meta_dat <- read.csv("data/meta_data_47.csv")

#Pathotype_fill_colors
Pathotype_cols <- RColorBrewer::brewer.pal(7, "Set3")
names(Pathotype_cols) <- meta_dat$Pathotypes %>% factor %>% levels

###tree OTU
Phylo_group <- list(A= meta_dat$Lineage[meta_dat$Phylogroup == "A"]%>% unique,
                    B1=meta_dat$Lineage[meta_dat$Phylogroup == "B1"]%>% unique,
                    B2=meta_dat$Lineage[meta_dat$Phylogroup == "B2"]%>% unique,
                    C=meta_dat$Lineage[meta_dat$Phylogroup == "C"]%>% unique,
                    D =meta_dat$Lineage[meta_dat$Phylogroup == "D"]%>% unique,
                    E =meta_dat$Lineage[meta_dat$Phylogroup == "E"]%>% unique,
                    `F`=meta_dat$Lineage[meta_dat$Phylogroup == "F"]%>% unique,
                    Shigella=meta_dat$Lineage[meta_dat$Phylogroup == "Shigella"]%>% unique)

Phylo_cols <- RColorBrewer::brewer.pal(8, "Dark2")
names(Phylo_cols) <- names(Phylo_group)

## plot tree
p_btuR_tree <- ggtree(btuR_tree) + geom_tiplab(align = T)
p_btuR_tree <- groupOTU(p_btuR_tree ,Phylo_group)+aes(color=group)  +
  scale_color_manual(values = c(Phylo_cols, "black"), na.value = "black", name = "Lineage",
                     breaks = c("A", "B1", "B2", "C", "D", "E", "F", "Shigella"),  guide="none")

p_btuR_tree <- p_btuR_tree +
  geom_strip('L29', 'L20', barsize=2, color=Phylo_cols[["B2"]],
             label="B2",offset =.01,  offset.text = 0.0015) +
  geom_strip('L28','L29', barsize=2, color=Phylo_cols[["A"]],
             label="A",offset =.01,  offset.text = 0.0015) +
  geom_strip('L15','L28', barsize=2, color=Phylo_cols[["B1"]],
             label="B1",offset =.01,  offset.text = 0.0015) +
  geom_strip('L45','L15', barsize=2, color=Phylo_cols[["Shigella"]],
             label="S.",offset =.01,  offset.text = 0.0015) +
  geom_strip('L36','L45', barsize=2, color=Phylo_cols[["B1"]],
             label="B1",offset =.01,  offset.text = 0.0015) +
  geom_strip('L30','L36', barsize=2, color=Phylo_cols[["Shigella"]],
             label="S.",offset =.01,  offset.text = 0.0015) +
  geom_strip('L39','L30', barsize=2, color=Phylo_cols[["B1"]],
             label="B1",offset =.01,  offset.text = 0.0015) +
  geom_strip('L40','L39', barsize=2, color=Phylo_cols[["C"]],
             label="C",offset =.01,  offset.text = 0.0015)  +
  geom_strip('L48','L40', barsize=2, color=Phylo_cols[["B1"]],
             label="B1",offset =.01,  offset.text = 0.0015) +
  geom_strip('L10','L48', barsize=2, color=Phylo_cols[["E"]],
             label="E",offset =.01,  offset.text = 0.0015) +
  geom_strip('L37','L10', barsize=2, color=Phylo_cols[["D"]],
             label="D",offset =.01,  offset.text = 0.0015)+
  geom_strip('L33','L37', barsize=2, color=Phylo_cols[["F"]],
             label="F",offset =.01,  offset.text = 0.0015)+
  geom_strip('L1','L33', barsize=2, color=Phylo_cols[["E"]],
             label="E",offset =.01,  offset.text = 0.0015)

##tree + meta data boxplots
p5B <- p_btuR_tree  +
  geom_treescale(x = 0,y = -1) +
  geom_fruit(data = aln,
             geom = geom_msa,
             end = 200,
             font = NULL,
             color = "Chemistry_NT",
             border = NA,
             consensus_views = T, 
             ref = "L38",
             pwidth = 3.5,
             offset = 0.3,
             axis.params = list(title = "Multiple Sequence Alignment of the btuR Gene",
                                title.height = 0.05,
                                title.size = 4.5,
                                axis = "x",
                                vjust = 1.1,
                                text.size = 3,
                                line.size = 1,
                                line.color = "black")) +
  new_scale_fill() +
  geom_fruit(mapping = aes(x = AMR_genes, y = Lineage, fill = MDR),
             data = meta_dat,
             geom = geom_boxplot,
             outlier.size = 0.5,
             pwidth=1,
             offset = 0.1,
             axis.params = list(title = "Antimicrobial Classes",
                                title.height = 0.05,
                                title.size = 4.5,
                                axis = "x",
                                vjust = 1.1,
                                text.size = 3,
                                line.size = 1,
                                line.color = "black")) +
  scale_fill_manual(values=c("Yes" = "#fcd7c5", "No" = "#dfdfdf")) +
  new_scale_fill() +
  geom_fruit(mapping = aes(x = virulence_genes, y = Lineage, fill = Pathotypes),
             data = meta_dat,
             geom = geom_boxplot,
             pwidth=1,
             offset = 0.05,
             axis.params = list(title = "Virulence Genes",
                                title.height = 0.05,
                                title.size = 4.5,
                                axis = "x",
                                vjust = 1.1,
                                text.size = 3,
                                line.size = 1,
                                line.color = "black"))+
  scale_fill_manual(values = Pathotype_cols)


p5 <- plot_list(gglist = list(p5A, p5B), ncol = 1, heights = c(0.3,0.7), tag_levels = 'A') 

ggsave("Fig5.png",p5, width = 18, height = 13)
ggsave("Fig5.pdf",p5, width = 18, height = 13)






