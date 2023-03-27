#Locurilla
library(tidyverse)
library(vegan)
library(clustsig)
library(Biostrings)
library(DECIPHER)
library(phyloseq)
library(ape)
library(dendextend)
library(heatmaply)
library(viridis)
library(htmlwidgets)
library(taxa)
library(metacoder)
library(ggvegan)
library(doParallel)
library(hrbrthemes)

data <- read.csv("D:/Documents/GitHub/CORE/Camilo/Analisis/PPE/HCS_PPE_Taxonomic_And_Abundance.csv", header = TRUE, sep = ";", dec = ".", skip = 0, row.names = 1)
environ <- read.csv("~/Documents/GitHub/CORE/Camilo/Analisis/PPE/Shared/HCS_Metal_Environment_NMDS.csv", header = TRUE, sep = ";", dec = ".", skip = 0, row.names = 1)
taxa <- read.csv("~/Documents/GitHub/CORE/Camilo/Analisis/PPE/HCS_PPE_TAXA.csv", header = TRUE, sep = ";", skip = 0, row.names = 1)

##SIMPROF
simp <- data[,c(3:65)]
rownames(simp) <- data$ASV
simp.hell <- decostand(as.data.frame(t(simp)), method = "hellinger")
simp.hc <- vegdist(simp.hell, method = "bray", diag = TRUE) %>% as.matrix() %>% as.data.frame()
sprof <- simprof(data = simp.hc, 
                 method.cluster = "average", 
                 silent = FALSE, 
                 increment = 100, 
                 undef.zero = TRUE)
splot <- simprof.plot(sprof)
sclust <- sprof$hclust
smpr <- add_column(simp.hell, 
                   .before = simp.hell[,1], 
                   Site = c(rep("Chanaral", 12), 
                            rep("Flamenco", 12), 
                            rep("Huasco", 12), 
                            rep("Pta.Choros", 5), 
                            rep("Quintero", 10), 
                            rep("LasCruces", 12)
                            ), 
                   .name_repair = "minimal")
smpr$Site <- factor(smpr$Site, levels = unique(smpr$Site))


##FASTA
nms <- paste(data$ASV, data$Species, sep = "_")
fasta <- data.frame(names = nms, sequences = all_of(data$Seq))
fastados <- data.frame(names = data$ASV, sequences = all_of(data$Seq))
seqRFLP::dataframe2fas(fasta, "D:/Documents/GitHub/CORE/Camilo/Otros/LOCURILLA/ppe_asv.fasta")
seqRFLP::dataframe2fas(fastados, "./Camilo/Otros/LOCURILLA/ppe_asv_dos.fasta")

##Multiple alignment with MUSCLE5 Super5 algorythm (muscle.exe -super5 ppe_asv.fasta -output ppe.afa)

##Alignment dendrogram
multialign <- "./Camilo/Otros/LOCURILLA/ppe.afa"
multidos <- "./Camilo/Otros/LOCURILLA/ppe_dos.afa"
dbConn <-dbConnect(SQLite(), ":memory:")
Seqs2DB(multialign, type = "FASTA", dbFile = dbConn, "")
Seqs2DB(multidos, type = "FASTA", dbFile = dbConn, "")
x <- dbGetQuery(dbConn, "select description from Seqs")$description
Add2DB(myData = data.frame(identifier = x, stringsAsFactors = FALSE), dbConn)
consensus <- IdConsensus(dbConn, threshold = 0.3, minInformation = 0.1)
distance.matrix <- DistanceMatrix(consensus, correction = "Jukes-Cantor", processors = NULL, verbose = TRUE)
dendrogram <- IdClusters(distance.matrix, method = "ML", showPlot = TRUE, type = "dendrogram", myXStringSet = consensus, processors = NULL, verbose = TRUE)
#The selected model was TN93+G4
dendros <- IdClusters(distance.matrix, method = "ML", showPlot = TRUE, type = "dendrogram", myXStringSet = consensus, processors = NULL, verbose = TRUE)
dbDisconnect(dbConn)

##Editing cluster for suitable representation
clust <- as.dendrogram(as.hclust(dendrogram)) %>% dendextend::set("branches_lwd", 0.3) %>% dendextend::ladderize(right = TRUE) 
plot(clust)
failo <- as.dendrogram(as.hclust(dendros)) %>% dendextend::set("branches_lwd", 0.3) %>% dendextend::ladderize(right = TRUE) 
plot(failo)
sclustS <- as.dendrogram(sclust) %>% dendextend::set("branches_lwd", 0.3) %>% dendextend::ladderize(right = FALSE)
plot(sclustS)

##Heatmap
colnames(simp.hell) <- nms
maligned <-readDNAMultipleAlignment(multialign, format = "fasta")
mstring <- as(maligned, "DNAStringSet")
BrowseSeqs(mstring, htmlFile = "./Camilo/Otros/LOCURILLA/ppe_multiple_alignment.html")
#Check if ID from fasta and abundance table match
smatrix <- as.matrix(maligned) %>% rownames() %>% as.vector()
print(all(colnames(simp.hell)%in%gtools::mixedsort(smatrix, decreasing = FALSE)))
#If TRUE they match
heatmap <- heatmaply(normalize(simp.hell), 
                     Colv = clust,
                     Rowv = sclustS, 
                     main = "HCS PPE ASVs with sequences clustered by maximum likelihood plus SIMPROF for sites",
                     color = viridis(n = 256, 
                                     alpha = 1, 
                                     begin = 0, 
                                     end = 1, 
                                     option = "rocket")
)
saveWidget(heatmap, file = "./Camilo/Otros/LOCURILLA/ppe_heatmap_ML_SIMPROF.html")

##Heatmap without SIMPROF
wos.heatmap <- heatmaply(normalize(simp.hell), 
                         Colv = clust, 
                         Rowv = NA, 
                         main = "HCS PPE ASVs clustered by maximum likelihood", 
                         color = viridis(n = 256, 
                                         alpha = 1, 
                                         begin = 0, 
                                         end = 1, 
                                         option = "rocket")
)
saveWidget(wos.heatmap, file = "./Camilo/Otros/LOCURILLA/ppe_heatmap.html")


##Phyloseq (Daniel Vaulot)
rownames(taxa) <- taxa$ASV
rownames(environ) <- environ$Samples

taxa <- as.matrix(taxa)
simp <- as.matrix(simp)

si2p <- add_column(simp, Total = rowSums(simp[c(1:63)])) %>% as.matrix()
OTU <- otu_table(si2p, taxa_are_rows = TRUE)
TAXA <- tax_table(taxa)
samples <- sample_data(environ)
filum <- as.phylo(failo)
hcs <- phyloseq(OTU, TAXA, samples, filum)
hcs

plot_richness(hcs, color = "Date", x = "Site")
chao1 <- vegan::estimateR(as.data.frame(t(simp))) %>% 
  t() %>% 
  as.data.frame() %>% 
  add_column(Site = c(rep("Chanaral", 12), 
                      rep("Flamenco", 12), 
                      rep("Huasco", 12), 
                      rep("Pta.Choros", 5), 
                      rep("Quintero", 10), 
                      rep("LasCruces", 12)
                      ), 
             .name_repair = "minimal")

chao1$Site <- factor(chao1$Site, levels = unique(chao1$Site))

chao1 %>% 
  rename("Chao1" = "S.chao1") %>% 
  ggplot(aes(x = Site, y = Chao1, fill = Site)) + 
  geom_violin(width = 1.2) + 
  geom_boxplot(width = 0.2, color = "grey", alpha = 0.2) + 
  scale_fill_viridis(discrete = TRUE) + 
  geom_jitter(color = "black", size = 0.5, alpha = 0.5) +
  theme_ipsum() + 
  theme(legend.position = "none") + 
  ggtitle("Chao1 richness")

chao1 %>% 
  rename("ACE" = "S.ACE") %>% 
  ggplot(aes(x = Site, y = ACE, fill = Site)) + 
  geom_violin(width = 1.1) + 
  geom_boxplot(width = 0.2, color = "grey", alpha = 0.3) + 
  scale_fill_viridis(discrete = TRUE) + 
  geom_jitter(color = "black", size = 0.5, alpha = 0.5) + 
  theme_ipsum() + 
  theme(legend.position = "none") + 
  ggtitle("ACE richness")

shannon <- vegan::diversity(as.data.frame(t(simp)), index = "shannon", equalize.groups = FALSE) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  add_column(Site = chao1$Site, .name_repair = "minimal") %>% 
  dplyr::rename("Shannon" = V1)

shannon %>% 
  ggplot(aes(x = Site, y = Shannon, fill = Site)) + 
  geom_violin(width = 1.1) + 
  geom_boxplot(width = 0.1, color = "white", alpha = 0.3) + 
  scale_fill_viridis(discrete = TRUE) + 
  theme_ipsum() + 
  theme(legend.position = "none") + 
  ggtitle("Shannon diversity")

simpson <- vegan::diversity(as.data.frame(t(simp)), index = "simpson", equalize.groups = FALSE) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  add_column(Site = chao1$Site, .name_repair = "minimal", .before = "V1") %>% 
  dplyr::rename("Simpson" = V1)

simpson %>% 
  ggplot(aes(x = Site, y = Simpson, fill = Site)) + 
  geom_violin(width = 1.2) + 
  geom_boxplot(width = 0.1, color = "white", alpha = 0.3) + 
  scale_fill_viridis(discrete = TRUE) + 
  theme_ipsum() + 
  theme(legend.position = "none") + 
  ggtitle("Simpson diversity")

invsimpson <- vegan::diversity(as.data.frame(t(simp)), index = "invsimpson", equalize.groups = FALSE) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  add_column(Site = chao1$Site, .name_repair = "minimal", .before = "V1") %>% 
  dplyr::rename("InvSimpson" = V1)

invsimpson %>% 
  ggplot(aes(x = Site, y = InvSimpson, fill = Site)) + 
  geom_violin(width = 1.3) + 
  geom_boxplot(width = 0.1, color = "white", alpha = 0.3) + 
  scale_fill_viridis(discrete = TRUE) + 
  geom_jitter(color = "black", size = 0.5, alpha = 0.3) + 
  theme_ipsum() + 
  theme(legend.position = "none") + 
  ggtitle("Inverted Simpson diversity")

rde <- estimate_richness(hcs, split = TRUE, measures = NULL) %>% 
  add_column(Site = chao1$Site, .name_repair = "minimal", .before = "Observed")

print(all(rde$Fisher%in%fisher.alpha(as.data.frame(t(simp)))))

rde %>% 
  ggplot(aes(x = Site, y = Fisher, fill = Site)) + 
  geom_violin(width = 1.1) + 
  geom_boxplot(width = 0.2, color = "white", alpha = 0.5) + 
  scale_fill_viridis(discrete = TRUE, option = "D") + 
  theme_ipsum() + 
  theme(legend.position = "none") + 
  ggtitle("Fisher diversity")


##UniFrac
library(doParallel)
ufc <- UniFrac(hcs, weighted = TRUE, normalized = TRUE, parallel = TRUE)
ufcw <- UniFrac(hcs, weighted = FALSE, normalized = FALSE, parallel = TRUE)
ufc.m <- as.data.frame(as.matrix(ufc)) %>% 
  as_tibble(rownames = "A") %>% 
  pivot_longer(-A, names_to = "B", values_to = "UniFrac")

simpr <- with(smpr, simper(simp.hell, Site, permutations = 999, trace = TRUE, ordered = TRUE))
summary(simpr)[1]
str(simpr)

sbriff <- vegdist(simp.hell, method = "bray", diag = TRUE) %>% as.matrix() %>% 
  as_tibble(rownames = "A") %>% 
  pivot_longer(-A, names_to = "B", values_to = "bray")

sjriff <- vegdist(simp.hell, method = "jaccard", diag = TRUE) %>% as.matrix() %>% 
  as_tibble(rownames = "A") %>% 
  pivot_longer(-A, names_to = "B", values_to = "jaccard")

combined <- inner_join(sbriff, sjriff, by = c("A", "B")) %>% 
  inner_join(., ufc.m, by = c("A", "B"))

combined %>% 
  filter(A < B) %>% 
  ggplot(aes(x = bray, y = UniFrac)) + 
  geom_point() +  
  geom_smooth(se = FALSE)

combined %>% 
  filter(A < B) %>% 
  ggplot(aes(x = bray, y = jaccard)) + 
  geom_point() + 
  geom_smooth(se = FALSE)

##MANTEL TEST (comparar matrices de distancia)
set.seed(1000)
brunif <- mantel(sbriff, ufc, method = "spearman")

set.seed(1000)
jaccunif <- mantel(sjriff, ufc, method = "spearman")

set.seed(1000)
brajacc <- mantel(sbriff, sjriff, method = "spearman")


##metacoder
axat <- parse_tax_data(tax_data = as_tibble(taxa), 
                       datasets = list(counts = as_tibble(simp, rownames = "ASV")), 
                       mappings = c("ASV" = "ASV"), 
                       class_cols = 3:10)
axat$data$type_abund <- calc_taxon_abund(axat, "counts", cols = environ$Samples, groups = environ$Site)

hell <- decostand(as.data.frame(t(simp)), method = "hellinger") %>% t() %>% as.data.frame()
axoat <- parse_tax_data(tax_data = as_tibble(taxa), 
                       datasets = list(Abundance = as_tibble(hell, rownames = "ASV")), 
                       mappings = c("ASV" = "ASV"), 
                       class_cols = 3:10)
axoat$data$type_abund <- calc_taxon_abund(axoat, "Abundance", cols = environ$Samples, groups = environ$Site)

set.seed(420)
axat %>% filter_taxa(Chanaral >= 1) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Chanaral, 
            node_color = Chanaral, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in Chañaral")


set.seed(420)
axoat %>% filter_taxa(Chanaral >= 1) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Chanaral, 
            node_color = Chanaral, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in Chañaral")

set.seed(420) 
axat %>% filter_taxa(Flamenco >= 1) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Flamenco, 
            node_color = Flamenco, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in Flamenco")

set.seed(420) 
axoat %>% filter_taxa(Flamenco >= 1) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Flamenco, 
            node_color = Flamenco, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in Flamenco")

set.seed(420) 
axoat %>% filter_taxa(Huasco >= 1) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Huasco, 
            node_color = Huasco, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in Huasco")

set.seed(420) 
axoat %>% filter_taxa(Pta.Choros >= 1) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Pta.Choros, 
            node_color = Pta.Choros, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in Pta.Choros")

set.seed(420) 
axoat %>% filter_taxa(Quintero >= 1) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Quintero, 
            node_color = Quintero, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold",  
            title = "Taxa in Quintero")

set.seed(420) 
axoat %>% filter_taxa(Las.Cruces >= 1) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Las.Cruces, 
            node_color = Las.Cruces, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in Las Cruces")

set.seed(2) 
axat %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Total, 
            layout = "kamada-kawai", 
            node_color = Total, 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in HCS")

identification <- add_column(simp, Total = rowSums(simp[1:63]))
write.csv2(identification, "~/Desktop/LOCURILLA/rare_taxa.csv")

plot_richness(hcs, color = "Date", x = "Samples", measures = "Chao1")
plot_richness(hcs, color = "Date", x = "Samples", measures = "ACE")
plot_richness(hcs, color = "Date", x = "Samples", measures = "Shannon")
plot_richness(hcs, color = "Date", x = "Samples", measures = "Simpson")
plot_richness(hcs, color = "Date", x = "Samples", measures = "InvSimpson")
plot_richness(hcs, color = "Date", x = "Samples", measures = "Fisher")

rltv.Otu.Table <- function(x) {
  x.Data.rltv <- NULL
  for (i in 1:dim(x)[1]) {
    x.Data.rltv <- rbind(x.Data.rltv, x[i,]/apply(x, 1, function(x) sum(x))[i])
  }
  rownames(x.Data.rltv) <- rownames(x)
  invisible(x.Data.rltv)
}

rltv <- rltv.Otu.Table(as.data.frame(t(simp))) %>% as.data.frame()
axoltl <- parse_tax_data(tax_data = as_tibble(taxa), 
                         datasets = list(Abundance = as_tibble(as.data.frame(t(rltv)), rownames = "ASV")), 
                         mappings = c("ASV" = "ASV"), 
                         class_cols = 3:10)
axoltl$data$type_abund <- calc_taxon_abund(axoltl, "Abundance", cols = environ$Samples, groups = environ$Site)

set.seed(420)
axoltl %>% filter_taxa(Chanaral >= 0) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Chanaral, 
            node_color = Chanaral, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in Chañaral")

set.seed(420) 
axoltl %>% filter_taxa(Flamenco >= 0) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Flamenco, 
            node_color = Flamenco, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in Flamenco")

set.seed(420) 
axoltl %>% filter_taxa(Huasco >= 0) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Huasco, 
            node_color = Huasco, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in Huasco")

set.seed(420) 
axoltl %>% filter_taxa(Pta.Choros >= 0) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Pta.Choros, 
            node_color = Pta.Choros, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in Pta.Choros")

set.seed(420) 
axoltl %>% filter_taxa(Quintero >= 0) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Quintero, 
            node_color = Quintero, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold",  
            title = "Taxa in Quintero")

set.seed(420) 
axoltl %>% filter_taxa(Las.Cruces >= 0) %>% 
  heat_tree(node_label = taxon_names, 
            node_size = Las.Cruces, 
            node_color = Las.Cruces, 
            layout = "kamada-kawai", 
            initial_layout = "fruchterman-reingold", 
            title = "Taxa in Las Cruces")

