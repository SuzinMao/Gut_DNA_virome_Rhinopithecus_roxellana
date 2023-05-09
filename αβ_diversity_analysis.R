###Import TPM value file of all viral contigs in all samples
TPM_all <- readxl::read_xlsx("43130_contigs_TPM.xlsx", col_names = TRUE)

#The data format is changed to dataframe Type
TPM_all_update <- data.frame(TPM_all)

###α diversity analysis
##Install the package
library(vegan)
library(picante)
library(ggpubr)
library(RColorBrewer)
library(dplyr)

##Use the first column as the row name
rownames(TPM_all_update) <- TPM_all_update[,1]
TPM_all_update=TPM_all_update[,-1] 

##Calculate Shannon,Simpson,Richness index
df <- TPM_all_update
Shannon <- diversity(df, index = "shannon", MARGIN = 2, base = exp(1))
Simpson <- diversity(df, index = "simpson", MARGIN = 2, base = exp(1))
Richness <- specnumber(df, MARGIN = 2) #spe.rich ==sobs

##Put the results of the three indices together
index <- as.data.frame(cbind(Shannon,Simpson,Richness))

##Calculate obs,chao,ace index
tdf <- t(df)
tdf <- ceiling(as.data.frame(tdf))
obs_chao_ace <- t(estimateR(tdf))
obs_chao_ace <- obs_chao_ace[row.names(index),]
index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$Sobs <- obs_chao_ace[,1]

##Variance statistics and plotting of diversity index
#Add the sample name to the last column
index$samples <- row.names(index)
#Add grouping information for samples and fix the order of samples in abscissa
index$group <- as.factor(c(rep("Wild",15),rep("Captive",10)))
index$group1 <- factor(index$group,levels = c("Wild","Captive"))

##Assign file 'index' to a new file 'df2'
df2 <- index

##Detect whether the data conforms to a normal distribution
shapiro.test(df2$Shannon)  

##Use the Bartlett test for homogeneity of variance
bartlett.test(Shannon~group, data = df2)

##Drawing
mycol <- c("#316879","#ff9a8d")
#Shannon
Shannon.plot <- ggboxplot(df2, x = "group1", y = "Shannon", color = "group1", palette = mycol,
                          legend = "none",
                          add = "jitter",
                          outlier.shape = NA,ylim = c(0.0,8.0))+
  stat_compare_means(comparisons = list(c('Wild','Captive')),
                     method = 't.test',label = "p.signif",hide.ns = TRUE)+
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title.y = element_text(vjust = 0.2, size = 12, face = "bold"))

compare_means(Shannon ~ group1, data = df2, method = 't.test')

#Simpson
Simpson.plot <- ggboxplot(df2, x = "group1", y = "Simpson", color = "group1", palette = mycol,
                          legend = "none",
                          add = "jitter",
                          outlier.shape = NA,ylim = c(0.8,1.2))+
  stat_compare_means(comparisons = list(c('Wild','Captive')),
                     method = 'wilcox.test',label = "p.signif",hide.ns = TRUE)+
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title.y = element_text(vjust = 0.2, size = 12, face = "bold"))

compare_means(Simpson ~ group1, data = df2, method = 'wilcox.test')

#Ace
Ace.plot <- ggboxplot(df2, x = "group1", y = "Ace", color = "group1", palette = mycol,
                          legend = "none",
                          add = "jitter",
                          outlier.shape = NA)+
  stat_compare_means(comparisons = list(c('Wild','Captive')),
                     method = 't.test',label = "p.signif",hide.ns = TRUE)+
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title.y = element_text(vjust = 0.2, size = 12, face = "bold"))

compare_means(Ace ~ group1, data = df2, method = 't.test')

##Save the above 3 figures
library("gridExtra")
library("cowplot")
plot_grid(Shannon.plot,Simpson.plot,Ace.plot, ncol = 3, nrow = 1)

ggsave("α_diversity.pdf",width = 10, height = 4)

###β diversity analysis:PCoA+PERMANOVA
data <- TPM_all_update
df <- t(data)

##Load the R package required to calculate distance
library(vegan)
#Calculate bray_curtis distance
df_dist <- vegdist(df, method="bray", binary=F)
df_pcoa <- cmdscale(df_dist, k=3, eig=T)
df_pcoa_points <- as.data.frame(df_pcoa$points)
sum_eig <- sum(df_pcoa$eig)
eig_percent <- round(df_pcoa$eig/sum_eig*100,1)

colnames(df_pcoa_points) <- paste0("PCoA", 1:3)

#Add group information
group_text <- c(rep("Wild",15),rep("Captive",10))
df_pcoa_points$group <- factor(group_text)

#Build the metadata information of 'df', which has only one column, that is, the grouping information
df.env <- data.frame(group_text)

df_pcoa_result <- df_pcoa_points
head(df_pcoa_result)

##Drawing
library(ggplot2)
ggplot(df_pcoa_result, aes(x=PCoA1, y=PCoA2, color=group)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  geom_point(size=4
  ) +
  theme_classic()
# + stat_ellipse(level=0.6) 

#PERMANOVA:test the significance of PCoA results
df.div <- adonis2(df ~ group_text, data = df.env, permutations = 999, method="bray")

df_adonis <- paste0("adonis R2: ",round(df.div$R2,2), "; P-value: ", df.div$`Pr(>F)`)
p <- ggplot(df_pcoa_result, aes(x=PCoA1, y=PCoA2, color=group, group = group)) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title=df_adonis) +
  geom_point(size=3) + 
  theme_classic()
#p <- p+stat_ellipse(level = 0.95) 

devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

# This is a wrapper function for multilevel pairwise comparison 
# using adonis() from package 'vegan'. 
# The function returns adjusted p-values using p.adjust(); Benjamini-Hochberg correction was used here
df.pairwise.adonis <- pairwise.adonis(x=df, factors=df.env$group_text, sim.function = "vegdist",
                                      sim.method = "bray",
                                      p.adjust.m = "BH",
                                      reduce = NULL,
                                      perm = 999)

##Stitch together the plots of PCoA results and PERMANOVA results
library(ggpubr)
library(patchwork)
tab2 <- ggtexttable(df.pairwise.adonis[,c("pairs","R2","p.value","p.adjusted")], rows = NULL, 
                    theme = ttheme("blank")) %>% 
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)  %>% 
  tab_add_hline(at.row = nrow(df.pairwise.adonis)+1, row.side = "bottom", linewidth = 1)  

p2 = p + tab2 
p2


p2 + plot_layout(design=c(area(1,1), area(2,1)))
