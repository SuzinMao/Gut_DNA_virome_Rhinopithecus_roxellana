###Comparison of relative abundance of each virus family in different groups

##At first, calculate TPM values for each viral family in each sample
#Import the taxonomy information file for the virus
DNA_family <- readxl::read_xlsx("family_tax_info.xlsx")

#Import the TPM values of all viruses in all samples
TPM_all_update <- readxl::read_xlsx("TPM_all_viruses.xlsx", col_names = TRUE)

#Change the file format to 'dataframe'
TPM_all_update <- data.frame(TPM_all_update)

#Merge file 'DNA_family' and file 'TPM_all_update';
#Put the TPM value and taxonomy information of the viruses in one file
TPM_famliy <- dplyr::left_join(DNA_family,TPM_all_update,by = "contig")

#Remove the first column of 'TPM_famliy' and use it to calculate species composition information
TPM_family_rm1 <- TPM_famliy[,2:27]
TPM_family_sorted <- TPM_family_rm1[order(TPM_family_rm1[,1]),]
TPM_family_sum <- aggregate(TPM_family_sorted[-1], TPM_family_sorted["family"], sum)

#Obtain TPM values for each viral family in each sample
write.table(TPM_family_sum,"TPM_family_sum.txt", sep = "\t")



##Then,relative abundance of viral family in each group
##Histogram
#Divide all samples into wild and captive groups and display the composition of viral family in each group
grouped <- data.frame(wild=rowSums(TPM_family_sum[,c(1:15)]),
                      captive=rowSums(TPM_family_sum[c(16:25)]))

#Change the base from the 6th power of 10 to 1
da <- apply(grouped,2,function(x) x/sum(x))
e <- sum(da[,1:2]) 

#Histogram
library(reshape2)
library(ggplot2)
df <- melt(da)
head(df)
mycol <- c(119,132,147,454,89,404,123,529,463,104,
           552,28,54,84,256,100,558,43,652,31,610,
           477,588,99,81,503,562,76,96,495)
mycol <- colors()[rep(mycol,20)]

p <- ggplot(df, aes(x=Var2, y=value)) +
  geom_bar(aes(fill=Var1), stat = "identity", width = 0.5)+
  scale_fill_manual(values = mycol)+
  guides(fill=guide_legend(reverse = T,title = NULL,ncol = 4))+
  labs(x="Sample",y="Relative abundance") +
  theme_classic()+
  theme(legend.position = "bottom") +
  theme(axis.title.x = element_text(vjust = 2, size = 10, face = "bold")) +  
  theme(axis.title.y = element_text(vjust = 2, size = 10, face = "bold")) +
  theme(axis.text.x=element_text(vjust=1,size=10,face = "bold")) +  
  theme(axis.text.y=element_text(vjust=0.5,size=10,face = "bold")) 

p



##At last, calculate whether there is a significant difference in the relative abundance of each viral family between different groups
##Boxplots
da <- apply(TPM_family_sum,2,function(x) x/sum(x))Transpose data
#Transpose data
da <- as.data.frame(t(da))

#Add group information
da$group <- as.factor(c(rep("Wild",15),rep("Captive",10)))
#Fixed horizontal coordinate order
da$group1 <- factor(da$group,levels = c("Wild","Captive"))

library(ggpubr)   
library(ggplot2)
library(vegan)

#Microviridae
shapiro.test(da$Microviridae)
mycol <- c("#316879","#ff9a8d")

Microviridae.plot <- ggboxplot(da, x = "group1", y = "Microviridae", color = "group1", palette = mycol,
                               legend = "none",
                               add = "jitter",
                               outlier.shape = NA)+
  geom_signif(                       
    comparisons=list(c("Wild","Captive")), 
    step_increase = 0.1,
    test="wilcox.test",                     
    map_signif_level=F                 
  )+
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title.y = element_text(vjust = 0.2, size = 12, face = "bold"))

compare_means(Microviridae ~ group1, data = da, method = 'wilcox.test')

#Siphoviridae
shapiro.test(da$Siphoviridae) 
bartlett.test(Siphoviridae~group1, data = da) 
Siphoviridae.plot <- ggboxplot(da, x = "group1", y = "Siphoviridae", color = "group1", palette = mycol,
                               legend = "none",
                               add = "jitter",
                               outlier.shape = NA)+
  geom_signif(                        
    comparisons=list(c("Wild","Captive")), 
    step_increase = 0.1,
    test="wilcox.test",                     
    map_signif_level=F                 
  )+
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title.y = element_text(vjust = 0.2, size = 12, face = "bold"))

compare_means(Siphoviridae ~ group1, data = da, method = 'wilcox.test')

#Caudovirales
#Merge values of several viral families belonging to the order Caudovirales
da$Caudovirales <- rowSums(da[,3:5]) 
shapiro.test(da$Caudovirales) 
bartlett.test(Caudovirales~group1, data = da)  
Caudovirales.plot <- ggboxplot(da, x = "group1", y = "Caudovirales", color = "group1", palette = mycol,
                              legend = "none",
                              add = "jitter",
                              outlier.shape = NA)+
  geom_signif(                         
    comparisons=list(c("Wild","Captive")), 
    step_increase = 0.1,
    test="wilcox.test",                     
    map_signif_level=F                 
  )+
  theme(panel.grid = element_line(color = NA),
        panel.grid.minor = element_line(color = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        axis.text.x = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = "black", face = "bold", hjust = 0.5),
        axis.title.y = element_text(vjust = 0.2, size = 12, face = "bold"))

compare_means(Caudovirales ~ group1, data = da, method = 'wilcox.test')

#Save plots
library("gridExtra")
library("cowplot")
plot_grid(Microviridae.plot,Siphoviridae.plot,Caudovirales.plot, ncol = 3, nrow = 1)
ggsave("3_family_boxplots.pdf",width = 15, height = 4)

