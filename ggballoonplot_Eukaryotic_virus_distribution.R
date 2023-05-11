##Perform balloon plotting of the distribution of all eukaryotic viruses detected in the all samples
#Install packages
library(ggplot2)
library(ggpubr)

#Import the file which contained TPM values of all eukaryotic viruses detected in the all samples
inputFile <- readxl::read_xlsx("TPM_values.xlsx")

#Use the first column as the row name
rownames(inputFile) <- inputFile[,1]
inputFile=inputFile[,-1]

#Draw a balloon chart
rt <- ggballoonplot(inputFile, fill = 'value') + gradient_fill(c('blue','white','red'))

#Export the chart
ggsave("balloonplot.pdf",width = 10, height = 4)


