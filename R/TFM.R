library(readr)
library(reshape2) #Para el uso de dcast()
library(ggplot2)
library(pheatmap)
library(DESeq2)
library(stringr) # Para usar scale_fill_discrete()
library(Heatplus)
library(gmodels)

resistomas <- as.data.frame(read_delim("C:/Users/dfern/TFM/archivos_preprocesados/resistomas",
                         "\t", escape_double = FALSE, trim_ws = TRUE))

metadata <- as.data.frame(read_csv("C:/Users/dfern/TFM/archivos_preprocesados/SraRunTable.txt"))

metadata$Isolation_source <- toupper(abbreviate(metadata$Isolation_source))

summary(resistomas$Sample)
table(resistomas[,3])
unique(resistomas$`ARO Accession`)
table(resistomas$Sample, resistomas$`ARO Accession`)

# Crear tabla de abundancias

tabla_ARGs <- as.data.frame(dcast(resistomas, `ARO Name`~Sample, value.var = "Relative abundance",sum)) # Se suman aquellos genes que se repiten en cada muestra
# Se quiere eliminar el título "sample" de la tabla
row.names(tabla_ARGs) <- tabla_ARGs$`ARO Name`
mapa_calor <- tabla_ARGs[,-1]
row.names(metadata) <- metadata$Run
pheat_row_ano <- resistomas[,c(11,14)]
pheat_row_ano <- pheat_row_ano[!duplicated(pheat_row_ano),]
row.names(pheat_row_ano) <- pheat_row_ano$`ARO Name`
pheat_row_ano <- pheat_row_ano[,-1, drop = FALSE]
names(pheat_row_ano)[1] <- "Familia génica AMR"
pheat_col_ano <- metadata[,-1, drop =FALSE]
names(pheat_col_ano)[19] <- "Tipo de muestra"
pheatmap(scale(as.matrix(mapa_calor)), 
         annotation_col = pheat_col_ano[,19, drop = FALSE], 
         annotation_row = pheat_row_ano, cellwidth = 30, 
         cellheight = 10, cluster_cols = TRUE, 
         cluster_rows = TRUE)


# Abundance chart
barras_abundancia <- as.data.frame(dcast(resistomas, Sample~`AMR Gene Family`, value.var = "Relative abundance",sum))
barras_abundancia <- melt(barras_abundancia, id.vars = "Sample", variable.name = "AMR Gene Family")
names(metadata)[1] <- "Sample"
barras_abundancia <- merge(barras_abundancia, metadata[,c(1,20)], by="Sample")
mean.df <- as.data.frame(aggregate(value ~ `AMR Gene Family` + Isolation_source, barras_abundancia, (function(x) c(mean = mean(x), sd = sd(x)))))
mean.df <- do.call(data.frame, mean.df)
names(mean.df)[3] <- "Mean"
names(mean.df)[4] <- "SD"
barras_abundancia_mean_sd <- merge(barras_abundancia, mean.df)
ggplot(barras_abundancia_mean_sd, aes(x = Isolation_source, y = Mean, fill = AMR.Gene.Family)) + 
        xlab("Tipo de muestra") + ylab("Abundacia relativa") + 
        geom_bar(stat = "identity") + guides(fill=guide_legend(title="Familia génica AMR")) + 
        geom_errorbar(aes(ymax = Mean + SD, ymin=Mean - SD)) +
        theme(legend.title = element_text(face='bold'), axis.title = element_text(face='bold'))


#ggplot(barras_abundancia, aes(x = Sample, y = value, fill = `AMR Gene Family`)) + 
        xlab("Tipo de muestra") + ylab("Abundacia relativa") + 
        geom_bar(stat = "identity") + guides(fill=guide_legend(title="Familia génica AMR")) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
              legend.title = element_text(face='bold'), axis.title = element_text(face='bold')) + 
        facet_grid(~Isolation_source, switch = "x", scales = "free_x", space = "free_x") + 
        theme(panel.spacing = unit(0, "lines"), 
              strip.background = element_blank(),
              strip.placement = "outside")

#ggplot(barras_abundancia, aes(x = Sample, y = value, fill = str_wrap(`AMR Gene Family`,20))) + xlab("Muetra") + ylab("Abundacia relativa") + geom_bar(stat = "identity") + scale_fill_discrete(label = function(x) stringr::str_trunc(x, 50)) + theme(legend.key.height=unit(1, "cm")) + guides(fill=guide_legend(title="Familia génica AMR")) + theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))


#Estadística
abundance_table <- dcast(resistomas, Sample~`ARO Accession`, value.var = "Relative abundance",sum)
test <- as.data.frame(metadata[,c(1,20)])
abundance_table$ARGsum <- rowSums(abundance_table[,2:47])

#¿Conservar las muestras sin ningún antibiótico?
test <- merge(test, abundance_table[,c(1,48)], by = "Sample", all = TRUE)
test[is.na(test)]<- 0

#Se aplica kurast-wallis o anova

kruskal.test(ARGsum~Isolation_source, data = test)


#Metadata para deseq2

sample <- as.factor(c("SRR13882924", "SRR13882925", "SRR13882926"))

treatment <- as.factor(c("treated", "untreated", "treated"))

df <- data.frame(sample, treatment)

#Preparamos los datos con los counts
deseq_table <- as.data.frame(dcast(resistomas, `ARO Accession`~Sample, value.var = "Counts",sum)) # Se suman aquellos genes que se repiten en cada muestra
deseq_table$SRR13882931 <- 0
deseq_table$SRR13882935 <- 0
head(deseq_table)
colnames(deseq_table)
deseq_table <- deseq_table[,c(1,2,3,4,5,6,7,8,10,9,11)]
View(deseq_table)
dds <- DESeqDataSetFromMatrix(countData=deseq_table, 
                              colData=metadata[-c(8,10),c(1,20)], 
                              design=~Isolation_source, tidy = TRUE)
dds
dds<-estimateSizeFactors(dds, type = "poscounts") #Evitar mensaje "every gene contains at least one zero, cannot compute log geometric means"
dds <- DESeq(dds)
res <- results(dds)
#Null hypothesis that there is no effect of the treatment on the gene
#and that the observed difference between treatment and control was merely caused by experimental
#variability
res
results(dds, tidy=TRUE)
summary(res)
vsdata <- varianceStabilizingTransformation(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Isolation_source") + labs(colour="Tipo de muestra") + theme(legend.title = element_text(face='bold'), axis.title = element_text(face='bold'))



