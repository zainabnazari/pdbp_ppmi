
rm(list = ls())
options("warn"=0)  #print max 10 warnings on screen
library(limma)      # for package 'limma'
library("preprocessCore")


##############################################  Parameters to set


				
work_path ="D:/disco_H/bandi_grants/Regione_Lazio_FILAS_2016_Confessore/E_LIFE_submitted_13dic2016/tesi_laurea/Zainab_Nazari/PPMI_analysis/R_limma"



name_file_input_PPMI="ir3_rna_step1_no_normalization_PzID_filtered_ver_02giu2024.txt"
name_file_input_PDBP="PDBP_bl_matrix_no_parents_with_PD_T.txt"

name_batch_factor_file_input="factor_PPMI_PDBP_counts_shared_genes_02giu2024.txt"


##############################################
   
setwd(work_path)


mydata_PPMI=read.table(file=name_file_input_PPMI, sep = "\t", quote = "\"",row.names=1,header=TRUE, fill=TRUE)  # re-read data into a dataframe with just numbers as real data

mydata_PDBP=read.table(file=name_file_input_PDBP, sep = "\t", quote = "\"",row.names=1,header=TRUE, fill=TRUE)  # re-read data into a dataframe with just numbers as real data

# rownames_PDBP=read.table(file="rownames_PDBP.txt", sep = "\t", quote = "\"",header=F, fill=F)  # re-read data into a dataframe with just numbers as real data

# rownames(mydata_PDBP)=unlist(rownames_PDBP)

# write.table(mydata_PDBP,file=name_file_input_PDBP, sep = "\t", row.names=TRUE,col.names=TRUE )  # re-read data into a dataframe with just numbers as real data

# write.table(rownames(mydata_PDBP), "rownames_PDBP.txt", sep=".",row.names=F, col.names=F)


myfactors=read.table(file=name_batch_factor_file_input, sep = "\t", quote = "\"",header=TRUE, fill=TRUE)  # re-read data into a dataframe 


shared_genes=intersect(rownames(mydata_PPMI), rownames(mydata_PDBP))


shared_PPMI_genes_index = match(shared_genes, rownames(mydata_PPMI))
shared_mydata_PPMI=mydata_PPMI[shared_PPMI_genes_index,]
shared_mydata_PPMI=shared_mydata_PPMI[order(rownames(shared_mydata_PPMI)),]


shared_PDBP_genes_index = match(shared_genes, rownames(mydata_PDBP))
shared_mydata_PDBP=mydata_PDBP[shared_PDBP_genes_index,]
shared_mydata_PDBP=shared_mydata_PDBP[order(rownames(shared_mydata_PDBP)),]


union_PPMI_PDBP_mydata=cbind(shared_mydata_PPMI,shared_mydata_PDBP)



write.table(union_PPMI_PDBP_mydata, "PPMI_PDBP_counts_shared_genes_02giu2024.txt", sep="\t",row.names=TRUE, col.names=TRUE)

fac = factor(myfactors$Diagnosis,levels=c("PD","CTR"))                  #  factor() discretizza i gruppi, creando dei fattori
sex_fac = factor(myfactors$Sex,levels=c("M","F"))




# Clinical_center_fac = factor(myfactors$Clinical_center)
# RIN_covariate=as.vector(myfactors$RIN)


dge<- DGEList(counts = union_PPMI_PDBP_mydata)  

		# disegnare/creare la matrice factor, with sex correction
design <- model.matrix(~0 + fac )       #  genera una matrice 'design' ,                               le righe corrispondono ai parametri da stimare, le colonne alla condizione sperimentale


# dge <- calcNormFactors(dge,method="TMM")


logCPM <- cpm(dge, log=TRUE, prior.count=2)  #The prior count is used here to avoid log(0). The logCPM values can then be used in any standard limma pipeline, using the trend=TRUE argument when running eBayes or treat. 
# logCPM_filtered=logCPM
plot(colSums(logCPM)))


logCPM_sex_batch_effect_removed = removeBatchEffect(logCPM, batch=sex_fac, design=design)

logCPM_sex_batch_effect_removed_QuantileNorm = normalize.quantiles(logCPM_sex_batch_effect_removed)

# plot(colSums(logCPM_sex_batch_effect_removed_QuantileNorm))
# boxplot(logCPM_sex_batch_effect_removed_QuantileNorm[,1:100])

colnames(logCPM_sex_batch_effect_removed_QuantileNorm)=colnames(logCPM)
rownames(logCPM_sex_batch_effect_removed_QuantileNorm)=rownames(logCPM)

write.table(logCPM_sex_batch_effect_removed_QuantileNorm, "PPMI_PDBP_shared_genes_Log2CPM_sex_effect_removed_QuantileNorm_02giu2024.txt", sep="\t",row.names=TRUE, col.names=TRUE)

write.table(logCPM_sex_batch_effect_removed_QuantileNorm[,1:545], "PPMI_subset_from_PPMI_PDBP_shared_genes_Log2CPM_sex_effect_removed_QuantileNorm_02giu2024.txt", sep="\t",row.names=TRUE, col.names=TRUE)

write.table(logCPM_sex_batch_effect_removed_QuantileNorm[,546:1752], "PDBP_subset_from_PPMI_PDBP_shared_genes_Log2CPM_sex_effect_removed_QuantileNorm_02giu2024.txt", sep="\t",row.names=TRUE, col.names=TRUE)














