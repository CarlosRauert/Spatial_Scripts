library(data.table)

DrugScreen <- fread("/data/cephfs-1/home/users/rauertc_c/work/DepMap/Data/20241016_Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv")
View(DrugScreen)

#LPS141=ACH001799
#93T449=ACH001794
#94T778=ACH001795
#95T1000=ACH001796
LPS_Gr <- c("ACH-001799","ACH-001794","ACH-001795","ACH-001796")
Out_Gr <- setdiff(colnames(DrugScreen),c("V1",LPS_Gr))

#Remove Compound names
log_fold_changes <- DrugScreen[,-1,with=FALSE]

#Apply Wilcoxon test for each compound
i=1
results <- lapply(1:nrow(log_fold_changes), function(i) {
  # Perform Wilcoxon test on two groups for each compound
  row=log_fold_changes[i,]
  group1_data <- as.numeric(row[,..LPS_Gr])
  group2_data <- as.numeric(row[,..Out_Gr])
  
  # Remove NA values from both groups
  group1_data <- group1_data[!is.na(group1_data)]
  group2_data <- group2_data[!is.na(group2_data)]
  
  # Perform Wilcoxon test only if both groups have enough non-missing values
  if (length(group1_data) > 0 & length(group2_data) > 0) {
    test <- wilcox.test(group1_data, group2_data, exact = FALSE)
    return(test$p.value)  # Return p-value
  } else {
    return(NA)  # Return NA if there are not enough valid observations
  }

  test <- wilcox.test(group1_data, group2_data, exact = FALSE)
  
  # Return p-value from Wilcoxon test
  test$p.value
})

pvals <- data.table(compound=DrugScreen$V1, p_val=results)
pvals[,p_fdr:=p.adjust(p_val, method="fdr")]
pvals[,p_bon:=p.adjust(p_val, method="bonferroni")]
pvals_int <- pvals[p_val<0.05]
pvals_int
fwrite(results, file="20241016_Wilcoxon_pVals_PRISM_Primary.csv")