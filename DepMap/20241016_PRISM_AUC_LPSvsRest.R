library(data.table)

DrugScreen <- fread("/data/cephfs-1/home/users/rauertc_c/work/DepMap/Data/20241016_Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_subsetted.csv")
View(DrugScreen)

#LPS067=ACH-001807
#LPS141=ACH001799
#93T449=ACH001794
#94T778=ACH001795
#95T1000=ACH001796
#SW872=ACH-002310
#LPS510=ACH-001804
#LPS27=ACH-001793
#LPS853=ACH-001802
#LPS6=ACH-001791
#KMLS1=ACH-001540

LPS_Gr <- c("ACH-001799","ACH-001794","ACH-001795","ACH-001796","ACH-001807","ACH-002310","ACH-001804","ACH-001793","ACH-001802","ACH-001791","ACH-001540")
Out_Gr <- setdiff(DrugScreen$V1,LPS_Gr)
#> length(DrugScreen$V1)
#[1] 480
#> length(Out_Gr)
#[1] 480

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
pvals_int <- pvals[p_fdr<0.05]
pvals_int
#Empty data.table (0 rows and 3 cols): compound,p_val,p_fdr