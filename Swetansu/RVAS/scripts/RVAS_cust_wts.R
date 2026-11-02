# Helper functions for burden testing

library(SKAT)
library(ggplot2)
library(pbapply)

# Generate SSD_INFO file
gen_ssd_info = function(dir_path, plink_fil_nam, ssdfil_name){
  File.Bed <- paste0(dir_path, "plink_fil_nam", ".bed")
  File.Bim <- paste0(dir_path, "plink_fil_nam", ".bim")
  File.Fam <- paste0(dir_path, "plink_fil_nam", ".fam")
  File.SetID <- paste0(dir_path, "ssdfil_name", ".SetID") #modified SetID or new file
  File.SSD <- paste0(dir_path, "ssdfil_name", ".SSD")
  File.Info <- paste0(dir_path, "ssdfil_name", ".SSD.info")
  
  Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
  
  SSD.INFO <- Open_SSD(File.SSD, File.Info)
  return(SSD.INFO)
}

##################################
#Flag cases with zero variants; retain cases with non-zero variants
# report the number variants in cases
# input files: plink FAM, SSD.Info
##################################

case_nonzero = function(Set_Index, SSD_INFO, fam_file){
  Z = Get_Genotypes_SSD(SSD_INFO, Set_Index, is_ID = TRUE)
  #add hits column to fam_file, yes for GT == 1 and no for GT == 0
  fam_file$hits = ifelse(fam_file$FID %in% fam_file$FID[which(rowSums(Z!= 0) == 1)], "yes", "no")
  #return gene or set index only if cases have non-zero variants; report the number of variants in cases
  if(sum(fam_file[fam_file$hits %in% "yes"]$Phenotype) > 0) {
    case_num = sum(fam_file[fam_file$hits %in% "yes"]$Phenotype)
    return(cbind.data.frame("set_ind" = Set_Index, "case_MAC" = case_num)) }
  else {return(NULL)}
}

# run example
# test_genes = lapply(seq(1, SSD.INFO$nSets, 1), function(x) case_nonzero(x, SSD_INFO, fam_file))

##################################
# Modify and update SSD.Info file
# input files: SSD.SetID, SSD.Info
# input args: set size
# add custom args to modify the function
# save as tab-separated output
##################################

update_ssdinfo = function(dir_path, ssdfil_name, SSD_INFO, case_var_fil, snp_set_size){
  SSD_INFO$SetInfo = SSD_INFO$SetInfo[SSD_INFO$SetInfo$SetIndex %in% case_var_fil$Set_Index,]
  SSD_INFO$SetInfo = SSD_INFO$SetInfo[SSD_INFO$SetInfo$SetSize > snp_set_size,]
  SSD_INFO$nSets = dim(SSD_INFO$SetInfo)[1]
  SSD_INFO$nSNPs = sum(SSD_INFO$SetInfo$SetSize)
  SSD_SetID = read.delim(paste0(dir_path, ssdfil_name, ".SetID"), sep = "\t", stringsAsFactors = F, header = F)
  SSD_SetID_mod = SSD_SetID[SSD_SetID$V1 %in% SSD_INFO$SetInfo$SetID,]
  write.table(SSD_SetID_mod, paste0(dir_path, ssdfil_name, "_mod.SetID"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  Close_SSD()
}

get_variant_CADD = function(dir_path, anno_file_setid, ssdfil_name_mod){
  tab_wt = read.delim(paste0(dir_path, anno_file_setid), sep = "\t", header = T, stringsAsFactors = F)
  
  CADD_wt_var = tab_wt[,c(1:3,10)]
  CADD_wt_var$CADD_PHRED = as.numeric(CADD_wt_var$CADD_PHRED)
  
  CADD_wt_var$CADD_PHRED[is.na(CADD_wt_var$CADD_PHRED)] <- 0
  #Assign Maximum score to PTV not scored by CADD
  CADD_wt_var$CADD_PHRED = ifelse(grepl("frameshift", CADD_wt_var$CADD_PHRED, ignore.case = T),
                                  max(CADD_wt_var$CADD_PHRED), CADD_wt_var$CADD_PHRED)
  #Harmonise CADD scores with updated ssd_info_mod file
  set_id_mod = read.delim(paste0(dir_path, ssdfil_name_mod, ".SetID"), sep = "\t", header = T, stringsAsFactors = F)
  set_id_mod$V3 = CADD_wt_var[match(set_id_mod$V2, CADD_wt_var$variation), 4]
  write.table(set_id_mod[,2:3], paste0(dir_path, "SNP_wt_CADD.txt"),
              sep = "\t", col.names = F, row.names = F, quote = F)
}


#Close_SSD() #after analysis is complete; close SSD.INFO.

##################################
#Input file preparation for SKAT
#################################

#Step.1 Generate SSD.INFO file
SSD_INFO = gen_ssd_info(dir_path = dir_path, plink_fil_nam = plink_fil_nam, ssdfil_name = ssdfil_name)

#Step.2 Generate FAM file with Covariates
# Covariate file (tab separated) contains ancestry PCs and other nuisance variables
# Ensure the cases are set to 1 and controls to 0 in FAM$Phenotype
File.Fam <- paste0(dir_path, "plink_fil_nam", ".fam")
File.Cov = paste0(dir_path, "cov_fil_nam", ".cov")
FAM = Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary = TRUE, cov_header = FALSE)

#Step.3 Apply filters to update SSD_INFO file
#3a. Retain genes with non-zero variants in cases
#3b. Retain genes with SNP set size greater than 1

burden_genes = lapply(seq(1, SSD.INFO$nSets, 1), function(x) case_nonzero(x, SSD_INFO, fam_file))
burden_genes_df = data.frame(do.call("rbind", burden_genes))

update_ssdinfo(dir_path = dir_path, ssdfil_name = ssdfil_name, SSD_INFO = SSD_INFO, 
               case_var_fil = burden_genes_df, snp_set_size = 1)

#Step.4 Generate the updated SSD_INFO file
ssdfil_name_mod = paste0(dir_path, ssdfil_name, "_mod.SetID")
SSD.INFO_mod = gen_ssd_info(dir_path = dir_path, plink_fil_nam = plink_fil_nam, ssdfil_name = ssdfil_name_mod)

#Step.5 Tests
obj = SKAT_Null_Model(Phenotype ~ Sex + Cov1 + Cov2, out_type = "D")

out.skat_burden_rob = SKATBinary_Robust.SSD.All(SSD.INFO_mod, obj, method = "Burden")

#Step.5a Tests with custom SNP weights

get_variant_CADD(dir_path = dir_path, anno_file_setid = anno_file_setid, ssdfil_name_mod = ssdfil_name_mod)

obj.SNP_wt_CADD = Read_SNP_WeightFile(paste0(dir_path,"SNP_wt_CADD.txt"))
cadd.skat_burden_rob = SKATBinary_Robust.SSD.All(SSD.INFO_mod, obj, obj.SNPWeight = obj.SNP_wt_CADD, method = "Burden")

#Step.6 Merge results with case_MAC
burden_filt = lapply(seq(1, SSD.INFO_mod$nSets, 1), function(x) case_nonzero(x, SSD.INFO_mod, fam_file))
burden_filt_df = data.frame(do.call("rbind", burden_filt))

df_res = cbind.data.frame(out.skat_burden_rob$results, "case_MAC" = burden_filt_df)
df_res_fin = df_res[order(df_res$P.value, decreasing = F),]

#Step.6a Merge CADD results with case_MAC
df_res_cadd = cbind.data.frame(cadd.skat_burden_rob$results, "case_MAC" = burden_filt_df)
df_res_cadd_fin = df_res_cadd[order(df_res_cadd$P.value, decreasing = F),]

