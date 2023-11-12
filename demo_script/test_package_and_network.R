#----------------- Set parameters -------------------
# load library
library(ChloroDBPFinder)
library(ISFrag)

# specify the path of models
binary_model_file <- "C:/Users/User/Tingting/2022-11-03-Cl_project/ChloroDBP Hunter/06-02/binary_with_noise_ntree_500_mtry_3.rds"
multi_model_file <- "C:/Users/User/Tingting/2022-11-03-Cl_project/ChloroDBP Hunter/06-02/multi_with_noise_ntree_500_mtry_5.rds"

# specify the path of raw lcms data
mzmldir <- "C:/Users/User/Desktop/testmzML"

# specify the format of the raw lcms data
lcmspattern <- ".mzML" # String: ".mzML" or ".mzXML"
mzMLfile <- list.files(pattern = lcmspattern, mzmldir)

# specify whether ues a customized feature table
# if true, change the path of the customized table
use_customized_table <- FALSE
customized_table <- 'C:/Users/User/Desktop/my_customized_feature_table.csv' 

# specify whether detect in-source fragments
isfrag <- TRUE # Boolean: TRUE or FALSE

# specify the path of the MS/MS spectra database, if users want to database search for compound annotation
Cl_db_path <- "C:/Users/User/Tingting/2022-11-03-Cl_project/ChloroDBP Hunter/06-02/Cl_compounds_in_NIST.csv"

# specify the path of the reference table
ref_path <- "C:/Users/User/Tingting/2022-11-03-Cl_project/ChloroDBP Hunter/06-02/reference_table.csv"

#---------------- Main program ---------------------------------
### --------- 1. extraction of Cl-containing features ----------

for(i in 1:length(mzMLfile)){
  
  # load customized table or extract chemical features
  if(use_customized_table){
    peaks <- read.csv(customized_table)
  }else{
    peaks <- extractPeak(mzMLdirectory = mzmldir, mzMLfile = mzMLfile[i] )
    write.csv(peaks, paste0(mzmldir,"/",strsplit(mzMLfile[i], split = lcmspattern)[[1]][1], "_",nrow(peaks), "_peaks_with_MS2.csv"),row.names = FALSE )
  }
  
  # determine chlorinated compounds
  binary_rf_model <- readRDS(binary_model_file)
  multi_rf_model <- readRDS(multi_model_file)
  xcmsrawlcms <- eicRawlcms(mzMLdirectory = mzmldir, mzMLfile = mzMLfile[i])
  cl_tb <- selectCl(mzMLdirectory = mzmldir, mzMLfile = mzMLfile[i], original_ft = peaks,
                    binary_model = binary_rf_model, multi_model = multi_rf_model)
  
  # identify in source fragment based on ISFrag package
  if(isfrag){
    customFT <- cl_tb
    customFT$Adduct <- 0
    customFT$isotope <- 0
    rownames(customFT) <- peaks$featureID
    if(grepl("mzXML", mzMLfile[i])){filename <- strsplit( mzMLfile[i], split=".mzXML")[[1]][1]}else{filename <- strsplit( mzMLfile[i], split=".mzML")[[1]][1]}
    ISFdirectory_name <- paste0(mzmldir,"/inSourceFrag_", filename)
    dir.create(ISFdirectory_name)
    file.copy(from = paste0(mzmldir, "/", mzMLfile[i]), to = ISFdirectory_name)
    featureTable <- ISFrag::ms2.assignment(MS2directory = ISFdirectory_name, customFT = customFT)
    featureTable <- featureTable[,-1]
    level3 <- ISFrag::find.level3(MS1directory = ISFdirectory_name,
                                  MS1.files = mzMLfile[i],
                                  featureTable = featureTable,
                                  type = "single")
    level2 <- ISFrag::find.level2(ISFtable = level3)
    level1 <- ISFrag::find.level1(ISF_putative = level2)
    results <- ISFrag::get.ISFrag.results(ISF_List = level1, featureTable = featureTable)
    result <- results$FeatureTable
    isf_featuerTable <- cbind(customFT[,1],result)
    colnames(isf_featuerTable)[1] <- "featureID"
    col_index <- which( colnames(isf_featuerTable) %in% c(colnames(peaks), "cl" ,"ISF_level"))
    result <- isf_featuerTable[,col_index]
    result <- result[result$cl != 0,]
    
    file.remove(paste0(ISFdirectory_name, "/",mzMLfile[i]))
    write.csv(result, paste0(ISFdirectory_name,"/isf_results.csv"), row.names = FALSE)
    
  }else{result <- 0}
  cl_tb_POS <- cl_tb[cl_tb$cl !=0,]
  
  # Identify salt adducts, isotopes
  cl_tb_cleaned <- ChloroDBPFinder::cleanFeature(peaks = peaks, chlorine_tb = cl_tb_POS,
                                                 rawlcms = xcmsrawlcms, rawfile_dir = mzmldir, lcmsfile = mzMLfile[i],
                                                 adducts = TRUE, isotopes = TRUE, inSourceFrag = isfrag,
                                                 ISFtable = result,
                                                 samNum = i)
  # Output the table of chlorine-containing features
  write.csv(cl_tb_cleaned, paste0(mzmldir,"/",strsplit(mzMLfile[i], split = lcmspattern)[[1]][1], "_", nrow(cl_tb_cleaned), "_cl.csv"),row.names = FALSE )
  high_quality_cl_tb <- cl_tb_cleaned[cl_tb_cleaned$Adduct == 0 & cl_tb_cleaned$isotope == 0 & cl_tb_cleaned$ISF_level == 0,]
  write.csv(high_quality_cl_tb, paste0(mzmldir,"/",strsplit(mzMLfile[i], split = lcmspattern)[[1]][1], "_", nrow(high_quality_cl_tb), "_cl_high_quality.csv"), row.names = FALSE)
}

### --------- 2. Alignment across samples  ----------------------
aligned_tb <- alignFeature(file_dir = mzmldir, filePattern = "_cl_high_quality.csv")
write.csv(aligned_tb, paste0(mzmldir, "/",nrow(aligned_tb),"_alignment.csv"), row.names = FALSE)

### --------- 3. Evidence-based missing value imputation  ---------------
filled_tb <- fillGap(file_dir = mzmldir, mzmlfiles_pattern = lcmspattern, aligned_tb = aligned_tb )
write.csv(filled_tb, paste0(mzmldir, "/gap_filled.csv"), row.names = FALSE )

###------------- 4. annotation -------------------------------
#### ----------- 4.1 Spectral database search ------------------------
# pls remember to specify the path of feature table that need annotation
# users need to manual change the path within read.csc() function below
table_need_annotation <- read.csv("C:/Users/User/Desktop/testmzML/TW APM noAscorbic_159_cl_high_quality.csv")

# Load database
Cl_db <-  read.csv(Cl_db_path)

# Compound annotation by spectral database searching
annotated_tb <- annonateFeature(featureTable = table_need_annotation, Cl_db, ion_mode = "P", ref_mz_tol =25, dp_score = 70, dp_num = 2)
# Output annotations
write.csv(annotated_tb, paste0(mzmldir, "/", nrow(annotated_tb[annotated_tb$score!=0,]),"_annotations.csv"), row.names = FALSE)

###------------- 4.2 molecular networking -------------------- 
# Load the reference table containing known compounds or precursors
#data('seed_demo')
reference_table <- read.csv(ref_path)

# load reaction pathway
data("reaction_pathways")

# construct the network
network <- createNetwork(featureTable = annotated_tb , seed_tb = reference_table,
                         reaction_pathways = reaction_pathways, reaction_pathway_mz_tol = 0.01,
                         nw_spectra_score=0.5, nw_spectra_match_num=3)

# Output reaction and spectral networks
all_nw <- network[[1]]
write.csv(all_nw, paste0(mzmldir, "/molecular_networks.csv"), row.names = FALSE)

# Output network with explainable connections which have high spectral connection and reaction pathway connection
integrated_nw <- network[[2]]
write.csv(integrated_nw, paste0(mzmldir, "/integrated_molecular_network.csv"), row.names = FALSE)
