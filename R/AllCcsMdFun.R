#' @title MdGet
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @description Calculate molecular descriptors (MD), and perform MD cleaning (i.e. m/z transformation, NA imputation & scaling)
#' @param mol_smiles a vector of smiles
#' @param mol_names a vector of names
#' @param thread the number of calculating cores. Default: 2
#' @param base_dir '.'
#' @return a list of MDs. "md_pos": Molecular descriptors in positive mode; "md_neg": Molecular descriptors in negative mode
#' @export
#' @examples
#' test <- read.csv(system.file("extdata", "demo_data.csv", package="AllCCS"), stringsAsFactors = F)
#’ test <- MdGet(mol_smiles = test$smiles, mol_names = test$id_allccs, base_dir = '.')

# test <- read.csv('./inst/extdata/demo_data.csv', stringsAsFactors = F)
# test <- MdGet(mol_smiles = test$smiles, mol_names = test$id_allccs, base_dir = '.')

setGeneric(name = 'MdGet',
           def = function(
             mol_smiles,
             mol_names,
             status_label,
             thread = 2,
             base_dir = '.'
           ){
             idx_valid <- (status_label == 'Valid') %>% which()
             idx_error <- (status_label == 'Error_1') %>% which()

             mol_name_valid <- mol_names[idx_valid]
             mol_name_error <- mol_names[idx_error]

             mol_smiles_valid <- mol_smiles[idx_valid]
             mol_smiles_error <- mol_smiles[idx_error]
             names(mol_smiles_valid) <- mol_name_valid
             names(mol_smiles_error) <- mol_name_error

             status_label_valid <- status_label[idx_valid]
             status_label_error <- status_label[idx_error]
             names(status_label_valid) <- mol_name_valid
             names(status_label_error) <- mol_name_error


             # check MD range to evaluate whether it is applicable
             cat('Check the mass range...\n\n')

             if (!('mass_valid.RData' %in% list.files(file.path(base_dir,
                                                               '00_intermediate_data')))) {
               mass_valid <- Smiles2Mass(mol_smiles = mol_smiles_valid,
                                         mol_names = mol_name_valid)

               dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
               save(mass_valid,
                    file = file.path(base_dir, '00_intermediate_data', 'mass_valid.RData'),
                    compress = 'gzip')

             } else {
               cat('Detected calculated mass_valid, and load it\n\n')
               load(file.path(base_dir, '00_intermediate_data', 'mass_valid.RData'))
             }


             if (!all(mass_valid >= 60 & mass_valid <= 1200)) {
               temp_idx <- which(mass_valid < 60 | mass_valid > 1200)


               # if (length(temp_idx)==nrow(result_md)) {
               #   stop('All compounds EM out of applicable range (60-1200)\n' )
               # }

               # warning(paste0(length(temp_idx), ' compounds EM out of applicable range (60-1200). These compounds were removed in the final table\n'))

               # result_md <- result_md[-temp_idx,,drop=F]

               mol_name_error <- c(mol_name_error, mol_name_valid[temp_idx])
               mol_name_valid <- mol_name_valid[-temp_idx]

               mol_smiles_error <- c(mol_smiles_error, mol_smiles_valid[temp_idx])
               mol_smiles_valid <- mol_smiles_valid[-temp_idx]
               names(mol_smiles_error) <- mol_name_error
               names(mol_smiles_valid) <- mol_name_valid

               status_label_error <- c(status_label_error,
                                       rep('Error_2', length(temp_idx)))
               status_label_valid <- status_label_valid[-temp_idx]
               names(status_label_error) <- mol_name_error
               names(status_label_valid) <- mol_name_valid

               # if all smiles are error, return meta_error
               if (length(status_label_error) == length(status_label)) {
                 cat('All compounds encounter errors\n\n' )
                 meta_error <- data.frame(name = mol_name_error,
                                          SMILES = mol_smiles_error,
                                          monoisotopic_mass = NA,
                                          adduct = NA,
                                          mz = NA,
                                          status = status_label_error,
                                          stringsAsFactors = F)

                 result <- list(md_pos = NULL,
                                md_neg = NULL,
                                meta_pos = NULL,
                                meta_neg = NULL,
                                meta_error = meta_error)

                 return(result)
               }
             }




             # MD calculation --------------------------------------------------
             # if exist calculated MD, direct load calculated MDs
             cat('Calculate molecular descriptors...\n\n')

             if (!('result_md.RData' %in% list.files(file.path(base_dir,
                                                                '00_intermediate_data')))) {
               if (length(mol_smiles_valid) > 1) {
                 temp_fun <- function(i, mol_smiles_valid, mol_name_valid, MdCalculate) {
                   MdCalculate(mol_smiles = mol_smiles_valid[i], mol_names = mol_name_valid[i])
                 }

                 result_md <- BiocParallel::bplapply(X = seq_along(mol_smiles_valid),
                                                     FUN = temp_fun,
                                                     BPPARAM = BiocParallel::SnowParam(workers = thread,
                                                                                       progressbar = TRUE),
                                                     mol_smiles_valid=mol_smiles_valid,
                                                     mol_name_valid=mol_name_valid,
                                                     MdCalculate=MdCalculate)

                 # result_md <- pbapply::pblapply(seq(length(mol_smiles)), function(i){
                 #   MdCalculate(mol_smiles = mol_smiles[i], mol_names = mol_names[i])
                 # })

                 result_md <- result_md %>% dplyr::bind_rows()

               } else {
                 result_md <- MdCalculation(mol_smiles = mol_smiles_valid,
                                            mol_names = mol_name_valid)
               }


               dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
               save(result_md,
                    file = file.path(base_dir, '00_intermediate_data', 'result_md.RData'))
             } else {
               cat('Detected calculated molecular descriptors, and load it\n\n')
               load(file.path(base_dir, '00_intermediate_data', 'result_md.RData'))

               result_md <- match(mol_name_valid, result_md$name) %>%
                 result_md[.,]
             }




             # transform to mz -------------------------------------------------
             cat('Calculate m/z...', '\n')
             result_md_trans_pos <- match(c('name', 'MW', colnames(md_hmdb_pos)), colnames(result_md)) %>%
               result_md[,.]
             result_md_trans_neg <- match(c('name', 'MW', colnames(md_hmdb_neg)), colnames(result_md)) %>%
               result_md[,.]

             rm(result_md);gc()

             result_md_trans_pos <- pbapply::pblapply(seq_along(result_md_trans_pos$name), function(i){
               result <- MzTransform(M = result_md_trans_pos$MW[i],
                                     adduct = c('[M+H]+', '[M+H-H2O]+', '[M+Na]+', '[M+NH4]+'))

               options(warn = -1)
               result <- data.frame(name=result_md_trans_pos$name[i],
                                    adduct=result$adduct,
                                    mz=result$mz,
                                    result_md_trans_pos[i, -1, drop=F],
                                    stringsAsFactors = F)

               rownames(result) <- NULL
               result
             }) %>%
               dplyr::bind_rows()


             result_md_trans_neg <- pbapply::pblapply(seq_along(result_md_trans_neg$name), function(i){
               result <- MzTransform(M = result_md_trans_neg$MW[i],
                                     adduct = c('[M-H]-', '[M+Na-2H]-', '[M+HCOO]-'))

               options(warn = -1)
               result <- data.frame(name=result_md_trans_neg$name[i],
                                    adduct=result$adduct,
                                    mz=result$mz,
                                    result_md_trans_neg[i, -1, drop=F],
                                    stringsAsFactors = F)

               rownames(result) <- NULL
               result
             }) %>%
               dplyr::bind_rows()

             meta_pos <- result_md_trans_pos %>%
               dplyr::mutate(smiles = match(result_md_trans_pos$name, names(mol_smiles_valid)) %>%
                               mol_smiles_valid[.],
                             status = match(result_md_trans_pos$name, names(status_label_valid)) %>%
                               status_label_valid[.]) %>%
               dplyr::select(name, smiles, MW, adduct, mz, status) %>%
               dplyr::rename(SMILES = smiles,
                             monoisotopic_mass = MW)

             meta_neg <- result_md_trans_neg %>%
               dplyr::mutate(smiles = match(result_md_trans_neg$name, names(mol_smiles_valid)) %>%
                               mol_smiles_valid[.],
                             status = match(result_md_trans_neg$name, names(status_label_valid)) %>%
                               status_label_valid[.]) %>%
               dplyr::select(name, smiles, MW, adduct, mz, status) %>%
               dplyr::rename(SMILES = smiles,
                             monoisotopic_mass = MW)

             if (length(status_label_error) > 0) {
               meta_error <- data.frame(name = mol_name_error,
                                        SMILES = mol_smiles_error,
                                        monoisotopic_mass = NA,
                                        adduct = NA,
                                        mz = NA,
                                        status = status_label_error,
                                        stringsAsFactors = F)
             } else {
               meta_error <- NULL
             }


             # transform to mz -------------------------------------------------
             cat('\n', 'Impute NAs...\n\n', sep = '')

             result_md_trans_pos <- quiet(
               NaImputate(result_md_trans = result_md_trans_pos,
                          polarity = 'pos')
             )

             result_md_trans_neg <- quiet(
               NaImputate(result_md_trans = result_md_trans_neg,
                          polarity = 'neg')
             )

             cat('Z-score scaling...\n\n')

             result_md_zscore_pos <- ZscoreScaling(result_md_trans = result_md_trans_pos,
                                                   polarity = 'pos')

             result_md_zscore_neg <- ZscoreScaling(result_md_trans = result_md_trans_neg,
                                                   polarity = 'neg')


             cat('\n', 'Molecular descriptors calculation has completed!\n\n', sep = '')
             result <- list(md_pos = result_md_zscore_pos,
                            md_neg = result_md_zscore_neg,
                            meta_pos = meta_pos,
                            meta_neg = meta_neg,
                            meta_error = meta_error)

             return(result)


           }
)




#' @title MdCalculate
#' @author Zhiwei Zhou
#' @param mol_smiles a smiles structure of molecule
#' @param mol_names the name of structure
#' @export
#' @examples
#' MdCalculate(mol_smiles = "COC1=CC=C(CC(O)=O)C=C1", mol_names = "4-Methoxyphenylacetic Acid")

# MdCalculate(mol_smiles = "COC1=CC=C(CC(O)=O)C=C1", mol_names = "4-Methoxyphenylacetic Acid")
setGeneric(name = 'MdCalculate',
           def = function(
             mol_smiles,
             mol_names
           ){

             if (is.na(mol_names) | is.null(mol_names)) {
               mol_names <- NA
             }

             if (is.na(mol_smiles) | is.null(mol_smiles)) {
               itp_smiles <- rcdk::parse.smiles("COC1=CC=C(CC(O)=O)C=C1")
               desc_names <- rcdk::get.desc.names(type = "all") #descriptors name
               descriptors <- rcdk::eval.desc(itp_smiles, desc_names)

               result <- rep(NA, 286) %>%
                 matrix(nrow = 1) %>%
                 as.data.frame()

               colnames(result) <- colnames(descriptors)

               result <- result %>%
                 dplyr::mutate(name = NA) %>%
                 dplyr::select(name, dplyr::everything())

               return(result)
             }

             itp_smiles <- rcdk::parse.smiles(mol_smiles)
             desc_names <- rcdk::get.desc.names(type = "all") #descriptors name
             descriptors <- rcdk::eval.desc(itp_smiles, desc_names)

             result <- data.frame(name=mol_names, descriptors, stringsAsFactors = F)
             rownames(result) <- NULL
             return(result)

           }
)



#' @title NaImpute
#' @author Zhiwei Zhou
#' @param result_md_trans
#' @param md_standardization_table
#' @param polarity Default: pos

# NaImputate(result_md_trans = result_md_trans_pos,
#            polarity = 'pos')

setGeneric(name = 'NaImputate',
           def = function(
             result_md_trans,
             polarity = c('pos', 'neg')
           ){

             polarity <- match.arg(polarity)

             switch (polarity,
                     'pos' = {
                       # load standardization table
                       # load(system.file("extdata",
                       #                  "md_optimized_pos_sd_table_190812.RData",
                       #                  package="AllCCS"))

                       md_hmdb <- md_hmdb_pos

                       md_standardization_table <- md_optimzed_pos_sd_table

                       # rm(md_hmdb_pos, md_optimzed_pos_sd_table);gc()
                       # rm(md_hmdb_pos, md_optimzed_pos_sd_table);gc()
                     },
                     'neg' = {
                       # load(system.file("extdata",
                       #                  "md_optimized_neg_sd_table_190812.RData",
                       #                  package="AllCCS"))

                       md_standardization_table <- md_optimzed_neg_sd_table

                       md_hmdb <- md_hmdb_neg
                       # rm(md_hmdb_neg, md_optimzed_neg_sd_table);gc()
                     }
             )



             select_md_name <- rownames(md_standardization_table)
             temp_md <- result_md_trans[,match(select_md_name[-1], colnames(result_md_trans))]
             idx_input <- seq(nrow(temp_md)) + nrow(md_hmdb)

             # integrate with HMDB for imputation
             temp_md <- md_hmdb %>% dplyr::bind_rows(temp_md)
             temp_md <- KnnImpute(raw_matrix = temp_md, label = seq(nrow(temp_md)))

             result_md_trans_imputed <- temp_md[idx_input, -1]
             rownames(result_md_trans_imputed) <- NULL
             result_md_trans_imputed <- data.frame(result_md_trans[,1:3],
                                                   result_md_trans_imputed,
                                                   stringsAsFactors = F)

           }
)




#' @title ZscoreScaling
#' @description z-score transformation for MDs of new structure
#' @author Zhiwei Zhou
#' @param result_md_trans Mds after mz transformation & NA imputation
#' @param polarity
#' @return a data.frame of MDs for new structures
#' @examples
#'

# ZscoreScaling(result_md_trans = result_md_trans_pos,
#               polarity = 'pos')

setGeneric(name = 'ZscoreScaling',
           def = function(
             result_md_trans,
             polarity = c('pos', 'neg')
           ){

             polarity <- match.arg(polarity)

             switch (polarity,
                     'pos' = {
                       # load standardization table
                       load(system.file("extdata",
                                        "md_optimized_pos_sd_table_190812.RData",
                                        package="AllCCS"))

                       md_standardization_table <- md_optimzed_pos_sd_table

                       rm(md_optimzed_pos_sd_table);gc()
                     },
                     'neg' = {
                       load(system.file("extdata",
                                        "md_optimized_neg_sd_table_190812.RData",
                                        package="AllCCS"))

                       md_standardization_table <- md_optimzed_neg_sd_table

                       rm(md_optimzed_neg_sd_table);gc()
                     }
             )


             temp_data <- match(rownames(md_standardization_table), colnames(result_md_trans)) %>%
               result_md_trans[,.]

             result <- pbapply::pblapply(seq(ncol(temp_data)), function(i){
               temp_mean <- md_standardization_table$mean[i]
               temp_sd <- md_standardization_table$sd[i]

               result <- (temp_data[,i]-temp_mean)/temp_sd
             })

             # result <- do.call(cbind, result)
             result <- result %>% dplyr::bind_cols()

             colnames(result) <- colnames(temp_data)

             result <- data.frame(name=result_md_trans$name,
                                  adduct=result_md_trans$adduct,
                                  result,
                                  stringsAsFactors = F)

             result
           }
)




