#' @title AllCcsPrediction
#' @author Zhiwei Zhou
#' @param mol_smiles
#' @param mol_names
#' @importFrom magrittr '%>%'
#' @importClassesFrom e1071 'svm'
#' @export
#' @examples

# test <- read.csv('./inst/extdata/demo_data_190523.csv', stringsAsFactors = F)
# test <- AllCcsPrediction(mol_smiles = test$smiles, mol_names = test$id, base_dir = 'F:/01 MetIMMS/00 data processing/190523 demo data of allccs')

setGeneric(name = 'AllCcsPrediction',
           def = function(
             mol_smiles,
             mol_names,
             base_dir = './AllCCS_result',
             is_output = TRUE
           ){

             library(e1071)
             # dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
             if (length(mol_smiles)!=length(mol_names)) {
               stop('Please check the number of mol_smiles and mol_names.\n')
             }

             # MD calculation
             # if exist calculated MD, direct load calculated MDs
             cat('Calculate molecular descriptors...\n\n')

            if (!file.exists(file.path(base_dir, '00_intermediate_data', 'calculated_md.RData'))) {
              if (length(mol_smiles) > 1) {
                md_result <- pbapply::pblapply(seq(length(mol_smiles)), function(i){
                  MdCalculation(mol_smiles = mol_smiles[i], mol_names = mol_names[i])
                })

                md_result <- do.call(rbind.data.frame, md_result)

              } else {
                md_result <- MdCalculation(mol_smiles = mol_smiles, mol_names = mol_names)
              }


              dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
              save(md_result,
                   file = file.path(base_dir, '00_intermediate_data', 'calculated_md.RData'))
            } else {
              cat('Detected calculated molecular descriptors, and load it\n\n')
              load(file.path(base_dir, '00_intermediate_data', 'calculated_md.RData'))
            }


             # check MD range to evaluate whether it is applicable
             cat('Check the range of molecular descriptors...\n\n')
             if (!all(md_result$MW >= 60 & md_result$MW <= 1200)) {
               temp_idx <- which(md_result$MW < 60 | md_result$MW > 1200)

               if (length(temp_idx)==nrow(md_result)) {
                 stop('All compounds EM out of applicable range (60-1200)\n' )
               }

               warning(paste0(length(temp_idx), ' compounds EM out of applicable range (60-1200). These compounds were removed in the final table\n'))

               md_result <- md_result[-temp_idx,,drop=F]

             }


             # calculate m/z
             trans_md_pos <- lapply(seq(nrow(md_result)), function(i){
               result <- mz_transform(M = md_result$MW[i],
                                     adduct = c('[M+H]+', '[M+H-H2O]+', '[M+Na]+', '[M+NH4]+'),
                                     polarity = 'pos')

               options(warn = -1)
               result <- data.frame(name=md_result$name[i],
                                    adduct=result$adduct,
                                    mz=result$mz,
                                    md_result[i, -1, drop=F],
                                    stringsAsFactors = F)

               rownames(result) <- NULL
               result
             })

             trans_md_neg <- lapply(seq(nrow(md_result)), function(i){
               result <- mz_transform(M = md_result$MW[i],
                                      adduct = c('[M-H]-', '[M+Na-2H]-', '[M+HCOO]-'),
                                      polarity = 'neg')

               options(warn = -1)
               result <- data.frame(name=md_result$name[i],
                                    adduct=result$adduct,
                                    mz=result$mz,
                                    md_result[i, -1, drop=F],
                                    stringsAsFactors = F)

               rownames(result) <- NULL
               result
             })

             trans_md_pos <- do.call(rbind.data.frame, trans_md_pos)
             trans_md_neg <- do.call(rbind.data.frame, trans_md_neg)

             cat('Predict CCS values...\n\n')
             # MD imputation and scaling ---------------------------------------
             # cat('Is: here1\n\n')
             load(system.file("extdata",
                              "md_optimized_pos_standardization_table_190425.RData",
                              package="AllCCS"))

             load(system.file("extdata",
                              "md_optimized_neg_standardization_table_190427.RData",
                              package="AllCCS"))

             # cat('Is: here2\n\n')
             trans_md_imputed_pos <- quiet(
               AllCcsImputation(trans_md = trans_md_pos,
                                md_standardization_table = md_optimzed_pos_standardization_table,
                                polarity = 'pos')
             )

             trans_md_imputed_neg <- quiet(
               AllCcsImputation(trans_md = trans_md_neg,
                                md_standardization_table = md_optimzed_neg_standardization_table,
                                polarity = 'neg')
             )


             trans_md_zscore_pos <- AllCcsZscore(trans_md_imputed = trans_md_imputed_pos,
                                                 md_standardization_table = md_optimzed_pos_standardization_table)
             trans_md_zscore_neg <- AllCcsZscore(trans_md_imputed = trans_md_imputed_neg,
                                                 md_standardization_table = md_optimzed_neg_standardization_table)

             # ---------------------------------------------------------------
             load(system.file("extdata",
                              "svr_reg_pos_model_10times_190430.RData",
                              package="AllCCS"))

             load(system.file("extdata",
                              "svr_reg_neg_model_10times_190430.RData",
                              package="AllCCS"))


             predicted_ccs_pos <- predict(svr_reg_pos, trans_md_zscore_pos[,-c(1:2)])
             predicted_ccs_neg <- predict(svr_reg_neg, trans_md_zscore_neg[,-c(1:2)])

             result_pos <- data.frame(name = trans_md_zscore_pos$name,
                                      exact_mass = rep(md_result$MW, each=4),
                                      adduct = trans_md_zscore_pos$adduct,
                                      mz = trans_md_pos$mz,
                                      pred_CCS = predicted_ccs_pos,
                                      stringsAsFactors = F)

             result_neg <- data.frame(name = trans_md_zscore_neg$name,
                                      exact_mass = rep(md_result$MW, each=3),
                                      adduct = trans_md_zscore_neg$adduct,
                                      mz = trans_md_neg$mz,
                                      pred_CCS = predicted_ccs_neg,
                                      stringsAsFactors = F)

             final_result <- rbind.data.frame(result_pos, result_neg)

             if (is_output) {
               readr::write_csv(x = final_result,
                                path = file.path(base_dir, 'pred_ccs_result.csv'))
             }

             cat('The prediction has completed !\n')

             return(final_result)
             # change to wide table

             # final_result_mz <- reshape2::dcast(final_result,
             #                                    name + exact_mass ~ adduct,
             #                                    value.var = c('mz'))
             #
             # final_result_ccs <- reshape2::dcast(final_result,
             #                                     name + exact_mass ~ adduct,
             #                                     value.var = c('pred_CCS'))
             #
             # final_result_ccs
           }
)


setGeneric(name = 'AllCcsImputation',
           def = function(
             trans_md,
             md_standardization_table,
             polarity = c('pos', 'neg')
           ){

             polarity <- match.arg(polarity)

             switch (polarity,
                     'pos' = {
                       # data('hmdb_mds_pos.rda', envir = environment())
                       # load(data('data/hmdb_mds_pos.rda', package = 'AllCCS'))
                       hmdb_mds <- hmdb_mds_pos
                       rm(hmdb_mds_pos)
                       gc()
                     },
                     'neg' = {
                       # data('hmdb_mds_neg.rda', envir = environment())
                       # load(data('data/hmdb_mds_neg.rda', package = 'AllCCS'))
                       hmdb_mds <- hmdb_mds_neg
                       rm(hmdb_mds_neg)
                       gc()
                     }
             )

             select_md_name <- rownames(md_standardization_table)
             temp_md <- trans_md[,match(select_md_name[-1], colnames(trans_md))]
             idx_input <- seq(nrow(temp_md)) + nrow(hmdb_mds)

             # integrate with HMDB for imputation
             temp_md <- rbind.data.frame(hmdb_mds, temp_md)
             temp_md <- ZZWImpute(raw_matrix = temp_md, label = seq(nrow(temp_md)))

             trans_md_imputed <- temp_md[idx_input, -1]
             rownames(trans_md_imputed) <- NULL
             trans_md_imputed <- data.frame(trans_md[,1:3], trans_md_imputed)

           }
)

setGeneric(name = 'AllCcsZscore',
           def = function(
             trans_md_imputed,
             md_standardization_table
           ){

             temp_idx <- match(rownames(md_standardization_table), colnames(trans_md_imputed))

             if (any(is.na(temp_idx))) {
               stop('The input md table is not matched with standardization table')
             }
             temp_data <- trans_md_imputed[, temp_idx]

             result <- lapply(seq(ncol(temp_data)), function(i){
               temp_mean <- md_standardization_table$mean[i]
               temp_sd <- md_standardization_table$sd[i]

               result <- (temp_data[,i]-temp_mean)/temp_sd
             })

             result <- do.call(cbind, result)
             colnames(result) <- colnames(temp_data)

             result <- data.frame(name=trans_md_imputed$name,
                                  adduct=trans_md_imputed$adduct,
                                  result,
                                  stringsAsFactors = F)

             result
           }
)

