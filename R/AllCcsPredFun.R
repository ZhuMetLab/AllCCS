#' @title CcsPredict
#' @description A function to predict CCS values for small molecules with the input of SMILES structures
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param mol_smiles a vector of SMILES structures
#' @param mol_names a vector of identifier
#' @param base_dir '.'
#' @param thread Default: 3
#' @param is_output Default: TRUE
#' @importFrom magrittr '%>%' '%$%'
#' @export
#' @examples
#' test <- read.csv(system.file("extdata", "demo_data.csv", package="AllCCS"), stringsAsFactors = F)
#' test <- MdGet(mol_smiles = test$smiles, mol_names = test$id_allccs, base_dir = '.')


# test <- read.csv('./inst/extdata/demo_data.csv', stringsAsFactors = F)
# test <- CcsPredict(mol_smiles = test$smiles,
#                    mol_names = test$id,
#                    thread = 2,
#                    base_dir = '.',
#                    is_output = F)

# test <- read.csv('./inst/extdata/demo_data_190902.csv', stringsAsFactors = F)
# test <- CcsPredict(mol_smiles = test$smiles,
#                    mol_names = test$id,
#                    thread = 1,
#                    base_dir = '.',
#                    is_output = F)

setGeneric(name = 'CcsPredict',
           def = function(
             mol_smiles,
             mol_names,
             base_dir = '.',
             thread = 3,
             is_output = TRUE
           ){

             library(e1071)
             # dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
             if (length(mol_smiles)!=length(mol_names)) {
               stop('Please check the number of mol_smiles and mol_names.\n')
             }

             if (length(unique(mol_names))!=length(mol_names)) {
               stop('Please check the inputed names (Replicates existed).\n')
             }

             cat('Check smiles validity...\n\n')
             if (!('status_label.RData' %in% list.files(file.path(base_dir,
                                                                  '00_intermediate_data')))) {
               status_label <- SmilesCheck(mol_smiles)

               dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
               save(status_label,
                    file = file.path(base_dir, '00_intermediate_data', 'status_label.RData'))
             } else {
               dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
               load(file.path(base_dir, '00_intermediate_data', 'status_label.RData'))
             }

             names(mol_smiles) <- mol_names
             names(status_label) <- mol_names

             if (sum(status_label == 'Error_1') == length(mol_names)) {
               temp <- data.frame(name = mol_names,
                                  SMILES = mol_smiles,
                                  monoisotopic_mass = NA,
                                  adduct = NA,
                                  mz = NA,
                                  ccs = NA,
                                  rss = NA,
                                  status = status_label,
                                  stringsAsFactors = F)
               return(temp)
             }


             # MD calculation
             # if exist calculated MD, direct load calculated MDs
             cat('------------------------------------------------------------\n')
             if (!('result_md_final.RData' %in% list.files(file.path(base_dir,
                                                                     '00_intermediate_data')))) {
               result_md_final <- MdGet(mol_smiles = mol_smiles,
                                        mol_names = mol_names,
                                        thread = thread,
                                        status_label = status_label,
                                        base_dir = base_dir)

               dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
               save(result_md_final,
                    file = file.path(base_dir, '00_intermediate_data', 'result_md_final.RData'))

             } else {
               dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
               load(file.path(base_dir, '00_intermediate_data', 'result_md_final.RData'))
             }




             result_md_pos <- result_md_final$md_pos
             result_md_neg <- result_md_final$md_neg

             meta_pos <- result_md_final$meta_pos
             meta_neg <- result_md_final$meta_neg

             meta_error <- result_md_final$meta_error

             rm(result_md_final);gc()

             if (length(meta_pos)<1) {
               final_result <- meta_error %>%
                 dplyr::mutate(pred_ccs = NA,
                               rss = NA) %>%
                 dplyr::select(name:mz, pred_ccs, rss, status) %>%
                 dplyr::arrange(name)

               return(final_result)
             }


             # ---------------------------------------------------------------
             cat('\n'); cat('------------------------------------------------------------\n')
             cat('Predict CCS values ...\n\n')
             load(system.file("extdata",
                              "svr_reg_pos_model_100times_190816.RData",
                              package="AllCCS"))

             load(system.file("extdata",
                              "svr_reg_neg_model_100times_190813.RData",
                              package="AllCCS"))


             pred_ccs_pos <- predict(svr_reg_pos, result_md_pos[,-c(1:2)])
             pred_ccs_neg <- predict(svr_reg_neg, result_md_neg[,-c(1:2)])

             rm(result_md_pos, result_md_neg); gc()


             cat('------------------------------------------------------------\n')
             if (!('rss_pos.RData' %in% list.files(file.path(base_dir,
                                                                '00_intermediate_data')))) {
               rss_pos <- RssCalculate(mol_smiles = match(meta_pos$name, mol_names) %>%
                                         unique() %>%
                                         mol_smiles[.],
                                       mol_names = match(meta_pos$name, mol_names) %>%
                                         unique() %>%
                                         mol_names[.],
                                       thread = thread,
                                       max_n = 5,
                                       polarity = 'pos',
                                       type = 'pubchem',
                                       method = 'tanimoto',
                                       base_dir = base_dir)

               dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
               save(rss_pos,
                    file = file.path(base_dir, '00_intermediate_data', 'rss_pos.RData'))
             } else {
               dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
               load(file.path(base_dir, '00_intermediate_data', 'rss_pos.RData'))
             }


             if (!('rss_neg.RData' %in% list.files(file.path(base_dir,
                                                             '00_intermediate_data')))) {
               rss_neg <- RssCalculate(mol_smiles = match(meta_neg$name, mol_names) %>%
                                         unique() %>%
                                         mol_smiles[.],
                                       mol_names = match(meta_neg$name, mol_names) %>%
                                         unique() %>%
                                         mol_names[.],
                                       thread = thread,
                                       max_n = 5,
                                       polarity = 'neg',
                                       type = 'pubchem',
                                       method = 'tanimoto',
                                       base_dir = base_dir)

               dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
               save(rss_neg,
                    file = file.path(base_dir, '00_intermediate_data', 'rss_neg.RData'))
             } else {
               dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
               load(file.path(base_dir, '00_intermediate_data', 'rss_neg.RData'))
             }

             rss_pos <- match(meta_pos$name, names(rss_pos)) %>%
               rss_pos[.] %>%
               round(digits = 4)

             rss_neg <- match(meta_neg$name, names(rss_neg)) %>%
               rss_neg[.] %>%
               round(digits = 4)

             cat('\n', '------------------------------------------------------------\n',
                 sep = '')
             cat('Export predicted result ...\n\n')

             result_pos <- meta_pos %>%
               dplyr::mutate(pred_ccs = pred_ccs_pos,
                             rss = rss_pos,
                             name = as.character(name)) %>%
               dplyr::select(name:mz, pred_ccs, rss, status)

             result_neg <- meta_neg %>%
               dplyr::mutate(pred_ccs = pred_ccs_neg,
                             rss = rss_neg,
                             name = as.character(name)) %>%
               dplyr::select(name:mz, pred_ccs, rss, status)

             final_result <- dplyr::bind_rows(result_pos, result_neg) %>%
               dplyr::arrange(name)

             if (length(meta_error) > 0) {
               result_error <- meta_error %>%
                 dplyr::mutate(pred_ccs = NA,
                               rss = NA,
                               name = as.character(name)) %>%
                 dplyr::select(name:mz, pred_ccs, rss, status) %>%
                 dplyr::arrange(name)

               final_result <- dplyr::bind_rows(result_pos, result_neg) %>%
                 dplyr::bind_rows(., result_error) %>%
                 dplyr::arrange(name)

             }


             if (is_output) {
               readr::write_csv(x = final_result,
                                path = file.path(base_dir, 'pred_ccs_result.csv'))
             }

             cat('The prediction has completed!\n')

             return(final_result)

           }
)


.onAttach <- function(libname, pkgname){
  packageStartupMessage("
More information can be found in http://imms.zhulab.cn.
If you have any questions, please send email to zhouzw@sioc.ac.cn or jiangzhu@sioc.ac.cn.
Authors: Zhiwei Zhou and Dr. Zhengjiang Zhu (jiangzhu@sioc.ac.cn).
Maintainer: Zhiwei Zhou
Version 0.1.4 (20190909)
--------------
o Optimize frame to save memory.
o Fix some bugs.")
}


