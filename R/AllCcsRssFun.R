#' @title RssCalculate
#' @description Calculate representative structure similarity
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param mol_smiles a vector of SMILES structures
#' @param mol_names a vector of identifier
#' @param thread Default: 2
#' @param max_n Top n for calculate RSS score
#' @param polarity 'pos', 'neg'
#' @param type The type of fingerprint. See \code{\link{get.fingerprint}}.
#' @param method The type of distance metric to use. See \code{\link{fp.sim.matrix}}.
#' @param base_dir
#' @return a vector of rss score
#' @export
#' @examples
#' test <- RssCalculate(mol_smiles = mol_smiles, mol_names = mol_names, thread = 2, max_n = 5, polarity = 'pos')

# test <- RssCalculate(mol_smiles = mol_smiles,
#                      mol_names = mol_names,
#                      thread = 2,
#                      max_n = 5,
#                      polarity = 'pos')

setGeneric(name = 'RssCalculate',
           def = function(
             mol_smiles,
             mol_names,
             thread = 2,
             max_n=5,
             polarity = c('pos', 'neg'),
             type = 'pubchem',
             method = 'tanimoto',
             base_dir = '.'){


             result_fps <- FpCalculate(mol_smiles = mol_smiles,
                                       mol_names = mol_names,
                                       type = type,
                                       thread = thread,
                                       base_dir = base_dir)

             switch (polarity,
                     'pos' = {
                       # load standardization table
                       load(system.file("extdata",
                                        "fps_training_pos_190819.RData",
                                        package="AllCCS"))

                       fplist <- fps_training_pos

                       rm(fps_training_pos);gc()
                     },
                     'neg' = {
                       load(system.file("extdata",
                                        "fps_training_neg_190819.RData",
                                        package="AllCCS"))

                       fplist <- fps_training_neg

                       rm(fps_training_neg);gc()
                     }
             )


             cat('Calculate Representive structure similarity ...\n\n')
             temp_fun <- function(i, result_fps, fplist, method){
               temp <- result_fps[i]
               result <- fingerprint::fp.sim.matrix(fplist = fplist, fplist2 = temp, method = method)

               result <- sort(result, decreasing = TRUE) %>%
                 .[1:max_n] %>%
                 mean()

               return(result)
             }

             rss_score <- BiocParallel::bplapply(X = seq_along(result_fps),
                                                 FUN = temp_fun,
                                                 BPPARAM = BiocParallel::SnowParam(workers = thread, progressbar = TRUE),
                                                 result_fps = result_fps,
                                                 fplist = fplist,
                                                 method = method)

             dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
             save(rss_score,
                  file = file.path(base_dir, '00_intermediate_data', 'rss_score.RData'))

             rss_score <- rss_score %>% unlist()

             # rss_score <- pbapply::pbsapply(seq_along(result_fps), function(i){
             #   # cat(i);cat(' ')
             #   temp <- result_fps[i]
             #   result <- fingerprint::fp.sim.matrix(fplist = fplist, fplist2 = temp, method = method)
             #
             #   result <- sort(result, decreasing = TRUE) %>%
             #     .[1:max_n] %>%
             #     mean()
             #
             #   return(result)
             # })

             names(rss_score) <- mol_names

             return(rss_score)
           }
)


#' @title FpCalcualte
#' @author Zhiwei Zhou
#' @param mol_smiles
#' @param mol_names
#' @param thread Default: 2
#' @param base_dir '.'
#' @export
#' @examples
#' result_fps <- FpCalculate(mol_smiles = mol_smiles,
#'                           mol_names = mol_names,
#'                           thread = 2)


# result_fps <- FpCalculate(mol_smiles = mol_smiles,
#                           mol_names = mol_names,
#                           thread = 2)
setGeneric(name = 'FpCalculate',
           def = function(
             mol_smiles,
             mol_names,
             type = 'pubchem',
             thread = 2,
             base_dir = '.'
           ){

             cat('Calculate molecular fingerprinting...\n\n')

             if (!('result_fps.RData' %in% list.files(file.path(base_dir,
                                                                '00_intermediate_data')))) {

               temp_fun <- function(i,
                                    smiles,
                                    type){
                 mols <- rcdk::parse.smiles(smiles[i])
                 fps <- rcdk::get.fingerprint(mols[[1]], type = type)
                 return(fps)
               }

               result_fps <- BiocParallel::bplapply(X = seq_along(mol_smiles),
                                                    FUN = temp_fun,
                                                    BPPARAM = BiocParallel::SnowParam(workers = thread,
                                                                                      progressbar = TRUE),
                                                    smiles = mol_smiles,
                                                    type = type)

               names(result_fps) <- mol_names

               dir.create(file.path(base_dir, '00_intermediate_data'), recursive = TRUE)
               save(result_fps,
                    file = file.path(base_dir, '00_intermediate_data', 'result_fps.RData'))
             } else {
               cat('Detected calculated fingerprints, and load it\n\n')
               load(file.path(base_dir, '00_intermediate_data', 'result_fps.RData'))
             }


             return(result_fps)
           }
)
