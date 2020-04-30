
#' @title CcsAnnotate
#' @author Zhiwei Zhou
#' @description
#' @param mz experimentally measured m/z.
#' @param ccs experimentally measured CCS.
#' @param candidate_list csv file. Default: NULL
#' @param rank_type 'ccs_scoring' or 'ccs_filtering'. 'ccs_scoring': annotate unknows by calculating CCS score and integrated score; 'ccs_filtering': annotate unknows by filtering. Default: "ccs_scoring"
#' @param filter_type 'percentage' or 'a_square'.
#' @param tolerance_ccs numeric. Range: 1-100
#' @param tolerance_min numeric, unit: %. Range: 0-10
#' @param tolerance_max numeric, unit: %. Range: 0-10
#' @param weight_ccs Default: 0.3
#' @param weight_msms Default: 0.7
#' @param base_dir Default: '.'
#' @param is_output TRUE
#' @param thread Default: 1
#' @export
#' @examples
#' # CCS scoring
#' candidate_list <- system.file("extdata", "demo_annotation.csv", package="AllCCS")
#' test <- CcsAnnotate(mz = 666.312,
#'                     ccs = 261.0582,
#'                     candidate_list = candidate_list,
#'                     rank_type = 'ccs_scoring',
#'                     tolerance_min = 2,
#'                     tolerance_max = 4,
#'                     weight_ccs = 0.3,
#'                     weight_msms = 0.7,
#'                     base_dir = '.',
#'                     is_output = TRUE)
#'
#' # CCS filtering
#' test <- CcsAnnotate(mz = 666.312,
#'                     ccs = 261.0582,
#'                     candidate_list = candidate_list,
#'                     rank_type = 'ccs_filtering',
#'                     filter_type = 'percentage',
#'                     tolerance_ccs = 4,
#'                     base_dir = '.',
#'                     is_output = TRUE)


# candidate_list <- 'F:/01 MetIMMS/00 data processing/191012_allccs_annotate_function_in_package/example/L0614/demo_data/demo_example.csv'

# CCS scoring
# candidate_list <- './inst/extdata/demo_annotation.csv'
# test <- CcsAnnotate(mz = 666.312,
#                     ccs = 261.0582,
#                     candidate_list = candidate_list,
#                     rank_type = 'ccs_scoring',
#                     tolerance_min = 2,
#                     tolerance_max = 4,
#                     weight_ccs = 0.3,
#                     weight_msms = 0.7,
#                     base_dir = '.',
#                     is_output = TRUE)

# # CCS filtering
# candidate_list <- './inst/extdata/demo_annotation.csv'
# test <- CcsAnnotate(mz = 666.312,
#                     ccs = 261.0582,
#                     candidate_list = candidate_list,
#                     rank_type = 'ccs_filtering',
#                     filter_type = 'percentage',
#                     tolerance_ccs = 4,
#                     base_dir = '.',
#                     thread = 2,
#                     is_output = TRUE)



setGeneric(name = 'CcsAnnotate',
           def = function(
             mz = NULL,
             ccs = NULL,
             candidate_list = NULL,
             rank_type = c('ccs_scoring', 'ccs_filtering'),
             filter_type = c('percentage', 'a_square'),
             tolerance_ccs = 4,
             tolerance_min = 2,
             tolerance_max = 4,
             weight_ccs = 0.3,
             weight_msms = 0.7,
             base_dir = '.',
             thread = 2,
             is_output = TRUE
           ){
             rank_type <- match.arg(rank_type)
             filter_type <- match.arg(filter_type)

             if (is.null(mz)) {
               stop('Please input mz\n')
             }

             if (is.null(ccs)) {
               stop('Please input ccs\n')
             }

             if (is.null(candidate_list)) {
               stop('Please input candidate_list\n')
             }

             ccs <- as.numeric(ccs)
             tolerance_ccs <- as.numeric(tolerance_ccs)
             tolerance_min <- as.numeric(tolerance_min)
             tolerance_max <- as.numeric(tolerance_max)
             weight_ccs <- as.numeric(weight_ccs)
             weight_msms <- as.numeric(weight_msms)

             cat('Read candidate list...\n\n')
             raw_data <- readr::read_csv(candidate_list)

             temp <- match(c('rank', 'name', 'smiles', 'inchikey', 'adduct', 'score'), colnames(raw_data))
             if (any(is.na(temp))) {
               stop('Please check the name of columns\n')
             }

             raw_data <- raw_data %>%
               dplyr::select(rank, name, smiles, inchikey, adduct, score) %>%
               dplyr::mutate(rank = as.character(rank)) %>%
               dplyr::filter(adduct %in% c('[M+H]+', '[M+H-H2O]+', '[M+Na]+', '[M+NH4]+',
                                           '[M-H]-', '[M+Na-2H]-', '[M+HCOO]-'))

             temp_path <- file.path(base_dir, 'ccs_annotate')
             dir.create(temp_path, recursive = T, showWarnings = F)

             cat('Predict CCS values for candidate list...\n\n')
             pred_ccs <- try(CcsPredict(mol_smiles = raw_data$smiles,
                                        mol_names = raw_data$rank,
                                        base_dir = temp_path,
                                        thread = thread,
                                        is_output = FALSE),
                             silent = TRUE)

             if (class(pred_ccs) == 'try-error') {
               stop('CCS prediction occurs errors\n')
             }

             pred_ccs <- pred_ccs %>%
               dplyr::mutate(id = paste0(name, '_', adduct)) %>%
               dplyr::select(id, pred_ccs)

             data_add_ccs <- raw_data %>%
               dplyr::mutate(id = paste0(rank, '_', adduct)) %>%
               dplyr::left_join(y = pred_ccs, by = 'id') %>%
               dplyr::select(-id)

             cat('Candidate filtering & scoring...\n\n')
             if (rank_type == 'ccs_filtering') {
               temp <- GetCcsRange(data = ccs,
                                   type = filter_type,
                                   tolerance = tolerance_ccs) %>%
                 as.numeric()

               result <- data_add_ccs %>%
                 dplyr::rename(msms_rank = rank,
                               msms_score = score) %>%
                 dplyr::filter(pred_ccs >= temp[1] & pred_ccs <= temp[2])

               result <- result %>%
                 dplyr::mutate(rank = seq(nrow(result))) %>%
                 dplyr::select(rank, name:msms_score, dplyr::everything())

               if (is_output) {
                 readr::write_csv(result,
                                  path = file.path(temp_path,
                                                   'annotate_result_filtering.csv'))
               }

             }

             if (rank_type == 'ccs_scoring') {
               temp <- GetCcsRange(data = ccs,
                                   type = 'percentage',
                                   tolerance = tolerance_max) %>%
                 as.numeric()

               result <- data_add_ccs %>%
                 dplyr::rename(msms_rank = rank,
                               msms_score = score) %>%
                 dplyr::mutate(msms_score = msms_score/max(msms_score)) %>%
                 dplyr::filter(pred_ccs >= temp[1] & pred_ccs <= temp[2])

               errors <- abs(result$pred_ccs - ccs)/ccs*100

               result$ccs_score <- GetTrapezoidalScore(delta = errors,
                                                       pf_range = tolerance_min,
                                                       range = tolerance_max)
               result$inte_score <- result$msms_score*weight_msms + result$ccs_score*weight_ccs

               result <- result %>% dplyr::arrange(desc(inte_score))

               result <- result %>%
                 dplyr::mutate(rank = seq(nrow(result))) %>%
                 dplyr::select(rank, name:msms_score, dplyr::everything())

               if (is_output) {
                 readr::write_csv(result,
                                  path = file.path(temp_path,
                                                   'annotate_result_scoring.csv'))
               }
             }

             cat('Completed!\n')


             return(result)

           }
)


#' @title GetCcsRange
#' @author Zhiwei Zhou
#' @param data Numeric.
#' @param type 'percentage', 'a_square'. Default: percentage
#' @param tolerance Numeric. CCS deviation tolerance. Default: 4
#' @return a table of range. column 1: minimum value; column 2: maxmium value.
#' @export
#' @examples
#' GetCcsRange(data = c(100, 200),
#'             type = 'percentage',
#'             tolerance_ccs = 4)
#' GetCcsRange(data = c(100, 200),
#'             type = 'a_square',
#'             tolerance_ccs = 4)


# GetCcsRange(data = c(100, 200),
#             type = 'percentage',
#             tolerance_ccs = 4)
#
# GetCcsRange(data = c(100, 200),
#             type = 'a_square',
#             tolerance_ccs = 4)


setGeneric(name = 'GetCcsRange',
           def = function(data,
                          type = c('percentage', 'a_square'),
                          tolerance = 4){

             result <- sapply(data, function(x){
               switch (type,
                       'percentage' = {
                         x*(1+c(-1, 1)*tolerance/100)
                       },
                       'a_square' = {
                         x + (c(-1, 1) * tolerance)
                       }
               )
             })

             result <- t(result)
             return(result)
           }
)



#' @title GetTrapezoidalScore
#' @author Zhiwei Zhou
#' @description Trapezodial scoring function
#' @param delta Numeric. The difference between experiment value and library value.
#' @param pf_range Numeric. Penalty free range.
#' @param range Numeric. The maxmium tolerance.
#' @export

setGeneric('GetTrapezoidalScore',
           def = function(delta,
                          pf_range,
                          range){
             delta <- abs(delta)
             delta[delta <= pf_range] <- pf_range
             score <- 1-((delta-pf_range)/(range-pf_range))
             score[score<0] <- 0
             return(score)
           })
