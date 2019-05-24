#' @title mz_transform
#' @description calculate the m/z of adducts
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @usage
#' \code{mz_transform(M, adduct, polarity=c('pos', 'neg'))}
#' @param M a vector of exact mass.
#' @param adduct a vector of adducts. These adducts need to be
#' @param polarity 'pos' or 'neg'
#' @param formula sdf
#' @export
#' @examples
#' mz_transform(M = 123.2324, adduct = c('[M+H]+', '[M+Na]+'), polarity = 'pos')
#' mz_transform(formula="C2H5OH", adduct = c('M-', '[M-H]-'), polarity = 'neg')


# transform exact mass to mz according adduct form
mz_transform <- function(M=NULL,
                         adduct=NULL,
                         formula=NULL,
                         polarity=c("pos", 'neg')){

  polarity <- match.arg(polarity)

  data('adduct.table.rda', envir = environment())
  # load(data('data/adduct.table.rda', package = 'AllCCS'))
  # load(data('data/adduct.table.rda', package = 'AllCCS'))

  if (length(M) > 1 | length(formula)>1) {
    stop('Please only input one compound\n')
  }

  if (all(is.null(M), is.null(formula))) {
    stop('Please input M or formula.')
  }


  if (!is.null(M)) {
    M <- as.numeric(M)
  } else {
    M <- Calcu_EM(formula)
  }


  if (polarity=="pos") {

    if (is.null(adduct)) {
      adduct <- adduct.table$pos$adduct
    } else {
      temp <- adduct %in% adduct.table$pos$adduct

      if (sum(temp)!=length(adduct)) {
        stop("Please check the adduct list.")
      }
    }

    adduct.table <- adduct.table$pos

  } else {
    if (is.null(adduct)) {
      adduct <- adduct.table$neg$adduct
    } else {
      temp <- adduct %in% adduct.table$neg$adduct

      if (sum(temp)!=length(adduct)) {
        stop("Please check the adduct list.")
      }
    }

    adduct.table <- adduct.table$neg

  }

  idx.adduct <- match(adduct, adduct.table$adduct)

  mz <- sapply(c(1:length(idx.adduct)), function(i){
    temp.idx <- idx.adduct[i]
    mz <- M + adduct.table$mz[temp.idx]
  })

  output <- data.frame(M=M,
                       adduct=adduct,
                       mz=mz,
                       stringsAsFactors = F)
  return(output)
}


#' @title Calcu_EM
#' @description Calculate exact mass using formula
#' @author Zhiwei Zhou
#' \email {zhouzw@@sioc.ac.cn}
#' @param formula
#' @example
#' Calcu_EM("C2H5OH")
#' @export

Calcu_EM <- function(formula="C2H5OH") {
  molecule <- Rdisop::getMolecule(formula)
  # getFormula(molecule)
  Rdisop::getMass(molecule)
}




#' @title ZZWImpute
#' @description impute NA using KNN algorthem
#' @author Zhiwei Zhou
#' @param raw_matrix row: feature; column: sample
#' @param label feature id
#' @export

setGeneric('ZZWImpute',
           def = function(
             raw_matrix,
             label,
             k=10,
             rowmax=0.5,
             colmax=0.8,
             maxp=1500,
             rng.seed=362436069
           ){
             n_inf <- apply(raw_matrix, 2, function(x){
               sum(is.infinite(x))
             })

             n_na <- apply(raw_matrix, 2, function(x){
               temp <- is.na(x)
               sum(temp)
             })

             n_nan <- apply(raw_matrix, 2, function(x){
               temp <- is.nan(x)
               sum(temp)
             })

             cat(paste0('The number of inf: ', sum(n_inf), '\n'),
                 paste0('The number of NA: ', sum(n_na), '\n'),
                 paste0('The number of NAN: ', sum(n_nan), '\n')
             )

             result <- impute::impute.knn(as.matrix(raw_matrix),
                                          k = k,
                                          rowmax = rowmax,
                                          colmax = colmax,
                                          maxp = maxp,
                                          rng.seed = rng.seed)

             result <- result[[1]]
             result <- data.frame(name=label,
                                  result,
                                  stringsAsFactors = F)

             return(result)
           }
)


quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
