#' @title MzTransform
#' @description calculate the m/z of adducts
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' \code{mz_transform(M, adduct, polarity=c('pos', 'neg'))}
#' @param M a vector of exact mass.
#' @param adduct a vector of adducts.
#' @param polarity 'pos' or 'neg'
#' @param formula chemical formula
#' @export
#' @examples
#' mz_transform(M = 123.2324, adduct = c('[M+H]+', '[M+Na]+'), polarity = 'pos')
#' mz_transform(formula="C2H5OH", adduct = c('M-', '[M-H]-'), polarity = 'neg')

# MzTransform(M = 123.2324, adduct = c('[M+H]+', '[M+Na]+'), polarity = 'pos')

MzTransform <- function(M=NULL,
                        adduct=NULL,
                        formula=NULL,
                        polarity=c("pos", 'neg')){

  polarity <- match.arg(polarity)

  if (length(M) > 1 | length(formula) >1) {
    stop('Please only input one compound\n')
  }

  if (all(is.null(M), is.null(formula))) {
    stop('Please input M or formula.')
  }


  if (!is.null(M)) {
    M <- as.numeric(M)
  } else {
    M <- ExactMassCalculate(formula)
  }

  adduct_table <- dplyr::bind_rows(adduct_table$pos, adduct_table$neg)

  if (length(adduct)<1) {
    adduct <- switch (polarity,
                      'pos' = {adduct_table$pos$adduct},
                      'neg' = {adduct_table$neg$adduct}
    )
  }

  idx_adduct <- match(adduct, adduct_table$adduct)

  mz <- sapply(c(1:length(idx_adduct)), function(i){
    temp_idx <- idx_adduct[i]

    if (is.na(temp_idx)) {
      mz <- NA
      return(mz)
    }

    if (stringr::str_detect(adduct_table$adduct[temp_idx], pattern = '2M')) {
      mz <- M*2 + adduct_table$mz[temp_idx]
    } else {
      mz <- M + adduct_table$mz[temp_idx]
    }

    mz

  })

  output <- data.frame(M=M,
                       adduct=adduct,
                       mz=mz,
                       stringsAsFactors = F)
  # cat('not here')
  return(output)
}


#' @title ExactMassCalculate
#' @name ExactMassCalculate
#' @description Calculate exact mass using formula
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param formula chemical formula
#' @examples
#' ExactMassCalculate("C2H5OH")

setGeneric(name = 'ExactMassCalculate',
           def = function(formula="C2H5OH"){
             molecule <- Rdisop::getMolecule(formula)
             # getFormula(molecule)
             Rdisop::getMass(molecule)
           }
)

# ExactMassCalculate <- function(formula="C2H5OH") {
#   molecule <- Rdisop::getMolecule(formula)
#   # getFormula(molecule)
#   Rdisop::getMass(molecule)
# }




#' @title KnnImpute
#' @description impute NA using KNN algorthem
#' @author Zhiwei Zhou
#' @param raw_matrix row: feature; column: sample
#' @param label feature id


setGeneric('KnnImpute',
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


#' @title Smiles2Mass
#' @author Zhiwei Zhou
#' @param mol_smiles See \code{\link{get.mol2formula}}
#' @param mol_names
#' @return a vector of mass
#' @export
#' @examples
#' test_smiles <- c("COC1=CC=C(CC(O)=O)C=C1", "CN(C)CCC1=CC=C(O)C=C1", "OC(=O)CCCCCCCC\\C=C\\C(O)=O")
#' test_name <- seq(length(test))
#' Smiles2Mass(mol_smiles = test, mol_names = test_name)

# test_smiles <- c("COC1=CC=C(CC(O)=O)C=C1", "CN(C)CCC1=CC=C(O)C=C1", "OC(=O)CCCCCCCC\\C=C\\C(O)=O")
# test_name <- seq(length(test))
# Smiles2Mass(mol_smiles = test, mol_names = test_name)

setGeneric(name = 'Smiles2Mass',
           def = function(
             mol_smiles,
             mol_names
           ){
             result <- pbapply::pbsapply(mol_smiles, function(x){
               temp <- rcdk::parse.smiles(x)[[1]]
               temp <- rcdk::get.mol2formula(molecule = temp)
               temp@mass
             })
             names(result) <- mol_names
             return(result)
           }
)




#' @title SmilesCheck
#' @author Zhiwei Zhou
#' @param mol_smiles See \code{\link{is.smiles}}
#' @return a vector of status, valid or error 1

setGeneric(name = 'SmilesCheck',
           function(mol_smiles) {
             status_label <- rep('Valid', length(mol_smiles))
             is_smiles <- pbapply::pbsapply(mol_smiles, function(x){
               webchem::is.smiles(x)
             })
             status_label[which(!is_smiles)] <- 'Error_1'
             return(status_label)
           }
)



#' @title IntegrityCheck
#' @name IntegrityCheck
#' @description Check hte integrity of input data
#' @author Zhiwei Zhou
#' @param mol_smiles
#' @param mol_names

setGeneric(name = 'IntegrityCheck',
           def = function(mol_smiles,
                          mol_names){

             if (length(mol_smiles)!=length(mol_names)) {
               stop('Please check the number of mol_smiles and mol_names.\n')
             }

             if (length(unique(mol_names))!=length(mol_names)) {
               stop('Please check the inputed names (Replicates existed).\n')
             }

           }
)


quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
