#' @title MdCalculation
#' @author Zhiwei Zhou
#' @param mol_smiles a smiles structure of molecule
#' @param mol_names the name of structure
#' @export
#' @examples
#' MdCalculation(mol_smiles = "COC1=CC=C(CC(O)=O)C=C1", mol_names = "4-Methoxyphenylacetic Acid")


setGeneric(name = 'MdCalculation',
           def = function(
             mol_smiles,
             mol_names
             # is_batch=FALSE,
             # is_plot=TRUE
           ){

             # cat("Start interprite molecular smiles.\n\n")

             ### interprite molecular smiles

             mol_smiles <- mol_smiles[1]
             mol_names <- mol_names[1]
             itp_smiles <- rcdk::parse.smiles(mol_smiles)


             ### plot structure---------------------------
             # cat("\n")
             # cat("Start plot structure.\n")
             # rcdkplot <- function(molecule, width=500, height=500){
             #   par(mar=c(0,0,0,0))
             #   temp1 <- view.image.2d(molecule, width, height)
             #   plot(NA,NA,xlim=c(1,10),ylim=c(1,10),xaxt='n',yaxt='n',xlab='',ylab='') # create an empty plot
             #   rasterImage(temp1,1,1,10,10)
             # }


             # if (Is.plot==T) {
             #   plot.address <- paste(location, "Structure Plot", sep = "/")
             #   dir.create(plot.address, recursive = T)
             #   plot.name <- gsub("\\:", "_", mol.names)
             #   plot.name <- gsub("\\/", "_", plot.name)
             #   plot.name <- paste(plot.name, "png", sep = ".")
             #   plot.name <- paste(plot.address, plot.name, sep = "/")
             #   png(filename = plot.name, width = 500, height = 500)
             #   rcdkplot(itp.smiles)
             #   dev.off()
             # }


             # cat("Start calculate descriptors of structure.\n\n")

             ### calculate descriptor
             desc_names <- rcdk::get.desc.names(type = "all") #descriptors name
             descriptors <- rcdk::eval.desc(itp_smiles, desc_names)

             descriptors <- data.frame(name=mol_names, descriptors, stringsAsFactors = F)
             rownames(descriptors) <- NULL
             # cat("MD calucation completed.\n")

             return(descriptors)

           }
)
