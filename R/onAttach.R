## Welcome message when package is loaded

.onAttach <- function(libname, pkgname) {
  pkg_version <- packageVersion("REMP")
  dev = as.integer(pkg_version[1L, 2L]) %% 2L == 1L

  packageStartupMessage(rep("*", 80))
  
  packageStartupMessage("REMP version ", pkg_version, if(dev) paste0 (" (Devel)"), "\n",
                        "To access full functionality of REMP, please make sure this version is current.")
  
  citation <- paste0("\nCitation:\n",
                     "  Zheng Y, Joyce BT, Liu L, Zhang Z, Kibbe WA, Zhang W, Hou L.\n",
                     "  Prediction of genome-wide DNA methylation in repetitive elements.\n",
                     "  Nucleic Acids Research. 2017;45(15):8697-711.\n",
                     "  PubMed PMID: 28911103; PMCID: PMC5587781. http://dx.doi.org/10.1093/nar/gkx587")
  
  packageStartupMessage(citation)

  packageStartupMessage(rep("*", 80))
}
