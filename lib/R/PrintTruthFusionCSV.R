args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("Two argument must be supplied.n", call.=FALSE)
}
infolder = args[1]
outfolder = args[2]
mapping_csv = args[3]

tx_list <- read.csv(mapping_csv)
load(paste0(infolder, "/GeneFusions.RData"))
simulated <- data.frame(tx5 = GeneFusions$trans1, 
                        tx3 = GeneFusions$trans2,
                        junction5 = GeneFusions$junction1, 
                        junction3 = GeneFusions$junction2,
                        gene5 = GeneFusions$gene1,
                        gene3 = GeneFusions$gene2,
                        breakpoint5 = GeneFusions$bp1,
                        breakpoint3 = GeneFusions$bp2)
simulated_hgnc <- simulated %>%
  dplyr::left_join(tx_list, by = c("tx5" = "transcript")) %>%
  dplyr::rename(gene5_hgnc = gene) %>%
  dplyr::left_join(tx_list, by = c("tx3" = "transcript")) %>%
  dplyr::rename(gene3_hgnc = gene) %>%
  dplyr::mutate(fusion = paste0(gene5_hgnc, '_', gene3_hgnc)) %>%
  dplyr::select(fusion, tx5, tx3, gene5, gene3, junction5, junction3, breakpoint5,
                breakpoint3)

write.csv(simulated_hgnc, paste0(outfolder, "simulated_fusions.csv"), quote = F, row.names = F)
