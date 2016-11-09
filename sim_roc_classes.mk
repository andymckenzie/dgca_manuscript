#need to be in the directory of sim_roc_classes_classes.R for this script to work
#run it via make -f sim_roc_classes_classes.mk
MAKEFLAGS += -j 6  # default number of parallel jobs

.PHONY: all

all: ddcor_vs_ebcoexpress_vs_discordant_100_samples1.tiff \
	ddcor_vs_ebcoexpress_vs_discordant_100_samples2.tiff \
	ddcor_vs_ebcoexpress_vs_discordant_100_samples3.tiff \
	ddcor_vs_ebcoexpress_vs_discordant_100_samples4.tiff \
	ddcor_vs_ebcoexpress_vs_discordant_100_samples5.tiff

ddcor_vs_ebcoexpress_vs_discordant_100_samples1.tiff: sim_roc_classes.R
	Rscript sim_roc_classes.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "medium" "FALSE" "FALSE" "FALSE"
	Rscript sim_roc_classes.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "strong" "FALSE" "FALSE" "FALSE"
ddcor_vs_ebcoexpress_vs_discordant_100_samples2.tiff: sim_roc_classes.R
	Rscript sim_roc_classes.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "medium" "FALSE" "FALSE" "FALSE"
	Rscript sim_roc_classes.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "strong" "FALSE" "FALSE" "FALSE"
ddcor_vs_ebcoexpress_vs_discordant_100_samples3.tiff: sim_roc_classes.R
	Rscript sim_roc_classes.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "medium" "FALSE" "FALSE" "FALSE"
	Rscript sim_roc_classes.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "strong" "FALSE" "FALSE" "FALSE"
ddcor_vs_ebcoexpress_vs_discordant_100_samples4.tiff: sim_roc_classes.R
	Rscript sim_roc_classes.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "medium" "FALSE" "FALSE" "FALSE"
	Rscript sim_roc_classes.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "strong" "FALSE" "FALSE" "FALSE"
ddcor_vs_ebcoexpress_vs_discordant_100_samples5.tiff: sim_roc_classes.R
	Rscript sim_roc_classes.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "medium" "FALSE" "FALSE" "FALSE"
	Rscript sim_roc_classes.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "strong" "FALSE" "FALSE" "FALSE"
