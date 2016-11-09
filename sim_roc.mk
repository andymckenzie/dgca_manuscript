#need to be in the directory of sim_roc.R for this script to work
#run it via make -f sim_roc.mk

.PHONY: all

all: ddcor_vs_ebcoexpress_vs_discordant_100_samples1.tiff \
	ddcor_vs_ebcoexpress_vs_discordant_100_samples2.tiff \
	ddcor_vs_ebcoexpress_vs_discordant_100_samples3.tiff \
	ddcor_vs_ebcoexpress_vs_discordant_100_samples4.tiff \
	ddcor_vs_ebcoexpress_vs_discordant_100_samples5.tiff \

ddcor_vs_ebcoexpress_vs_discordant_100_samples1.tiff: sim_roc.R
	Rscript sim_roc.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "FALSE" "FALSE" "FALSE"
ddcor_vs_ebcoexpress_vs_discordant_100_samples2.tiff: sim_roc.R
	Rscript sim_roc.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "FALSE" "FALSE" "FALSE"
ddcor_vs_ebcoexpress_vs_discordant_100_samples3.tiff: sim_roc.R
	Rscript sim_roc.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "FALSE" "FALSE" "FALSE"
ddcor_vs_ebcoexpress_vs_discordant_100_samples4.tiff: sim_roc.R
	Rscript sim_roc.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "FALSE" "FALSE" "FALSE"
ddcor_vs_ebcoexpress_vs_discordant_100_samples5.tiff: sim_roc.R
	Rscript sim_roc.R "100" "random" "TRUE" "TRUE" "TRUE" "TRUE" "FALSE" "FALSE" "FALSE"
