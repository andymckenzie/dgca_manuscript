
library(genefilter)
library(MASS)
library(pROC)
library(docopt)
library(bayesbio)

setwd("/Users/amckenz/Dropbox/zhang/diffcorr")

###############################
# steps to perform in the program
test = FALSE

# hardcode the parameters during a test
if(test){
	n_samples = 200
	fraction_pairs = 1
	seed = floor(runif(1, 0, 100)) #42

	ddcor_test = TRUE
	DiffCoEx_test = TRUE
	dicer_test = FALSE
}

#get in parameters from the command line when in "production" use
if(!test){
	args = commandArgs(trailingOnly = TRUE)
	print(args)

	n_samples = as.numeric(args[1])
	if(args[2] == "random"){
		seed = runif(1, 0, 100000)
		print("seed")
		print(seed)
	} else {
		seed = as.numeric(args[2])
	}
	fraction_pairs = as.numeric(args[3])
	ddcor_test = as.logical(args[4])
	DiffCoEx_test = as.logical(args[5])
	dicer_test = as.logical(args[6])
}

corr_classes = "strong"

#reproducibility
set.seed(seed)

#################################
# functions

#sets the non-diagonals of a matrix to the number x.
non_diag <- function(matrix, x){
	for(i in 1:nrow(matrix)){
	  for(j in 1:ncol(matrix)){
			if(j >= i){
				next
			} else {
				matrix[i, j] = sample(c(0, x), 1, prob = c(1-fraction_pairs, fraction_pairs))
			}
	  }
	}
	#make the upper triangle the same as the lower triangle
	matrix = makeMatSym(matrix, replaceUpper = FALSE)
	return(matrix)
}

#changes the correlation structure of a particular submatrix of a larger matrix
add_sigma_to_matrix <- function(total_matrix, row_start, nrows, x){
	extract_mat = total_matrix[row_start:(row_start + nrows),
		row_start:(row_start + nrows)]
	mat = non_diag(extract_mat, x)
	total_matrix[row_start:(row_start + nrows),
		row_start:(row_start + nrows)] = mat
	return(total_matrix)
}

#find gene pair names in the relevant submatrix
find_gene_pair_names <- function(total_dc_genes, n_genes){
	res = apply(combn(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes + 1)], m = 2), 2, paste, collapse = " ")
	return(res)
}

#updates the corrrelation matrix and adds the gene pairs that were modified to a class for subsequent extraction
update_corr_structure <- function(sigma, n_genes, variance, corr_factor,
	dcor = NULL, dcor_specific = NULL){
	total_dc_genes = sigma[["total_dc_genes"]]
	sigma[["sigma_tot"]] = add_sigma_to_matrix(sigma[["sigma_tot"]],
		row_start = total_dc_genes + 1,
		nrows = n_genes, x = variance * corr_factor)
	if(!is.null(dcor)){
		if(dcor == "goc"){
			sigma[["goc_gene_pairs"]] = c(sigma[["goc_gene_pairs"]],
				find_gene_pair_names(total_dc_genes, n_genes))
		}
		if(dcor == "loc"){
			sigma[["loc_gene_pairs"]] = c(sigma[["loc_gene_pairs"]],
				find_gene_pair_names(total_dc_genes, n_genes))
		}
	}
	if(!is.null(dcor_specific)){
		if(dcor_specific == "+/0"){
			sigma[["pos_none_gene_pairs"]] = c(sigma[["pos_none_gene_pairs"]],
				find_gene_pair_names(total_dc_genes, n_genes))
		}
		if(dcor_specific == "+/-"){
			sigma[["pos_neg_gene_pairs"]] = c(sigma[["pos_neg_gene_pairs"]],
				find_gene_pair_names(total_dc_genes, n_genes))
		}
		if(dcor_specific == "0/+"){
			sigma[["none_pos_gene_pairs"]] = c(sigma[["none_pos_gene_pairs"]],
				find_gene_pair_names(total_dc_genes, n_genes))
		}
		if(dcor_specific == "0/-"){
			sigma[["none_neg_gene_pairs"]] = c(sigma[["none_neg_gene_pairs"]],
				find_gene_pair_names(total_dc_genes, n_genes))
		}
		if(dcor_specific == "-/+"){
			sigma[["neg_pos_gene_pairs"]] = c(sigma[["neg_pos_gene_pairs"]],
				find_gene_pair_names(total_dc_genes, n_genes))
		}
		if(dcor_specific == "-/0"){
			sigma[["neg_none_gene_pairs"]] = c(sigma[["neg_none_gene_pairs"]],
				find_gene_pair_names(total_dc_genes, n_genes))
		}
	}
	sigma[["total_dc_genes"]] = total_dc_genes + n_genes
	return(sigma)
}

#for gene name generation, i.e. "aa", "ab", "ac", etc
combine <- function (vecs, number){
  combn(vecs, number, paste, collapse = "")
}

letters_unique = combine(c(letters[1:26], letters[1:26]), 2)
letters_unique = sort(unique(letters_unique))

#################################
# parameters for simulation data

n_non_expr_genes = 200
n_housekeeping_genes = 100
n_activated_genes = 300
n_tot_genes = n_non_expr_genes + n_housekeeping_genes + n_activated_genes

low_mean = 2
high_mean = 50
low_var = 50
high_var = 100
positive_correlation = 0.9
negative_correlation = 0 #if this is set too low, then the matrix will be non-positive definite

######################################
#create overall correlation matrices and classes to store them
sigma_tot_a <- structure(list(
	sigma_tot = matrix(0, nrow = n_tot_genes, ncol = n_tot_genes),
	total_dc_genes = 0
	), class = "sigma")

sigma_tot_b <- structure(list(
	sigma_tot = matrix(0, nrow = n_tot_genes, ncol = n_tot_genes),
	total_dc_genes = 0,
	loc_gene_pairs = vector(),
	goc_gene_pairs = vector(),
	pos_none_gene_pairs = vector(),
	pos_neg_gene_pairs = vector(),
	none_pos_gene_pairs = vector(),
	none_neg_gene_pairs = vector(),
	neg_pos_gene_pairs = vector(),
	neg_none_gene_pairs = vector()
	), class = "sigma")

if(corr_classes == "strong"){
	#A+, B- = loss of correlation
	n_apbm = 29
	sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_apbm,
		variance = high_var, corr_factor = positive_correlation, dcor = NULL)
	sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_apbm,
		variance = high_var, corr_factor = negative_correlation, dcor = "loc",
		dcor_specific = "+/-")

	#A-, B+ = gain of correlation
	n_ambp = 29
	sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_ambp,
		variance = high_var, corr_factor = negative_correlation, dcor = NULL)
	sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_ambp,
		variance = high_var, corr_factor = positive_correlation, dcor = "goc",
		dcor_specific = "-/+")
}

sigma_tot_a_mat = sigma_tot_a[["sigma_tot"]]
sigma_tot_b_mat = sigma_tot_b[["sigma_tot"]]

########################################
#set the means and variances of the genes

#set the variance
var_diag = c(rep(high_var, n_activated_genes),
	rep(low_var, n_housekeeping_genes),
	rep(high_var, n_non_expr_genes))

#set the mean for the high expression genes
mu = rnbinom(n = (n_activated_genes + n_housekeeping_genes), mu = high_mean, size = 0.5)

#set the mean for the low expression genes
mu = c(mu, rnbinom(n = n_non_expr_genes, mu = low_mean, size = 0.5))

#if mu is less than one, set it to one
mu[mu < 1] = 1

diag(sigma_tot_a_mat) = var_diag
diag(sigma_tot_b_mat) = var_diag

#generate the simulation matrix
sim_data_a = as.matrix(mvrnorm(n_samples, mu = mu, Sigma = sigma_tot_a_mat, tol = 0.4))
sim_data_b = as.matrix(mvrnorm(n_samples, mu = mu, Sigma = sigma_tot_b_mat, tol = 0.4))

true_cases = c(sigma_tot_b$loc_gene_pairs, sigma_tot_b$goc_gene_pairs)

####################################
#run the ddcor differential correlation pipeline on the simulated data

if(ddcor_test){

library(DGCA)

print("starting the ddcor pipeline")

sim_data_a = t(sim_data_a)
sim_data_b = t(sim_data_b)

sim_data_merge = cbind(sim_data_a, sim_data_b)
rownames(sim_data_merge) = letters_unique[1:nrow(sim_data_merge)]

conditions = c(rep("cond_a", ncol(sim_data_a)),
  rep("cond_b", ncol(sim_data_b)))
design_mat = model.matrix(~ conditions + 0)
labels = c("cond_a", "cond_b")
colnames(design_mat) = labels
npairs = (nrow(sim_data_merge)^2)/2 - nrow(sim_data_merge)

ddcor_start = Sys.time()
ddcor_res = ddcorAll(inputMat = sim_data_merge, design = design_mat, compare = labels, corr_cutoff = 0.95, adjust = "perm", corrType = "pearson", nPerm = 10, heatmapPlot = FALSE, nPairs = "all")
print(paste0("time taken to run the ddcor pipeline with ", corr_classes, n_samples, " samples..."))
ddcor_time = difftime(ddcor_start, Sys.time(), units = "secs")
print(difftime(ddcor_start, Sys.time()))

ddcor_res$combos = paste(ddcor_res$Gene1, ddcor_res$Gene2, sep = " ")
ddcor_res$true = ddcor_res$combos %in% true_cases

library(MEGENA)
library(doMC)

n.perm = 100
hub.pval = 0.05
module.pval = 1
min.size = 15
max.size = 1000
n.cores = 4
save.output = FALSE
doPar = TRUE

ddcor_res_sig_05 = ddcor_res[ddcor_res$pValDiff_adj < 0.05, ]

ddcor_res_megena = ddcor_res_sig_05[ , colnames(ddcor_res_sig_05) %in% c("Gene1", "Gene2", "zScoreDiff")]
ddcor_res_megena$zScoreDiff = abs(ddcor_res_megena$zScoreDiff)

pfn_res = MEGENA::calculate.PFN(ddcor_res_megena, doPar = doPar, num.cores = n.cores)

#normalize to follow the convention that weights are in [0, 1]
pfn_res$weight = (pfn_res$weight/max(pfn_res$weight)) * 0.999999999

g = igraph::graph.data.frame(pfn_res, directed = FALSE)

MEGENA.output = MEGENA::do.MEGENA(g, mod.pval = module.pval, hub.pval = hub.pval, remove.unsig = TRUE,
 min.size = min.size, max.size = max.size, doPar = doPar, num.cores = n.cores, n.perm = n.perm,
 save.output = save.output)

MEGENA_modules = MEGENA.output$module.outpu$modules

dgcamegena_modules =
	data.frame(Genes = unlist(MEGENA_modules), Modules = rep(names(MEGENA_modules), sapply(MEGENA_modules, length)))

if(is.null(dgcamegena_modules)){
	grey_module = data.frame(Genes = rownames(sim_data_merge), Modules = rep("grey", nrow(sim_data_merge)))
	dgcamegena_modules = grey_module
} else {
	genes_in_grey_module = rownames(sim_data_merge)[!rownames(sim_data_merge) %in% dgcamegena_modules$Genes]
	grey_module = data.frame(Genes = genes_in_grey_module, Modules = rep("grey", length(genes_in_grey_module)))
	dgcamegena_modules = rbind(dgcamegena_modules, grey_module)
}

}

###################################
#run the DiffCoEx differential correlation pipeline on the simulated data

if(DiffCoEx_test){

	print("starting the DiffCoEx pipeline")

	library(flashClust)

	beta1 = 6 #user defined parameter for soft thresholding
	AdjMatC1 = sign(cor(t(sim_data_a), method="spearman"))*(cor(t(sim_data_a), method="spearman"))^2
	AdjMatC2 = sign(cor(t(sim_data_b), method="spearman"))*(cor(t(sim_data_b), method="spearman"))^2
	diag(AdjMatC1) = 0
	diag(AdjMatC2) = 0
	collectGarbage()

	dissTOMC1C2 = TOMdist((abs(AdjMatC1-AdjMatC2)/2)^(beta1/2))
	collectGarbage()

	geneTreeC1C2 = flashClust(as.dist(dissTOMC1C2), method = "average");

	dynamicModsHybridC1C2 = cutreeDynamic(dendro = geneTreeC1C2, distM = dissTOMC1C2,method="hybrid",cutHeight=0.996,deepSplit = T, pamRespectsDendro = FALSE, minClusterSize = 15);

	dynamicColorsHybridC1C2 = labels2colors(dynamicModsHybridC1C2)

	#module merging fails if all of the colors are "grey" (i.e., there is only one module), so avoid this possibility
	if(any(!dynamicColorsHybridC1C2 == "grey")){
		#the next step merges clusters which are close (see WGCNA package documentation)
		mergedColorC1C2 = mergeCloseModules(rbind(t(sim_data_a),t(sim_data_b)),dynamicColorsHybridC1C2,cutHeight=.2)$color
	} else {
		mergedColorC1C2 = dynamicColorsHybridC1C2
	}

	#We write each module to an individual file containing affymetrix probeset IDs
	diffcoex_modules =
		data.frame(Genes = letters_unique[1:nrow(sim_data_a)], Modules = mergedColorC1C2)

}

###########################################

if(dicer_test){

	print("starting the dicer pipeline")

	library(bayesbio)

	ID = 1:ncol(sim_data_merge)
	sim_df = rbind(ID, sim_data_merge)
	write.delim(sim_df,
		paste0("dicer_sim/sim_gnxp.tsv"), row.names = TRUE, col.names = FALSE)

	sim_classes = data.frame(a = c(0, ncol(sim_data_a)), b = c(ncol(sim_data_a) + 1, (ncol(sim_data_a) + ncol(sim_data_b))), c = c("cond_a", "cond_b"))
	write.delim(sim_classes, "dicer_sim/sim_classes.tsv", col.names = FALSE)

	#dicer.jar needs to be in the same directory
	system("java -jar dicer.jar dicer_sim/sim_gnxp.tsv dicer_sim/sim_classes.tsv 0 dicer_sim/sim_res.txt")

	dicer_modules = tryCatch(read.table("dicer_sim/sim_res.txt", header = FALSE), error=function(e) NULL)
	if(!is.null(dicer_modules)){
		colnames(dicer_modules) = c("Genes", "Modules")
		dicer_line = c(seed, fraction_pairs, n_samples, nrow(dicer_modules))
	} else {
		print("no dicer modules detected")
		dicer_line = c(seed, fraction_pairs, n_samples, 0)
		genes_in_grey_module = rownames(sim_data_merge)
		grey_module = data.frame(Genes = genes_in_grey_module, Modules = rep("grey", length(genes_in_grey_module)))
		dicer_modules = rbind(grey_module)
	}

	write.table(as.matrix(t(dicer_line)), file = "dicer_sim/dicer_sim_module_res.txt",
		sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE)

}

#############################
# compare the modules detected

module_goc = unique(unlist(strsplit(sigma_tot_b$goc_gene_pairs, " ")))
module_loc = unique(unlist(strsplit(sigma_tot_b$loc_gene_pairs, " ")))

#find the modules that best overlap for each, based on jaccard index
calculate_module_jaccard <- function(modules, true_module){

	unique_module_names = unique(modules$Modules)
	jaccard_stats_true_module = vector()
	for(i in 1:length(unique_module_names)){
		jaccard_stats_true_module[i] = jaccardSets(true_module,
			modules[modules$Modules == unique_module_names[i], "Genes"])
	}
	modules_max_jaccard = max(jaccard_stats_true_module)
	return(modules_max_jaccard)

}

calculate_module_sensitivity <- function(modules, true_module){

	unique_module_names = unique(modules$Modules)
	sensitivity_stats_true_module = vector()
	for(i in 1:length(unique_module_names)){
		if(length(modules[modules$Modules == unique_module_names[i], "Genes"]) > 50){
			sensitivity_stats_true_module[i] = 0
		} else {
			sensitivity_stats_true_module[i] = length(intersect(true_module,
				modules[modules$Modules == unique_module_names[i], "Genes"]))/length(true_module)
		}
	}
	modules_max_sensitivity = max(sensitivity_stats_true_module)
	return(modules_max_sensitivity)

}

find_true_module_jaccard_sensitivity <- function(modules, true_module, data_set, module_type){

	max_jaccard = calculate_module_jaccard(modules, true_module)
	line_out = c(seed, fraction_pairs, n_samples, data_set, module_type, max_jaccard)
	write.table(as.matrix(t(line_out)), file = "module_DC_jaccard.txt",
		sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)

	max_sensitivity = calculate_module_sensitivity(modules, true_module)
	line_out = c(seed, fraction_pairs, n_samples, data_set, module_type, max_sensitivity)
	write.table(as.matrix(t(line_out)), file = "module_DC_sensitivity.txt",
		sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)


}

if(!is.null(dgcamegena_modules)){
	find_true_module_jaccard_sensitivity(dgcamegena_modules, module_goc, "DGCA_MEGENA", "GOC")
	find_true_module_jaccard_sensitivity(dgcamegena_modules, module_loc, "DGCA_MEGENA", "LOC")
}
if(!is.null(diffcoex_modules)){
	find_true_module_jaccard_sensitivity(diffcoex_modules, module_goc, "diffcoex", "GOC")
	find_true_module_jaccard_sensitivity(diffcoex_modules, module_loc, "diffcoex", "LOC")
}
if(!is.null(dicer_modules)){
	find_true_module_jaccard_sensitivity(dicer_modules, module_goc, "DICER", "GOC")
	find_true_module_jaccard_sensitivity(dicer_modules, module_loc, "DICER", "LOC")
}
