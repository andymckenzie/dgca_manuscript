
library(ggplot2)

jaccard_run = FALSE

if(jaccard_run){
  jaccard = read.table("module_DC_jaccard.txt")
} else {
  jaccard = read.table("module_DC_sensitivity.txt")
}

#eliminate duplicate seeds
jaccard = jaccard[!duplicated(paste0(jaccard$V1, jaccard$V2, jaccard$V3, jaccard$V4, jaccard$V5)), ]

colnames(jaccard) = c("seed", "fraction", "n_samples", "method", "type", "jaccard")

#make sure that only those seeds which include data points from all two of the methods are included
full_seeds = names(table(jaccard$seed))[table(jaccard$seed) >= 2]
jaccard = jaccard[jaccard$seed %in% full_seeds, ]
jaccard = jaccard[jaccard$n_samples %in% c(100, 200, 300, 400), ]
jaccard$method = gsub("DGCA_MEGENA", "DGCA/MEGENA", jaccard$method)
jaccard$method = gsub("diffcoex", "DiffCoEx", jaccard$method)
jaccard$fraction = paste0("κ = ", jaccard$fraction)

#for each n samples point, get the means and standard errors
std_error <- function(x) sd(x)/sqrt(length(x))
jaccard_ag = aggregate(jaccard ~ method+fraction+n_samples, jaccard, function(x) c(mean = mean(x), se = std_error(x)))
jaccard_ag$mean = jaccard_ag$jaccard[,1]
jaccard_ag$se = jaccard_ag$jaccard[,2]

# manual offsets to prevent overplotting
jaccard_ag[jaccard_ag$method == "DGCA/MEGENA", ]$n_samples = jaccard_ag[jaccard_ag$method == "DGCA/MEGENA", ]$n_samples - 10
jaccard_ag[jaccard_ag$method == "DiffCoEx", ]$n_samples = jaccard_ag[jaccard_ag$method == "DiffCoEx", ]$n_samples + 10
limits = aes(ymax = jaccard_ag$mean + jaccard_ag$se,
  ymin = jaccard_ag$mean - jaccard_ag$se)

jaccard_ag_plot = ggplot(data = jaccard_ag, aes(x = n_samples, y = mean, color = method)) +
  geom_point()  + #,
  geom_errorbar(data = jaccard_ag, aes(x = n_samples, colour = method,
    ymin = mean - se, ymax = mean + se, fill = method, group = method), width = 0.2) +
  theme_bw() + scale_x_continuous(breaks = c(100, 200, 300, 400)) +
  ylab("Sensitivity") + xlab("Number of Simulated Samples") + ylim(c(0, 1)) +
  scale_colour_manual(values = c("black", "blue", "red"),
    guide = guide_legend(title = "Method")) + facet_wrap( ~ fraction, nrow = 2) #+ geom_line()

if(jaccard_run){
  jaccard_ag_plot = jaccard_ag_plot + geom_hline(yintercept = 30/600, colour = "black", linetype = 2) +
    ylab("Jaccard Index")
  ggsave(jaccard_ag_plot, file = "jaccard_ag_plot.tiff", width = 7, height = 4.5)
} else {
  ggsave(jaccard_ag_plot, file = "sensitivity_ag_plot.tiff", width = 7, height = 4.5)
}
