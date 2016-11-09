
library(ggplot2)
library(reshape2)

auc = read.table("auc_results.txt")
auc_df = data.frame(matrix(NA, nrow = length(auc$V1)/5, ncol = 5))
for(i in 1:length(auc$V1)){
  col = ifelse(!(i %% 5 == 0), i %% 5, 5)
  auc_df[ceiling(i/5), col] = auc$V1[i]
}
colnames(auc_df) = c("method", "classes", "n_samples", "seed", "time", "auc")
auc_df$time = as.numeric(auc_df$time)
auc_df$auc = as.numeric(auc_df$auc)
auc_df$n_samples = as.numeric(auc_df$n_samples)

#make sure that only those seeds which include data points from all four of the methods are included
full_seeds = names(table(auc_df$seed))[table(auc_df$seed) >= 4]
auc_full = auc_df[auc_df$seed %in% full_seeds, ]
auc_full$method = gsub("ddcor", "DGCA", auc_full$method)
auc_full$method = gsub("discordant", "Discordant", auc_full$method)
auc_full$method = gsub("Discordant_t", "Discordant (FT)", auc_full$method)
auc_full$method = gsub("ebcoexpress", "EBcoexpress", auc_full$method)
auc_full$time = abs(auc_full$time)

#for each n samples point, get the means and standard errors
std_error <- function(x) sd(x)/sqrt(length(x))
auc_ag = aggregate(auc ~ n_samples+method, auc_full, function(x) c(mean = mean(x), se = std_error(x)))
auc_ag$mean = auc_ag$auc[,1]
auc_ag$se = auc_ag$auc[,2]
limits = aes(ymax = mean + se, ymin = mean - se)

auc_ag_plot = ggplot(data = auc_ag, aes(x = n_samples, y = mean, colour = method)) + geom_point() +
  geom_errorbar(limits, width = 0.2) + theme_bw() + scale_x_continuous(breaks = c(10, 30, 50, 70, 90, 100)) +
  ylab("AUC") + xlab("Number of Simulated Samples") + ylim(c(0.5, 1)) +
  scale_colour_manual(values = c("black", "red", "orange", "blue"),
    guide = guide_legend(title = "Method")) + geom_line()
