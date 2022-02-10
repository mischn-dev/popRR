library(ggplot2)
library(data.table)
library(plyr)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

input = args[1] # path to file


# read the files
rr = fread(input, header = T)

### generate some plots 

## 1. density plot for all samples 

# melt the dt
rr4 = melt(rr, id.vars = c("Group", "Chr", "pos"))
# calc mean value for each sample 
mv <- ddply(rr4, "variable", summarise, grp.mean=mean(value, na.rm = T))

p1 = ggplot(rr5, aes(x=value, fill=variable)) + geom_density(alpha=0.4) + theme_minimal() + theme(legend.position = c(0.7, 0.7)) + geom_vline(data=mv2, aes(xintercept=grp.mean, color=variable),linetype="dashed") + 
  xlab("recombination rate in window [cM/MB]") +  ylab("Count") + ggtitle("Recombination rate distribution") + labs(fill="Sample", color="Sample") 

p2 = ggplot(rr5, aes(x=pos, y=value, color = variable)) + facet_wrap(.~Chr, ncol=1) + geom_line() + geom_point(size=1.1) + labs(color = "Sample") + ylab("Recombination rate [cM/MB]") + xlab("Physical position on chromosome") +
  theme_minimal() + theme(legend.position = "none", axis.ticks.x=element_blank(), axis.text.x=element_blank()) 


png("RecombinationRate_plot.png", width = 600, height = 250, res = 300, units = "mm")
grid.arrange(p1,p2, ncol=1)
dev.off()

print(paste("Figure was create in", getwd()))