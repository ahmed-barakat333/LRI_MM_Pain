# setting working directory 
setwd("~/Work/PhD/Projects/MM Pain-Bone/Computational Analysis/Sensory Data/Neuron/hmMA")

# load files 
HS_Sensory_Neuron_MA <- read.csv("data/HS_Sensory_Neuron_MA.csv")
colnames(HS_Sensory_Neuron_MA)[2] <-"Human_Score" 

HS_Mapped_Sensory_Neuron_MA <- read.csv("data/HS_Mapped_Sensory_Neuron_MA.csv")
colnames(HS_Mapped_Sensory_Neuron_MA)[1] <-"Name" 
colnames(HS_Mapped_Sensory_Neuron_MA)[2] <-"Mouse_Score" 

# inner merge files 
HM_Sensory_Neuron_overlapped_overlapped <- merge(HS_Sensory_Neuron_MA, HS_Mapped_Sensory_Neuron_MA, by = "Name")

# write file
write.csv(HM_Sensory_Neuron_overlapped_overlapped,file = "results/HM_Sensory_Neuron_overlapped_overlapped.csv", row.names = FALSE)

# plot distribution
par(mfrow = c(1, 2))
hist(HM_Sensory_Neuron_overlapped$Human_Score)
hist(HM_Sensory_Neuron_overlapped$Mouse_Score)

#  overlayed  bar plot (needs filtering to zoom in)
par(mfrow = c(1, 1))
require(ggplot2)
require(reshape2)

HM_Sensory_Neuron_overlapped <- HM_Sensory_Neuron_overlapped[order(HM_Sensory_Neuron_overlapped$Human_Score), ]
HM_Sensory_Neuron_overlapped_filtered <- HM_Sensory_Neuron_overlapped[301:400,]

HM_Sensory_Neuron_overlapped_filtered_m <- melt(HM_Sensory_Neuron_overlapped_filtered)

ggplot(HM_Sensory_Neuron_overlapped_filtered_m, aes(Name, value, fill = variable)) + 
  geom_bar(stat = "identity", position = "identity", alpha = 0.75) + 
  xlab("Gene Symbol") + ylab("Score") 

# line plot
HM_Sensory_Neuron_overlapped_filtered$Name <- as.numeric(as.factor(HM_Sensory_Neuron_overlapped_filtered$Name))

HM_Sensory_Neuron_overlapped_filtered_m <- melt(HM_Sensory_Neuron_overlapped_filtered ,  id.vars = 'Name', variable.name = 'series')

ggplot(HM_Sensory_Neuron_overlapped_filtered_m, aes(Name, value)) +
  geom_line(aes(colour = series)) 

# spearman rank correlation 
HM_Sensory_Neuron_overlapped$human_rank <- rank(HM_Sensory_Neuron_overlapped$Human_Score)
HM_Sensory_Neuron_overlapped$mouse_rank <- rank(HM_Sensory_Neuron_overlapped$Mouse_Score)
cor.test(HM_Sensory_Neuron_overlapped$human_rank, HM_Sensory_Neuron_overlapped$mouse_rank, method= "spearman")
