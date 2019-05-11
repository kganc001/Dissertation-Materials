###############
#CONTORL PANEL#
###############
rm(list=ls())

wd <- getwd()
setwd(wd)
rm(wd)
###############


###################
#DATA SYNTHESIZING#
###################
data1 <- read.csv("simulation_output1.csv")[-1]

for(i in seq(4135,10000,1)){
  
  try(df <- read.csv(paste("simulation_output",i,".csv", sep=""))[-1])
  try(data1 <- rbind(data1,df))
  
  data1 <- unique(data1)
  
  try(file.remove(paste("simulation_output",i,".csv", sep="")))
}

write.csv(data1,"simulation_final.csv")
#print("saved file successfully")

minval <- with(data1, data1[sse == min(sse),])
write.csv(minval, "simualtion_minval.csv")
###################