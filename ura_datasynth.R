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
data1 <- read.csv("output_mood10.csv")[-1]

for(i in seq(11,20000,1)){
  
  try(df <- read.csv(paste("output_mood",i,".csv", sep=""))[-1])
  try(data1 <- rbind(data1,df))
  
  data1 <- unique(data1)
  
  try(file.remove(paste("output_mood",i,".csv", sep="")))
}

write.csv(data1,"final_output_mood.csv")
print("saved file successfully")

minval <- with(data1, data1[sse == min(sse),])
write.csv(minval, "minval_mood_apr5.csv")
###################