#install.packages("R.matlab")

#install.packages("aods3")

setwd("/Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Large_scale_monkey_dopamine_WM/Paper_draft/large_scale_DA_WM_draft_Neuron_submission/macaque_hierarchy/src")

library(R.matlab)
# load in data from matlab
#glm.input <- readMat('glm_input.mat', maxLength=NULL, fixNames=TRUE)
bbfit.input <- readMat('../processed_data/beta_binomial_data.mat', maxLength=NULL, fixNames=TRUE)

# convert to dataframe to make R happy
bbfit.df <- as.data.frame(bbfit.input)


# perform the logistic regression
glm1 <- glm(cbind(supragranular.vec,infragranular.vec) ~ . - 1, data = bbfit.df, family = binomial(link="logit"))

# beta-binomial version
library("aods3")
# fit beta-binomial model
fm1 <- aodml(cbind(supragranular.vec,infragranular.vec) ~ . - 1, data = bbfit.df, family = "bb")
hierVals <- coef(fm1)

writeMat('../processed_data/betaBinHierVals.mat',hier_vals = hierVals)


