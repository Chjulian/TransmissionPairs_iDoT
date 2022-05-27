#   /.)
#  /)\|
# // /    Binomial Lasso regression
# /'" "
#######################################################
#Load libraries and data
#######################################################
library(glmnet)
library(pROC)
library(ggplot2)
library(tidyverse)
library(matchmaker)
asinTransform <- function(p) asin(sqrt(p)) #arcsine square root transformation

#Load data
data <- read.csv("direction_public/data/pairdata.csv")

#Set the approach (TRUE for using the identity of the transmitter, FALSE otherwise)
pairs.identity <- FALSE 

#Set threshold for consistent/inconsistent
l1 <- l2 <-  0.5

#Create a label to identify models
mylabel<- paste0('binomial_l1_', l1, '_identity',pairs.identity)

#######################################################
#Set the model formulas to test
#######################################################
if(pairs.identity){
        epi <- c("status.at.transmission.source","risk")
        sam <-c("min.n.hap.cat","hap.n.diff","divergence.time")
        gen<- c("alignment.length", "nuc.div", "pFr1.cat")
        top <- c("topology","MinRootToTip","PD.UEH", "identities","min.inter.TBL")
        formulas <- list(
                'E'=c("response",epi),
                'S'=c("response",sam),
                'G'=c("response",gen),
                'P'=c("response",top),
                'ES'=c("response",epi,sam),
                'EG'=c("response",epi,gen),
                'EP'=c("response",epi,top),
                'SG'=c("response",sam,gen),
                'SP'=c("response",sam,top),
                'GP'=c("response",gen,top),
                'ESG'=c("response",epi,sam,gen),
                'ESP'=c("response",epi,sam,top),
                'EGP'=c("response",epi,gen,top),
                'SGP'=c("response",sam,gen,top),
                'ESGP'=c("response",epi,sam,gen,top)
        )
}else{
        epi <- c("risk")
        sam <-c("min.n.hap.cat","hap.n.diff")
        gen<- c("alignment.length", "nuc.div")
        top <- c("topology","MinRootToTip","PD.UEH", "identities","min.inter.TBL")
        formulas <- list(
                'S'=c("response",sam),
                'G'=c("response",gen),
                'P'=c("response",top),
                'ES'=c("response",epi,sam),
                'EG'=c("response",epi,gen),
                'EP'=c("response",epi,top),
                'SG'=c("response",sam,gen),
                'SP'=c("response",sam,top),
                'GP'=c("response",gen,top),
                'ESG'=c("response",epi,sam,gen),
                'ESP'=c("response",epi,sam,top),
                'EGP'=c("response",epi,gen,top),
                'SGP'=c("response",sam,gen,top),
                'ESGP'=c("response",epi,sam,gen,top))
}

#######################################################
#Loop each of the formulas
#######################################################
#Create lists to save plots and additional details of the models
my.plots <- m0.list <- m0.roc.list <- m0.matriX.list <- m0.matriY.list <- list()

for(index in 1:length(formulas)){
        mydata <- data
        
        #######################################################
        #Set response using the threshold
        #######################################################
        mydata$response <- sapply(1:nrow(mydata), function(a) ifelse(mydata$p.ancestral.ER[a] <= l1, 'inconsistent', 'consistent'))
        
        #######################################################
        #Set the covariates
        #######################################################
        mydata$divergence.time <- as.integer((abs(mydata$sampling.calc.source.mid) + abs(mydata$sampling.calc.recipient.mid)))
        mydata$min.n.hap <- sapply(1:nrow(mydata), function(a) min(c(mydata$hap.n.source[a], mydata$hap.n.rec[a])))
        mydata$min.n.hap.cat <- sapply(1:nrow(mydata), function(a) ifelse(mydata$min.n.hap[a] < 10, 'low', 'high'))
        mydata$pFr1.cat <- sapply(1:nrow(mydata), function(a) ifelse(mydata$pFr1[a] < 0.75, 'multiple', 'single'))
        mydata$pFr1 <- asinTransform(mydata$pFr1)
        mydata$hap.n.diff <- (mydata$hap.n.source - mydata$hap.n.rec)
        mydata$PD.UEH <- (mydata$PD.UEH.source - mydata$PD.UEH.recipient)
        mydata$nuc.div <- (mydata$source.nuc.div - mydata$rec.nuc.div)
        mydata$MinRootToTip <- (mydata$tree.MinRootToTip.source - mydata$tree.MinRootToTip.recipient)
        if(!pairs.identity){
                mydata$hap.n.diff <- abs(mydata$hap.n.diff)
                mydata$PD.UEH <- abs(mydata$PD.UEH)
                mydata$nuc.div <- abs(mydata$nuc.div)
                mydata$MinRootToTip <- abs(mydata$MinRootToTip)
                identidades <- rep('both', nrow(mydata))
                for(identidad in 1:nrow(mydata)){
                        if(mydata$identities[identidad]!="both") {
                                if(mydata$p.ancestral.ER[identidad]>l1) {
                                        if(mydata$identities[identidad]=="source") identidades[identidad] <- "agree" else identidades[identidad] <- "disagree"
                                } else{
                                        if(mydata$p.ancestral.ER[identidad]<=l2){
                                                if(mydata$identities[identidad]=="source") identidades[identidad] <- "disagree" else identidades[identidad] <- "agree"
                                        }
                                }
                        }
                }
                mydata$identities <- identidades
                mydata$identities <- factor(mydata$identities,levels = c('both', 'agree', 'disagree'))
        }
        
        #######################################################
        #Subset data based on formula
        #######################################################
        mydata <- mydata[,formulas[[index]]]
        
        #######################################################
        #Lasso regression
        #######################################################
        #Design a matrix for the covariates
        x <- model.matrix(response~., mydata)[,-1]
        #Set the response as 0s or 1s
        y <- ifelse(mydata$response == "consistent", 1, 0)
        #K-fold leave-one-out cross-validation for Lambda
        lambdas <- 10^seq(2, -3, by = -.1)
        m0lambda <- cv.glmnet(x, y, family = "binomial", lambda = lambdas,
                              grouped=FALSE, nfolds = nrow(mydata))
        #Fit with lambda that gives minimum mean cross-validated error
        m0 <- glmnet(x, y, family = "binomial", 
                     alpha = 1, lambda = m0lambda$lambda.min)
        coef.m0 <-coef(m0)
        
        #######################################################
        #Proceed if at least one covariate per class was included 
        #######################################################
        included.classes <- character()
        terms <- coef.m0@Dimnames[[1]][coef.m0@i+1][-1]
        if(any(terms%in%c("status.at.transmission.sourcechronic","riskmsm"))) included.classes <- c(included.classes,'E')
        if(any(terms%in%c("min.n.hap.catlow","hap.n.diff","divergence.time"))) included.classes <- c(included.classes,'S')
        if(any(terms%in%c("alignment.length", "nuc.div", "pFr1.catsingle"))) included.classes <- c(included.classes,'G')
        if(any(terms%in%c("topologyPM","topologyPP", "MinRootToTip","PD.UEH", 
                          "identitiesrecipient", "identitiessource","identitiesagree", 
                          "identitiesdisagree","min.inter.TBL"))) included.classes <- c(included.classes,'P')
        if(all(unlist(strsplit(names(formulas[index]), split = ""))%in% included.classes)){
                
                #Append information to lists
                m0.element <- list(m0)
                names(m0.element) <- names(formulas)[index]
                m0.list <- c(m0.list, m0.element)
                
                m0.x.element <- list(x)
                names(m0.x.element) <- names(formulas)[index]
                m0.matriX.list <- c(m0.matriX.list, m0.x.element)
                
                m0.y.element <- list(y)
                names(m0.y.element) <- names(formulas)[index]
                m0.matriY.list <- c(m0.matriY.list, m0.y.element)
                
                #######################################################
                #Predict using the model
                #######################################################
                probabilities <- predict(object=m0, newx = x, 
                                         type="response")
                mydata$prediction <- ifelse(probabilities > 0.5, "consistent", "inconsistent")
                predicted <- sum(mydata$prediction == mydata$response)
                predicted.nominal <- mydata$prediction[mydata$prediction == mydata$response]
                non.predicted.nominal <- mydata$prediction[mydata$prediction != mydata$response]
                
                #######################################################
                #Build a ROC-AUC with confidence intervals
                #######################################################
                roc.m0 <- roc(y, as.numeric(probabilities), ci=T)
                
                roc.m0.element <- list(roc.m0)
                names(roc.m0.element) <- names(formulas)[index]
                m0.roc.list <- c(m0.roc.list,roc.m0.element)
                
                roc.m0.coords <- pROC::coords(roc.m0, x="best", best.method = "youden",transpose = F)
                roc.m0.cis <- ci(roc.m0)
                
                #Plot ROC-AUC
                plot.roc(roc.m0, reuse.auc=TRUE, print.auc=T, ci=T, 
                         legacy.axes=T, print.thres=T, 
                         main=names(formulas)[index])        
                my.plots[[names(formulas)[index]]] <- recordPlot()
                
                #######################################################
                #Save model details to a csv file
                #######################################################
                write.table(data.frame(
                        "type" =names(formulas)[index],
                        "non0Terms" = paste(coef.m0@Dimnames[[1]][coef.m0@i+1][-1], collapse = "+"),
                        "thresholds" = paste(l1, l2, sep = ","),
                        "n.covariates" = length(coef.m0@Dimnames[[1]][coef.m0@i+1][-1]),
                        "auc" = round(as.numeric(roc.m0$auc),3),
                        "auc.min" = round(as.numeric(min(roc.m0.cis)),3),
                        "auc.max" = round(as.numeric(max(roc.m0.cis)),3),
                        "lamdba" = round(m0lambda$lambda.min,3)
                ), file = paste0("direction_public/models/", mylabel, ".csv"), 
                sep=',', row.names = F, append = T,
                col.names = !file.exists(paste0("direction_public/models/", mylabel, ".csv")))
        }
}

#######################################################
#Save ROC-AUC plots (Supplementary Figure 3) and additional data
#######################################################
saveRDS(m0.list, paste0("direction_public/models/", mylabel, "_models.rds"))
saveRDS(m0.roc.list, paste0("direction_public/models/", mylabel, "_rocs.rds"))
saveRDS(m0.matriX.list, paste0("direction_public/models/", mylabel, "_matriX.rds"))
saveRDS(m0.matriY.list, paste0("direction_public/models/", mylabel, "_matriY.rds"))
pdf(paste0("direction_public/models/", mylabel, ".pdf"), onefile=TRUE)
for (my.plot in my.plots) {
        suppressWarnings(replayPlot(my.plot))
}
graphics.off()
