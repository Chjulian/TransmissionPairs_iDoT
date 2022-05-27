#   /.)
#  /)\|
# // /    Multinomial lasso regression
# /'" "
#######################################################
#Load libraries and data
#######################################################
library(ordinalNet)
library(multiROC)
library(dummies)
library(ggplot2)
library(tidyverse)
library(matchmaker)
asinTransform <- function(p) asin(sqrt(p)) #arcsine square root transformation

#Load data
data <- read.csv("data/pairdata.csv")

#Set the approach
#TRUE for using the identity of the transmitter, FALSE otherwise:
pairs.identity <- FALSE 
#TRUE for using ML, FALSE for using MPR:
probabilistic <- TRUE

#Set threshold for consistent, equivocal and inconsistent
l1 <- 0.4 #Either 0.4 OR 0.05
l2 <- 1-l1 

#create label to identify models
if(!probabilistic) l1 <- l2 <- 0.5
mylabel <- paste0('ordinal_l_', l1, '-', l2, '_prob',probabilistic,'_identity',pairs.identity)

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
#list to save details of the models
my.plots <- m0.list <- m0.roc.list <- m0.matriX.list <- m0.matriY.list <- list()

for(index in 1:length(formulas)){
        mydata <- data

        #######################################################
        #Set response using the thresholds
        #######################################################
        if(probabilistic) {
                mydata$response <- sapply(1:nrow(mydata), function(a) ifelse(mydata$p.ancestral.ER[a] <= l1, 'inconsistent', ifelse(mydata$p.ancestral.ER[a] >= l2, 'consistent', 'equivocal')))
        } else {
                dictionary <- data.frame('from'=c(0,0.5,1),
                                         'to'= c('inconsistent', 'equivocal','consistent'))
                mydata$response <- match_vec(mydata$p.ancestral.pars, dictionary)
                l1 <- l2 <- 0.5
        }
        
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
        mydata<- mydata[,formulas[[index]]]
        
        #######################################################
        #Lasso regression
        #######################################################
        #Design a matrix for covariates
        x <- model.matrix(response~., mydata)[,-1]
        #Set response as ordinal categorical
        y.ord <-factor(mydata$response, ordered = T,
                       levels = c("consistent","equivocal","inconsistent"))
        #K-fold leave-one-out cross-validation for lambda
        lambdas <- 10^seq(2, -3, by = -.1)
        print('Tuning lamdba')
        m0.ord.lambda <- ordinalNetTune(x, y.ord,
                                      lambdaVals = lambdas,
                                      nFolds = nrow(mydata),
                                      printProgress=F)
        bestLambdaIndex <- which.max(rowMeans(m0.ord.lambda$loglik))
        #Fit with lambda that gives higher log likelihood
        m0.ord <- ordinalNet(x, y.ord, alpha=1,
                             lambdaVals = lambdas[bestLambdaIndex])
        coef.m0.ord <-coef(m0.ord)

        #######################################################
        #proceed if at least one covariate per class was included 
        #######################################################
        included.classes <- character()
        terms <- names(coef.m0.ord)[coef.m0.ord!=0]
        terms<- terms[-grep("Intercept", terms)]
        if(any(terms%in%c("status.at.transmission.sourcechronic","riskmsm"))) included.classes <- c(included.classes,'E')
        if(any(terms%in%c("min.n.hap.catlow","hap.n.diff","divergence.time"))) included.classes <- c(included.classes,'S')
        if(any(terms%in%c("alignment.length", "nuc.div", "pFr1.catsingle"))) included.classes <- c(included.classes,'G')
        if(any(terms%in%c("topologyPM","topologyPP", "MinRootToTip","PD.UEH", 
                               "identitiesrecipient", "identitiessource","identitiesagree", 
                               "identitiesdisagree","min.inter.TBL"))) included.classes <- c(included.classes,'P')
        if(all(unlist(strsplit(names(formulas[index]), split = ""))%in% included.classes)){
                
                #Append information to lists
                m0.element <- list(m0.ord)
                names(m0.element) <- names(formulas)[index]
                m0.list <- c(m0.list, m0.element)
                
                m0.x.element <- list(x)
                names(m0.x.element) <- names(formulas)[index]
                m0.matriX.list <- c(m0.matriX.list, m0.x.element)
                
                m0.y.element <- list(y.ord)
                names(m0.y.element) <- names(formulas)[index]
                m0.matriY.list <- c(m0.matriY.list, m0.y.element)
                
                #######################################################
                #Predict using the model
                #######################################################
                true_y <- suppressWarnings(dummies::dummy(y.ord, sep = "."))
                true_y <- data.frame(true_y)
                colnames(true_y) <- gsub(".*?\\.", "", colnames(true_y))
                colnames(true_y) <- paste(colnames(true_y), "_true", sep="")

                probabilities.ord <- predict(object=m0.ord, newx = x, 
                                             type="response")
                probabilities.ord <-data.frame(probabilities.ord)
                colnames(probabilities.ord) <- levels(y.ord)
                colnames(probabilities.ord) <- paste(colnames(probabilities.ord), 
                                                     "_pred_Lasso", sep="")
                probabilities.ord <- cbind(true_y, probabilities.ord)
                
                #######################################################
                #Build a (Macro) ROC-AUC with confidence intervals
                #######################################################
                roc.ord <- suppressWarnings(multi_roc(probabilities.ord, force_diag=T))
                
                roc.m0.element <- list(roc.ord)
                names(roc.m0.element) <- names(formulas)[index]
                m0.roc.list <- c(m0.roc.list,roc.m0.element)
                
                boolFalse <- FALSE
                print('using bootstrap to generate confidence intervals for ROC-AUC')
                while(boolFalse==F) {
                        tryCatch({
                                roc.m0.cis <- suppressWarnings(roc_ci(probabilities.ord))
                                boolFalse<-TRUE
                        },error=function(e){
                        },finally={})
                }
                
                #Plot ROC-AUC
                roc.plot <- plot_roc_data(roc.ord)
                roc.plot <- roc.plot[roc.plot$Group%in%c('inconsistent','equivocal',
                                                         'consistent',
                                                         'Macro'),]
                print(ggplot(roc.plot, aes(x = 1-Specificity, y=Sensitivity)) + 
                        geom_path(aes(color = Group)) + ggtitle(names(formulas)[index]) +
                        geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour='grey', linetype = 'dotdash') + 
                        theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                                           legend.justification=c(1, 0), 
                                           legend.position=c(.95, .05), 
                                           legend.title=element_blank()) +
                        scale_color_manual(name = "", values = c("consistent" = "#76A08A",
                                                                "equivocal"= "#C4CFD0",
                                                                "inconsistent" = "#FDDDA4",
                                                                'Macro'='black')))
                my.plots[[names(formulas)[index]]] <- recordPlot()
                
                #######################################################
                #Save model details to a csv file
                #######################################################
                write.table(data.frame(
                        "type" = names(formulas)[index],
                        "non0Terms" = paste(terms, collapse = "+"),
                        "thresholds" = paste(l1, l2, sep = ","),
                        "n.covariates" = length(terms),
                        "auc" = round(as.numeric(roc.ord$AUC$Lasso$macro),3),
                        "auc.min" = round(roc.m0.cis$basic[4],3),
                        "auc.max" = round(roc.m0.cis$basic[5],3),
                        "lamdba" = round(lambdas[bestLambdaIndex],3)
                ), file = paste0("models/", mylabel, ".csv"), 
                sep=',', row.names = F, append = T,
                col.names = !file.exists(paste0("models/", mylabel, ".csv")))
                
        }

}

#######################################################
#Save ROC-AUC plots (Supplementary Figure 3) and additional data
#######################################################
saveRDS(m0.list, paste0("models/", mylabel, "_models.rds"))
saveRDS(m0.roc.list, paste0("models/", mylabel, "_rocs.rds"))
saveRDS(m0.matriX.list, paste0("models/", mylabel, "_matriX.rds"))
saveRDS(m0.matriY.list, paste0("models/", mylabel, "_matriY.rds"))
pdf(paste0("models/", mylabel, ".pdf"), onefile=TRUE)
for (my.plot in my.plots) {
        suppressWarnings(replayPlot(my.plot))
}
graphics.off()

