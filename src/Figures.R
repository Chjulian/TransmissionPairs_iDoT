#   /.)
#  /)\|
# // /    Figures
# /'" "
#######################################################
#Load libraries and data
#######################################################
library(ggplot2)
library(ggpubr)
library(corrplot)
library(RColorBrewer)
library(tidyverse)
library(broom)
library(glmnet)
library(pROC)
library(ordinalNet)
library(matchmaker)


#######################################################
#Figure 2. The probability for each transmission pair, i, 
#that the transmitting partner is correctly identified 
#using ML ancestral state reconstruction.
#######################################################
data <- read.csv("direction_public/data/pairdata.csv")

ggplot(data, aes(x=p.ancestral.ER, fill=topology)) + 
        geom_vline(xintercept = c(0.05,0.4,0.6,0.95), linetype = "dashed", alpha=0.3) +
        geom_histogram() +
        theme_minimal() + xlab("pi") + ylab("Frequency") + theme(legend.position="bottom") +
        guides(fill=guide_legend(title="")) + theme(text=element_text(size=14)) +
        scale_fill_manual(values=c("#003c67", "#0073c2", "#efc001"), 
                          labels = c("monophyletic-monophyletic", 
                                     "paraphyletic-monophyletic", 
                                     "paraphyletic-polyphyletic")) +
        xlab("Pi (Probability of the transmitter being\ninferred as the character state at the root)") +
        theme(axis.text.x = element_text(angle = 45))

#######################################################
#Figure 3. Model results. (Panels A and C)
#######################################################
auc.data_T <- read.csv("direction_public/models/binomial_l1_0.5_identityTRUE.csv")
auc.data_T$type <- factor(auc.data_T$type, levels = c('E','S','G','P','ES','EG','EP','SG','SP','GP','ESG','ESP','EGP','SGP','ESGP'))
auc.data_T$auc.equivalent <- (max(auc.data_T$auc)-auc.data_T$auc)<0.05


auc.data_F <- read.csv("direction_public/models/binomial_l1_0.5_identityFALSE.csv")
auc.data_F$type <- factor(auc.data_F$type, levels = c('E','S','G','P','ES','EG','EP','SG','SP','GP','ESG','ESP','EGP','SGP','ESGP'))
auc.data_F$auc.equivalent <- (max(auc.data_F$auc)-auc.data_F$auc)<0.05

pA <- ggplot(auc.data_T,aes(x=auc, y=type, xmin=auc.min, xmax=auc.max)) +
        geom_linerange(position = position_dodge(0.5)) +
        geom_point(aes(size=n.covariates, color=auc.equivalent),position = position_dodge(0.5)) + 
        # geom_label(aes(label=type))+
        theme_minimal() + xlim(c(0.5,1)) + ggtitle('A') +
        scale_color_manual(values=c('TRUE'="#41ab5d",
                                    'FALSE'="#000000")) + coord_flip() +
        xlab("") + ylab("AUC") + guides(colour = "none") +
        scale_size_continuous(range = c(1, 5), breaks=c(2,4,6,8)) +
        theme(axis.text.x = element_text(angle = 45))
pC <- ggplot(auc.data_F,aes(x=auc, y=type, xmin=auc.min, xmax=auc.max)) +
        geom_linerange(position = position_dodge(0.5)) +
        geom_point(aes(size=n.covariates, color=auc.equivalent),position = position_dodge(0.5)) + 
        # geom_label(aes(label=type))+
        theme_minimal() + xlim(c(0.5,1)) + ggtitle('C') +
        scale_color_manual(values=c('TRUE'="#41ab5d",
                                    'FALSE'="#000000")) + coord_flip() +
        xlab("") + ylab("AUC") + guides(colour = "none") +
        scale_size_continuous(range = c(1, 5), breaks=c(2,4,6,8)) +
        theme(axis.text.x = element_text(angle = 45))
ggarrange(pA,pC, ncol = 1)

#######################################################
#Figure 4. Predicting the success of inferring the
#direction of transmission. (Panels A and B) The binary
#‘P’ model with three predictor covariates. 
#######################################################
#specify how many points to use for continuous covariates
length.out <- 30

#Load data of models
data <- read.csv("direction_public/data/pairdata.csv")
auc.data <- read.csv("direction_public/models/binomial_l1_0.5_identityFALSE.csv")
m0s <- readRDS("direction_public/models/binomial_l1_0.5_identityFALSE_models.rds")
roc.m0s <- readRDS("direction_public/models/binomial_l1_0.5_identityFALSE_rocs.rds")
matriX <- readRDS("direction_public/models/binomial_l1_0.5_identityFALSE_matriX.rds")
matriY <- readRDS("direction_public/models/binomial_l1_0.5_identityFALSE_matriY.rds")

#Pick the best-fit model 
auc.index <- auc.data[(max(auc.data$auc)-auc.data$auc)<0.05,] 
if(nrow(auc.index)>1){
        auc.index <- auc.index[which(auc.index$n.covariates==min(auc.index$n.covariates)),]
        if(nrow(auc.index)>1){
                auc.index <- auc.index[which(auc.index$auc==max(auc.index$auc)),]
        }
}
model.name <- auc.index$type
m0 <- m0s[[model.name]]
coef.m0 <- coef(m0)

#Set the terms of the model
terms <- coef.m0@Dimnames[[1]][coef.m0@i+1][-1]
terms.list <- list()
if(any(terms=="riskmsm")) terms.list <- c(terms.list, list(risk=c('het','msm')))
if(any(terms=="min.n.hap.catlow")) terms.list <- c(terms.list, list(min.n.hap.cat=c('high','low')))
if(any(terms=="hap.n.diff")) terms.list <- c(terms.list, list(hap.n.diff=seq(0,50,length.out=length.out)))
if(any(terms=="alignment.length")) terms.list <- c(terms.list, list(alignment.length=seq(200,9100,length.out=length.out)))
if(any(terms=="nuc.div")) terms.list <- c(terms.list, list(nuc.div=seq(4e-05,0.12,length.out=length.out)))
if(any(terms=="pFr1.catsingle")) terms.list <- c(terms.list, list(pFr1.cat=c('multiple','single')))
if(any(terms%in%c("topologyPM","topologyPP"))) {
        topology<- character()
        if(any(terms==("topologyPM"))) topology <- c(topology, 'PM')
        if(any(terms==("topologyPP"))) topology <- c(topology, 'PP')
        terms.list <- c(terms.list, list(topology=c('MM',topology)))
        rm(topology)
}
if(any(terms=="MinRootToTip")) terms.list <- c(terms.list, list(MinRootToTip=seq(1e-09,0.13,length.out=length.out)))
if(any(terms=="PD.UEH")) terms.list <- c(terms.list, list(PD.UEH=seq(5e-05,0.52,length.out=length.out)))
if(any(terms%in%c("identitiesagree","identitiesdisagree"))) {
        identities<- character()
        if(any(terms==("identitiesagree"))) identities <- c(identities, 'agree')
        if(any(terms==("identitiesdisagree"))) identities <- c(identities, 'disagree')
        identities <- c('both',identities)
        identities <- factor(identities, levels = c('both','agree','disagree'))
        terms.list <- c(terms.list, list(identities=identities))
        rm(identities)
}
if(any(terms=="min.inter.TBL")) terms.list <- c(terms.list, list(min.inter.TBL=seq(2e-09,0.11,length.out=length.out)))

#Set the model for prediction (fit w/out 0-shrinked covariates)
x <- matriX[[model.name]]; y <- matriY[[model.name]]
x <- x[,terms]
m0 <- glmnet(x, y, family = "binomial", 
             alpha = 1, lambda = m0$lambda)

#Prepare the input data
x2 <- model.data <- cross_df(terms.list)
x2 <- model.matrix(~., x2)[,-1]

#Make prediction
model.data <- as.data.frame(model.data)
model.data$probabilities <- as.vector(predict(m0, newx = x2, type="response"))

#Config the Plot
toplot <- model.data
model.data[1:(ncol(model.data)-1)] <- lapply(model.data[1:(ncol(model.data)-1)],as.factor)
model.data_l  <- pivot_longer(model.data,cols=names(model.data)[-length(names(model.data))],
                              names_to='covariate')
values <- c('het','msm',
            '200','9100',
            '4e-05','0.12',
            'multiple','single',
            'low','high',
            '0','50',
            'MM','PM','PP',
            'agree','both','disagree',
            '1e-09','0.13',
            '5e-05','0.52',
            '2e-09','0.11')
dictionary = data.frame('from'=c('Sexual risk exposure group',				
                                 'Sequence alignment length',					
                                 'Intra-host nucleotide diversity difference',
                                 'Multiplicity of infection',					
                                 'Sample size',								
                                 'Sample size difference',					
                                 'Topology class',							
                                 'Identity of the most basal tip',			
                                 'Root-to-tip difference',					
                                 'Phylogenetic diversity difference',			
                                 'Inter-host patristic distance'),
                        'to'=c('#ffffff', '#800026',
                               '#d73027', '#f46d43',
                               '#fdae61', '#fee090',
                               '#053061', '#4393c3',
                               '#92c5de', '#d7eef4',
                               '#edf3f4'))
#filter continuous covariates
model.data_l <- model.data_l[model.data_l$value%in%values,]
model.data_l$value <- factor(model.data_l$value, levels = values)
model.data_l$covariate <- factor(model.data_l$covariate, 
                                 levels=c("risk", 
                                          "alignment.length","nuc.div","pFr1.cat",
                                          "min.n.hap.cat","hap.n.diff",
                                          "topology","identities","MinRootToTip","PD.UEH","min.inter.TBL"), 
                                 labels=dictionary$from)

mypalette <- setNames(as.character(match_vec(as.character(unique(model.data_l$covariate)), dictionary)),
                      as.character(unique(model.data_l$covariate)))

#Plot
pA <- ggplot(model.data_l,aes(x=probabilities, y=value))+
        geom_violin(aes(fill=covariate))+
        stat_summary(fun = "mean",geom = "point", colour = "white") +
        theme_minimal() + theme(legend.position="bottom") +
        xlim(c(0,1)) +  
        scale_y_discrete(limits=rev) +
        scale_fill_manual(values = mypalette) +
        labs(x="Probability of correctly identifying the transmitter", 
             fill="Covariates", y="") + ggtitle(model.name) +
        facet_grid(covariate~.,scales = "free_y") + theme(strip.text.y = element_blank(),
                                                          panel.spacing = unit(2, "lines"))
pB <- ggplot()+
        geom_smooth(data=toplot[toplot$topology=='MM',],
                    aes(x=PD.UEH , y=probabilities, color=identities))+
        geom_smooth(data=toplot[toplot$topology=='PM',],
                    aes(x=PD.UEH , y=probabilities, color=identities))+
        geom_point(data=toplot[toplot$topology=='MM',]
                   ,aes(x=PD.UEH , y=probabilities, color=identities, shape=topology), size=3)+
        geom_point(data=toplot[toplot$topology=='PM',]
                   ,aes(x=PD.UEH , y=probabilities, color=identities, shape=topology), size=3)+
        theme_minimal() + ylim(0,1.0) + theme(legend.position="bottom") +
        scale_color_manual(values=c(both='#C4CFD0','agree'='#1C366B','disagree'='#F24D29'))
ggarrange(pA,pB)

#######################################################
#Figure 4. Predicting the success of inferring the
#direction of transmission. (panels C and D) The ordinal
#‘SP’ model (with a relaxed threshold for the direction
#of transmission classification t=0.6) has four predictor covariates. 
#######################################################
#Specify the approach
#TRUE for using the identity of the transmitter, FALSE otherwise:
pairs.identity <- FALSE 
#TRUE for using ML, FALSE for using MPR:
probabilistic <- TRUE
#Set threshold for consistent, equivocal and inconsistent
l1 <- 0.4 #Either 0.4 OR 0.05
l2 <- 1-l1 
#create label to load models
if(!probabilistic) l1 <- l2 <- 0.5
mylabel <- paste0('ordinal_l_', l1, '-', l2, '_prob',probabilistic,'_identity',pairs.identity)

#specify how many points to use for continuous covariates
length.out <- 20

#Load data of models
data <- read.csv("direction_public/data/pairdata.csv")
auc.data <- read.csv(paste0("direction_public/models/",mylabel,".csv"))
m0s <- readRDS(paste0("direction_public/models/",mylabel,"_models.rds"))
roc.m0s <- readRDS(paste0("direction_public/models/",mylabel,"_rocs.rds"))
matriX <- readRDS(paste0("direction_public/models/",mylabel,"_matriX.rds"))
matriY <- readRDS(paste0("direction_public/models/",mylabel,"_matriY.rds"))

#Pick the best-fit model 
auc.index <- auc.data[(max(auc.data$auc)-auc.data$auc)<0.05,] #diff 0.05
if(nrow(auc.index)>1){
        auc.index <- auc.index[which(auc.index$n.covariates==min(auc.index$n.covariates)),]
        if(nrow(auc.index)>1){
                auc.index <- auc.index[which(auc.index$auc==max(auc.index$auc)),]
        }
}
model.name <- auc.index$type
m0.ord <- m0s[[model.name]]
coef.m0 <- coef(m0.ord)

#Set the terms of the model
terms <- coef(m0.ord)[coef(m0.ord)!=0] 
terms <- names(terms[-grep("Intercept", names(terms))])
terms.list <- list()
if(any(terms=="riskmsm")) terms.list <- c(terms.list, list(risk=c('het','msm')))
if(any(terms=="min.n.hap.catlow")) terms.list <- c(terms.list, list(min.n.hap.cat=c('high','low')))
if(any(terms=="hap.n.diff")) terms.list <- c(terms.list, list(hap.n.diff=seq(0,50,length.out=length.out)))
if(any(terms=="alignment.length")) terms.list <- c(terms.list, list(alignment.length=seq(200,9100,length.out=length.out)))
if(any(terms=="nuc.div")) terms.list <- c(terms.list, list(nuc.div=seq(0.00004,0.12,length.out=length.out)))
if(any(terms=="pFr1.catsingle")) terms.list <- c(terms.list, list(pFr1.cat=c('multiple','single')))
if(any(terms%in%c("topologyPM","topologyPP"))) {
        topology<- character()
        if(any(terms==("topologyPM"))) topology <- c(topology, 'PM')
        if(any(terms==("topologyPP"))) topology <- c(topology, 'PP')
        terms.list <- c(terms.list, list(topology=c('MM',topology)))
        rm(topology)
}
if(any(terms=="MinRootToTip")) terms.list <- c(terms.list, list(MinRootToTip=seq(1e-09,0.13,length.out=length.out)))
if(any(terms=="PD.UEH")) terms.list <- c(terms.list, list(PD.UEH=seq(0.00005,0.52,length.out=length.out)))
if(any(terms%in%c("identitiesagree","identitiesdisagree"))) {
        identities<- character()
        if(any(terms==("identitiesagree"))) identities <- c(identities, 'agree')
        if(any(terms==("identitiesdisagree"))) identities <- c(identities, 'disagree')
        identities <- c('both',identities)
        identities <- factor(identities, levels = c('both','agree','disagree'))
        identities <- droplevels(identities)
        terms.list <- c(terms.list, list(identities=identities))
        rm(identities)
}
if(any(terms=="min.inter.TBL")) terms.list <- c(terms.list, list(min.inter.TBL=seq(2e-09,0.11,length.out=length.out)))

#Set the model for prediction (fit w/out 0-shrinked covariates)
x <- matriX[[model.name]]; y.ord <- matriY[[model.name]]
x <- x[,terms]
m0.ord <- ordinalNet(x[,terms], y.ord,alpha=1,
                     lambdaVals = m0.ord$lambdaVals)

#Prepare the input data
x2 <- model.data <- cross_df(terms.list)
x2 <- model.matrix(~., x2)[,-1]

#Make prediction
model.data <- as.data.frame(model.data)
probabilities.ord <- predict(object=m0.ord, newx = x2, type="response")
probabilities.ord<-data.frame(probabilities.ord)
colnames(probabilities.ord) <- levels(y.ord)
model.data <- cbind(model.data,probabilities.ord); rm(probabilities.ord)
model.data$equivocal <- NULL

#Config the Plot
toplot <- model.data
model.data[1:(ncol(model.data)-2)] <- lapply(model.data[1:(ncol(model.data)-2)],as.factor)
model.data_l  <- pivot_longer(model.data,cols=names(model.data)[1:(ncol(model.data)-2)],
                              names_to='covariate')
values <- c('het','msm',
            '200','9100',
            '4e-05','0.12',
            'multiple','single',
            'low','high',
            '0','50',
            'MM','PM','PP',
            'agree','both','disagree',
            '1e-09','0.13',
            '5e-05','0.52',
            '2e-09','0.11')
dictionary = data.frame('from'=c('Sexual risk exposure group',				
                                 'Sequence alignment length',					
                                 'Intra-host nucleotide diversity difference',
                                 'Multiplicity of infection',					
                                 'Sample size',								
                                 'Sample size difference',					
                                 'Topology class',							
                                 'Identity of the most basal tip',			
                                 'Root-to-tip difference',					
                                 'Phylogenetic diversity difference',			
                                 'Inter-host patristic distance'),
                        'to'=c('#ffffff', '#800026',
                               '#d73027', '#f46d43',
                               '#fdae61', '#fee090',
                               '#053061', '#4393c3',
                               '#92c5de', '#d7eef4',
                               '#edf3f4'))

#filter continuous covariates
model.data_l <- model.data_l[model.data_l$value%in%values,]
model.data_l$value <- factor(model.data_l$value, levels = values)
model.data_l$covariate <- factor(model.data_l$covariate, 
                                 levels=c("risk", 
                                          "alignment.length","nuc.div","pFr1.cat",
                                          "min.n.hap.cat","hap.n.diff",
                                          "topology","identities","MinRootToTip","PD.UEH","min.inter.TBL"), 
                                 labels=dictionary$from)

mypalette <- setNames(as.character(matchmaker::match_vec(as.character(unique(model.data_l$covariate)), dictionary)),
                      as.character(unique(model.data_l$covariate)))

#Plot
pc <- ggplot(model.data_l,aes(x=consistent, y=value))+
        geom_violin(aes(fill=covariate))+
        # geom_jitter(alpha=0.3, size=0.5) +
        stat_summary(fun = "mean",geom = "point", colour = "white") +
        theme_minimal() + xlim(c(0,1)) +  
        scale_y_discrete(limits=rev) +
        scale_fill_manual(values = mypalette)+
        labs(x="Probability of correctly\nidentifying the transmitter", 
             fill="Covariates", y="") +
        facet_grid(covariate~.,scales = "free_y") + theme(strip.text.y = element_blank(),
                                                          panel.spacing = unit(2, "lines"))
pi <- ggplot(model.data_l,aes(x=inconsistent, y=value))+
        geom_violin(aes(fill=covariate))+
        # geom_jitter(alpha=0.3, size=0.5) +
        stat_summary(fun = "mean",geom = "point", colour = "white") +
        theme_minimal() + xlim(c(0,1)) +  
        scale_y_discrete(limits=rev, position = "right") +
        scale_fill_manual(values = mypalette)+
        labs(x="Probability of incorrectly\nidentifying the transmitter", 
             fill="Covariates", y="") +
        facet_grid(covariate~.,scales = "free_y") + theme(strip.text.y = element_blank(),
                                                          panel.spacing = unit(2, "lines"))
pC <- ggarrange(pc, pi,common.legend = T, legend='bottom')

toplot$min.n.hap.cat<- factor(toplot$min.n.hap.cat, levels =c('low','high'))
pD <- ggplot()+
        geom_smooth(data=toplot[toplot$topology=='MM'&toplot$min.n.hap.cat=='high',],
                    aes(x=hap.n.diff , y=consistent, color=identities))+
        geom_smooth(data=toplot[toplot$topology=='PM'&toplot$min.n.hap.cat=='high',],
                    aes(x=hap.n.diff , y=consistent, color=identities))+
        geom_smooth(data=toplot[toplot$topology=='MM'&toplot$min.n.hap.cat=='low',],
                    aes(x=hap.n.diff , y=consistent, color=identities))+
        geom_smooth(data=toplot[toplot$topology=='PM'&toplot$min.n.hap.cat=='low',],
                    aes(x=hap.n.diff , y=consistent, color=identities))+
        geom_point(data=toplot[toplot$topology%in%c('MM','PM'),],
                   aes(x=hap.n.diff , y=consistent, color=identities, shape=topology, size=min.n.hap.cat))+
        theme_minimal() + theme(legend.position="bottom") + ylim(0,1.0) + 
        scale_color_manual(values=c('disagree'='#F24D29',both='#C4CFD0'))

ggarrange(annotate_figure(pci,  top = text_grob(model.name)), pD)


#######################################################
#Supplementary Figure 1. Distribution of the values of the covariates.
#Supplementary Figure 2.  Correlation matrix of the quantitative covariates.
#######################################################
data <- read.csv("direction_public/data/pairdata.csv")

#Set the covariates
data$response <- sapply(1:nrow(data), function(a) ifelse(data$p.ancestral.ER[a] <= 0.5, 'inconsistent', 'consistent'))
data$divergence.time <- as.integer((abs(data$sampling.calc.source.mid) + abs(data$sampling.calc.recipient.mid)))
data$min.n.hap <- sapply(1:nrow(data), function(a) min(c(data$hap.n.source[a], data$hap.n.rec[a])))
data$min.n.hap.cat <- sapply(1:nrow(data), function(a) ifelse(data$min.n.hap[a] < 10, 'low', 'high'))
data$pFr1.cat <- sapply(1:nrow(data), function(a) ifelse(data$pFr1[a] < 0.75, 'multiple', 'single'))
asinTransform <- function(p) asin(sqrt(p)) #arcsine square root transformation
data$pFr1 <- asinTransform(data$pFr1)
data$hap.n.diff <- (data$hap.n.source - data$hap.n.rec)
data$PD.UEH <- (data$PD.UEH.source - data$PD.UEH.recipient)
data$nuc.div <- (data$source.nuc.div - data$rec.nuc.div)
data$MinRootToTip <- (data$tree.MinRootToTip.source - data$tree.MinRootToTip.recipient)

#Supplementary Figure 1.
toplot <- data[,c('response', 'risk','status.at.transmission.source', 'pFr1.cat',
                  'min.n.hap.cat', 'topology', 'identities')]
toplot2 <- data[,c('response', 'hap.n.diff','divergence.time', 'alignment.length',
                   'nuc.div', 'MinRootToTip', 'PD.UEH', 'min.inter.TBL')]

p1 <-ggplot(toplot) + geom_bar(aes(x=status.at.transmission.source), fill="#999999") + theme_minimal() + ylab("Frequency") + xlab("Recency of the transmitter infection\n ")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))
p2 <-ggplot(toplot) + geom_bar(aes(x=risk), fill="#999999") + theme_minimal()+ ylab("Frequency") + xlab("Sexual risk exposure group\n ")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))
p3 <-ggplot(toplot2) + geom_histogram(aes(x=alignment.length), fill="#FC4E07") + theme_minimal()+ ylab("Frequency") + xlab("Sequence alignment length\n(number of base pairs)")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))
p4 <-ggplot(toplot2) + geom_histogram(aes(x=nuc.div), fill="#FC4E07") + theme_minimal()+ ylab("Frequency") + xlab(" Intra-host nucleotide diversity difference\n(substitutions per site)")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))
p5 <-ggplot(toplot) + geom_bar(aes(x=pFr1.cat), fill="#FC4E07") + theme_minimal()+ ylab("Frequency") + xlab("Multiplicity of infection\n ")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))
p6 <-ggplot(toplot) + geom_bar(aes(x=min.n.hap.cat), fill="#E7B800") + theme_minimal()+ ylab("Frequency") + xlab("Sample size\n(number of unique sampled sequences)")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))
p7 <-ggplot(toplot2) + geom_histogram(aes(x=hap.n.diff), fill="#E7B800") + theme_minimal()+ ylab("Frequency") + xlab("Sample size difference\n")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))
p8 <-ggplot(toplot2) + geom_histogram(aes(x=divergence.time), fill="#E7B800") + theme_minimal()+ ylab("Frequency") + xlab("Time from transmission\n(days)")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))
p9 <-ggplot(toplot) + geom_bar(aes(x=topology), fill="#31688EFF") + theme_minimal()+ ylab("Frequency") + xlab("Topology class\n ")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))
p10 <-ggplot(toplot) + geom_bar(aes(x=identities), fill="#31688EFF") + theme_minimal() + ylab("Frequency") + xlab("Most basal tip identity\n ")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))
p11 <-ggplot(toplot2) + geom_histogram(aes(x=MinRootToTip), fill="#31688EFF") + theme_minimal()+ ylab("Frequency") + xlab("Root-to-Tip difference\n(substitutions per site)")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))
p12 <-ggplot(toplot2) + geom_histogram(aes(x=PD.UEH), fill="#31688EFF") + theme_minimal()+ ylab("Frequency") + xlab("Phylogenetic diversity difference\n(substitutions per site)")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))
p13 <-ggplot(toplot2) + geom_histogram(aes(x=min.inter.TBL), fill="#31688EFF") + theme_minimal()+ ylab("Frequency") + xlab("Inter-host patristic distance\n(substitutions per site)")+theme(axis.text=element_text(size=7),axis.title=element_text(size=8))

ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,
          ncol = 3, nrow = 5, align = "v")

#setting up data for Supplementary Figure 2
toPlot <- data[, c("nuc.div",
                   "hap.n.diff",
                   "MinRootToTip",
                   "PD.UEH",	
                   "min.inter.TBL")]
names(toPlot) <- c('Intra-host nucleotide\ndiversity difference',
                   'Sample size difference',
                   'Root-to-tip difference',
                   'Phylogenetic diversity\ndifference',
                   'Inter-host patristic\ndistance')

#Supplementary Figure 2. Panel A.
testRes = cor.mtest(toPlot, conf.level = 0.95)
corrplot(cor(toPlot), type="upper", 
         col=brewer.pal(n=8, name="RdYlBu"), 
         diag=FALSE,
         p.mat = testRes$p, addCoef.col ='black',
         tl.col = 'black')

#Supplementary Figure 2. Panel B.
toPlotabs <- as.data.frame(sapply(toPlot, abs)) 
testResabs = cor.mtest(toPlotabs, conf.level = 0.95)
corrplot(cor(toPlotabs), type="upper", 
         col=brewer.pal(n=8, name="RdYlBu"), 
         diag=FALSE,
         p.mat = testResabs$p, addCoef.col ='black',
         tl.col = 'black')


#######################################################
#Supplementary Figure 4. Ordinal models outcomes when using routinely-available data. 
#######################################################
o0.6 <- read.csv("direction_public/models/ordinal_l_0.4-0.6_probTRUE_identityFALSE.csv")
o0.95 <- read.csv("direction_public/models/ordinal_l_0.05-0.95_probTRUE_identityFALSE.csv")
oPar <- read.csv("direction_public/models/ordinal_l_0.5-0.5_probFALSE_identityFALSE.csv")

o0.6$auc.equivalent <- (max(o0.6$auc)-o0.6$auc)<0.05
o0.95$auc.equivalent <- (max(o0.95$auc)-o0.95$auc)<0.05
oPar$auc.equivalent <- (max(oPar$auc)-oPar$auc)<0.05

o0.6$threshold <- '0.4_0.6'
o0.95$threshold <- '0.05_0.95'
oPar$threshold <- 'MP'

ord <- rbind(o0.6,o0.95)
ord <- rbind(ord,oPar)

ord$type <- factor(ord$type, levels = c('E','S','G','P','ES','EG','EP','SG','SP','GP','ESG','ESP','EGP','SGP','ESGP'))

ggplot(ord,aes(x=auc, y=type, xmin=auc.min, xmax=auc.max, color=auc.equivalent))+
        geom_point(aes(size=n.covariates),position = position_dodge(0.5)) + 
        geom_linerange(position = position_dodge(0.5)) +
        theme_minimal() + xlim(c(0.4,1)) + 
        scale_color_manual(values=c('TRUE'="#41ab5d",
                                    'FALSE'="#000000")) + coord_flip() +
        facet_grid(~threshold,scales = "free_x") + theme(panel.spacing = unit(2, "lines"),
                                                         legend.position="bottom") +
        theme(axis.text.x = element_text(angle = 45)) + xlab('AUC') +
        guides(colour = "none", size=guide_legend(title="Number of covariates")) +
        scale_size_continuous(range = c(1, 5), breaks=c(2,4,6,8))

