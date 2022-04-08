#######################################
### Author: Dr. Nikolaos Trasanidis ### 
#######################################

#SURVIVAL ANALYSES FOR FIGURE1

library("survival")
library("survminer")
library(boot)

#############################################
### STEP1. chr1q-WIDE UNIVARIATE ANALYSES ###
#############################################

#1A. CNV analysis WES data #
############################################################

CNV_1q<-read.table("1q_CNV_WES_amplif_SURVIVAL.zip", header=TRUE)

Acovariates<-colnames(CNV_1q[,-(1:3), drop=FALSE])
univ_formulas_C <- sapply(Acovariates,
                          function(x) as.formula(paste('Surv(Time, Flag)~', x)))

univ_models_C <- lapply( univ_formulas_C, function(x){coxph(x, data = CNV_1q)})
# Extract data 
univ_results_C <- lapply(univ_models_C,
                         function(x){ 
                           x <- summary(x)
                           wald.p<-signif(x$wald["pvalue"], digits=3)
                           wald.test<-signif(x$wald["test"], digits=3)
                           logrank.p<-signif(x$sctest["pvalue"], digits=3)
                           logrank.test<-signif(x$sctest["test"], digits=3)
                           likel.test<-signif(x$logtest["test"], digits=3)
                           likel.p<-signif(x$logtest["pvalue"], digits=3)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR <- paste0(HR, " (", 
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res_CNV<-c(beta, HR, wald.test, wald.p, logrank.test, logrank.p, likel.test, likel.p)
                           names(res_CNV)<-c("beta", "HR (95% CI for HR)", "wald.test", "wald.p", "logrank.test", "logrank.p", "likel.test", "likel.p")
                           return(res_CNV)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
res_CNV <- t(as.data.frame(univ_results_C, check.names = FALSE))
res_CNV2<-as.data.frame(res_CNV)
write.table(res_CNV2,file="Survival_results_CNV.txt", sep="\t")

###### Select only genes displaying p<0.05 in logrank test
res_CNV2$names<-rownames(res_CNV2)
res_CNV2$logrank.p<-as.numeric(as.character( res_CNV2$logrank.p ))
res_CNV3<-subset(res_CNV2,logrank.p < 0.05)

Amcovariates <- res_CNV3$names
Amcovariates2<-do.call(paste,c(as.list(Amcovariates), sep="+"))

multiv_results_CNV<-as.formula(paste('coxph(Surv(Time, Flag) ~', Amcovariates2,', data=CNV_1q)'))




#1B. CNV analysis WGS data #
#####################################################################

CNVG_1q<-read.table("1q_CNV_WGS_amplif_SURVIVAL.zip", header=TRUE)

AGcovariates<-colnames(CNVG_1q[,-(1:3), drop=FALSE])
univ_formulas_CG <- sapply(AGcovariates,
                           function(x) as.formula(paste('Surv(Time, Flag)~', x)))

univ_models_CG <- lapply( univ_formulas_CG, function(x){coxph(x, data = CNVG_1q)})
# Extract data 
univ_results_CG <- lapply(univ_models_CG,
                          function(x){ 
                            x <- summary(x)
                            wald.p<-signif(x$wald["pvalue"], digits=3)
                            wald.test<-signif(x$wald["test"], digits=3)
                            logrank.p<-signif(x$sctest["pvalue"], digits=3)
                            logrank.test<-signif(x$sctest["test"], digits=3)
                            likel.test<-signif(x$logtest["test"], digits=3)
                            likel.p<-signif(x$logtest["pvalue"], digits=3)
                            beta<-signif(x$coef[1], digits=2);#coeficient beta
                            HR <-signif(x$coef[2], digits=2);#exp(beta)
                            HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                            HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                            HR <- paste0(HR, " (", 
                                         HR.confint.lower, "-", HR.confint.upper, ")")
                            res_CNVG<-c(beta, HR, wald.test, wald.p, logrank.test, logrank.p, likel.test, likel.p)
                            names(res_CNVG)<-c("beta", "HR (95% CI for HR)", "wald.test", "wald.p", "logrank.test", "logrank.p", "likel.test", "likel.p")
                            return(res_CNVG)
                            #return(exp(cbind(coef(x),confint(x))))
                          })
res_CNVG <- t(as.data.frame(univ_results_CG, check.names = FALSE))
res_CNVG2<-as.data.frame(res_CNVG)
write.table(res_CNVG2,file="Survival_results_CNV_WGS.txt", sep="\t")

###### Select only genes displaying p<0.05 in logrank test
res_CNVG2$names<-rownames(res_CNVG2)
res_CNVG2$logrank.p<-as.numeric(as.character( res_CNVG2$logrank.p ))
res_CNVG3<-subset(res_CNVG2,logrank.p < 0.05)

AGmcovariates <- res_CNVG3$names
AGmcovariates2<-do.call(paste,c(as.list(AGmcovariates), sep="+"))

multiv_results_CNVG<-as.formula(paste('coxph(Surv(Time, Flag) ~', AGmcovariates2,', data=CNVG_1q)'))



#1C. RNA analysis using TPM scores #
######################################################################

TPM_1q<-read.table("1q_TPM_expression_SURVIVAL.zip", header=TRUE)

Ecovariates<-colnames(TPM_1q[,-(1:3), drop=FALSE])
univ_formulas <- sapply(Ecovariates,
                        function(x) as.formula(paste('Surv(Time, Flag)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = TPM_1q)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         wald.p<-signif(x$wald["pvalue"], digits=3)
                         wald.test<-signif(x$wald["test"], digits=3)
                         logrank.p<-signif(x$sctest["pvalue"], digits=3)
                         logrank.test<-signif(x$sctest["test"], digits=3)
                         likel.test<-signif(x$logtest["test"], digits=3)
                         likel.p<-signif(x$logtest["pvalue"], digits=3)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res_TPM<-c(beta, HR, wald.test, wald.p, logrank.test, logrank.p, likel.test, likel.p)
                         names(res_TPM)<-c("beta", "HR (95% CI for HR)", "wald.test", "wald.p", "logrank.test", "logrank.p", "likel.test", "likel.p")
                         return(res_TPM)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res_TPM <- t(as.data.frame(univ_results, check.names = FALSE))
res_TPM2<-as.data.frame(res_TPM)
write.table(res_TPM2,file="Survival_results_TPM.txt", sep="\t")

###### Select only genes displaying p<0.05 in logrank test
res_TPM2$names<-rownames(res_TPM2)
res_TPM2$logrank.p<-as.numeric(as.character( res_TPM2$logrank.p ))
res_TPM3<-subset(res_TPM2,logrank.p < 0.05)

##### Multivariate analysis###
Emcovariates <- res_TPM3$names
Emcovariates2<-do.call(paste,c(as.list(Emcovariates), sep="+"))

multiv_results<-as.formula(paste('coxph(Surv(Time, Flag) ~', Emcovariates2,', data=TPM_1q)'))






###############################################
### STEP2. chr1q-WIDE MULTIVARIATE ANALYSES ###
###############################################

# First identify other genetic factors which are significantly associated with poor prognosis

#2A: Univariate analysis on WGS_SeqFISH parameters (IgG tranlocations)
WGS_SeqFISH<-read.table("MMRF_WGS_SeqFISH_CNV.txt", header=TRUE, fill=TRUE)
SeqFISH_cov<-colnames(WGS_SeqFISH[,-(1), drop=FALSE])
SeqFISH_cov1<-SeqFISH_cov[53:108]
CNVG_1q_Multiv<-merge.data.frame(CNVG_1q,WGS_SeqFISH,by="Patient_ID")


univ_SeqFISH_formulas_CG <- sapply(SeqFISH_cov,
                           function(x) as.formula(paste('Surv(Time, Flag)~', x)))

univ_SeqFISH_models_CG <- lapply( univ_SeqFISH_formulas_CG, function(x){coxph(x, data = CNVG_1q_Multiv)})
# Extract data 
univ_SeqFISH_results_CG <- lapply(univ_SeqFISH_models_CG,
                          function(x){ 
                            x <- summary(x)
                            wald.p<-signif(x$wald["pvalue"], digits=3)
                            wald.test<-signif(x$wald["test"], digits=3)
                            logrank.p<-signif(x$sctest["pvalue"], digits=3)
                            logrank.test<-signif(x$sctest["test"], digits=3)
                            likel.test<-signif(x$logtest["test"], digits=3)
                            likel.p<-signif(x$logtest["pvalue"], digits=3)
                            beta<-signif(x$coef[1], digits=2);#coeficient beta
                            HR <-signif(x$coef[2], digits=2);#exp(beta)
                            HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                            HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                            HR <- paste0(HR, " (", 
                                         HR.confint.lower, "-", HR.confint.upper, ")")
                            res_SeqFISH_CNVG<-c(beta, HR, wald.test, wald.p, logrank.test, logrank.p, likel.test, likel.p)
                            names(res_SeqFISH_CNVG)<-c("beta", "HR (95% CI for HR)", "wald.test", "wald.p", "logrank.test", "logrank.p", "likel.test", "likel.p")
                            return(res_SeqFISH_CNVG)
                            #return(exp(cbind(coef(x),confint(x))))
                          })
res_SeqFISH_CNVG <- t(as.data.frame(univ_SeqFISH_results_CG, check.names = FALSE))
res_SeqFISH_CNVG2<-as.data.frame(res_SeqFISH_CNVG)
write.table(res_SeqFISH_CNVG2,file="Survival_results_SeqFISH_WGS.txt", sep="\t")

res_SeqFISH_CNVG2$names<-rownames(res_SeqFISH_CNVG2)
res_SeqFISH_CNVG2$logrank.p<-as.numeric(as.character( res_SeqFISH_CNVG2$logrank.p ))
res_SeqFISH_CNVG3<-subset(res_SeqFISH_CNVG2,logrank.p < 0.05)

## Select sign. SeqFISH parameters 
SeqFISH_cov2<-rownames(res_SeqFISH_CNVG3)[1:15]


####2B: Univariate analysis on RNAseq_Ig_Translocations (IgG tranlocations called based on Oncogenes Overexpression)
WGS_IgTransl<-read.table("MMRF_RNAseq_Ig_Translocations.txt", header=TRUE)
WGS_IgTransl_cov<-colnames(WGS_IgTransl[,-(1), drop=FALSE])
CNVG_1q_Multiv2<-merge(CNVG_1q,WGS_IgTransl,by="Patient_ID")


univ_IgTransl_formulas_CG <- sapply(WGS_IgTransl_cov,
                                   function(x) as.formula(paste('Surv(Time, Flag)~', x)))

univ_IgTransl_models_CG <- lapply(univ_IgTransl_formulas_CG, function(x){coxph(x, data = CNVG_1q_Multiv2)})
# Extract data 
univ_IgTransl_results_CG <- lapply(univ_IgTransl_models_CG,
                                  function(x){ 
                                    x <- summary(x)
                                    wald.p<-signif(x$wald["pvalue"], digits=3)
                                    wald.test<-signif(x$wald["test"], digits=3)
                                    logrank.p<-signif(x$sctest["pvalue"], digits=3)
                                    logrank.test<-signif(x$sctest["test"], digits=3)
                                    likel.test<-signif(x$logtest["test"], digits=3)
                                    likel.p<-signif(x$logtest["pvalue"], digits=3)
                                    beta<-signif(x$coef[1], digits=2);#coeficient beta
                                    HR <-signif(x$coef[2], digits=2);#exp(beta)
                                    HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                                    HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                                    HR <- paste0(HR, " (", 
                                                 HR.confint.lower, "-", HR.confint.upper, ")")
                                    res_IgTransl_CNVG<-c(beta, HR, wald.test, wald.p, logrank.test, logrank.p, likel.test, likel.p)
                                    names(res_IgTransl_CNVG)<-c("beta", "HR (95% CI for HR)", "wald.test", "wald.p", "logrank.test", "logrank.p", "likel.test", "likel.p")
                                    return(res_IgTransl_CNVG)
                                    #return(exp(cbind(coef(x),confint(x))))
                                  })
res_IgTransl_CNVG <- t(as.data.frame(univ_IgTransl_results_CG, check.names = FALSE))
res_IgTransl_CNVG2<-as.data.frame(res_IgTransl_CNVG)
write.table(res_IgTransl_CNVG2,file="Survival_results_RNA_IgTransl_WGS.txt", sep="\t")

res_IgTransl_CNVG2$names<-rownames(res_IgTransl_CNVG2)
res_IgTransl_CNVG2$logrank.p<-as.numeric(as.character( res_IgTransl_CNVG2$logrank.p ))
res_IgTransl_CNVG3<-subset(res_IgTransl_CNVG2,logrank.p < 0.1)

## Select sign. SeqFISH parameters 
WGS_IgTransl_cov2<-rownames(res_IgTransl_CNVG3)



##2C: Multiple multivariate analyses
CNVG_1q_Multiv3<-merge(CNVG_1q_Multiv, WGS_IgTransl, by="Patient_ID")
SeqFISH_n_IgTransl_cov<-append(SeqFISH_cov2,WGS_IgTransl_cov2)

ALLcovariates<-do.call(paste,c(as.list(SeqFISH_n_IgTransl_cov), sep="+"))
ALLcovariates2<-paste("+",ALLcovariates)  


Multiv_formulas_IgFISH <- sapply(AGcovariates,
                             function(x) as.formula(paste('Surv(Time, Flag)~', x, ALLcovariates2)))
multiv_models_IgFISH <- lapply(Multiv_formulas_IgFISH, function(x){coxph(x, data = CNVG_1q_Multiv3)})
# Extract data 
multiv_results_IgFISH <- lapply(multiv_models_IgFISH,
                            function(x){ 
                              x <- summary(x) })

res_multiv_CNVG_IgFISH<-as.data.frame(multiv_results_IgFISH)
options(max.print=999999)
capture.output(multiv_results_IgFISH, file="Survival_results_CNV_WGS_Ig_n_FISH_multiv.txt")





##############################################################
### STEP3. SIMULATIONS ON chr1q-WIDE MULTIVARIATE ANALYSES ###
##############################################################

library("dplyr")
library(data.table)
library(boot)
library("survival")
library("survminer")
library(tidyr)

# Censboot integrated function
data=CNVG_1q_Multiv3
data[is.na(data)] <- 0

# Bootstrapping parameters settings
# R = No of iterations
# Sample = number of patients to subsample

R=10000 
Npatient=780

Pvalues<-data.frame(sample(0,2215,rep=TRUE))
HRs<-data.frame(sample(0,2215,rep=TRUE))
data0<- data.table(data)
#test<-data.frame()
for (i in 1:R){
  set.seed(i)
  data1<-data0[sample(.N, as.integer(Npatient+i/100000))]
  #data1<-sample_n(data1,860)
  Multiv_formulas_ALL_MC <- sapply(colnames(data1[,4:2218]),
                                   function(x) as.formula(paste('Surv(Time, Flag)~', x, ALLcovariates2)))
  Multiv_models_ALL_MC <- lapply(Multiv_formulas_ALL_MC, function(x){coxph(x, data1)})
  Multiv_results_ALL_MC <- lapply(Multiv_models_ALL_MC,
                                  function(x){ 
                                    x <- summary(x)})
  SummarizedTable<-data.frame()
  for (x in rownames(summary(Multiv_results_ALL_MC))){ 
    y <- as.data.frame(Multiv_results_ALL_MC[[x]]$coefficients)
    SummarizedTable <- rbind(SummarizedTable, y[1,])
  }
  #return(SummarizedTable) #no need
  Pvalues<-cbind(Pvalues, i=as.data.frame(SummarizedTable$`Pr(>|z|)`))
  rownames(Pvalues)<-rownames(SummarizedTable)
  HRs<-cbind(HRs,i=as.data.frame(SummarizedTable$`exp(coef)`))
  rownames(HRs)<-rownames(SummarizedTable)
  #SummarizedTable2<-cbind(SummarizedTable$`Pr(>|z|)`,SummarizedTable$`exp(coef)`)
  #test1<-data1[1,]
  #test<-rbind(test,test1)
}

Pvalues<-Pvalues[,-c(1)]
HRs<-HRs[,-c(1)]
colnames(Pvalues)<-c(1:R)
colnames(HRs)<-c(1:R)

Pvaluest<-as.data.frame(t(Pvalues))
Pvalues_overview <- gather(Pvaluest, factor_key=TRUE)
Poverview <- Pvalues_overview%>% group_by(key)%>%
  summarise(mean= mean(value), median= median(value), sd= sd(value), cv=sd(value)/mean(value), max = max(value),min = min(value))

HRst<-as.data.frame(t(HRs))
HRs_overview<- gather(HRst, factor_key=TRUE)
HRoverview<-HRs_overview%>% group_by(key)%>%
  summarise(mean= mean(value), median= median(value), sd= sd(value), cv=sd(value)/mean(value), max = max(value),min = min(value))

write.table(Pvalues,"WGS_simulations_Pvalues.txt", sep="\t")
write.table(Poverview,"WGS_simulations_Pvalues_overview.txt", sep="\t")
write.table(HRs,"WGS_simulations_HRs.txt", sep="\t")
write.table(HRoverview,"WGS_simulations_HRs_overview.txt", sep="\t")



