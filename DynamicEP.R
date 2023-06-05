rm(list=ls())
# Set the directory to the path of the current script
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Functions.R")
seed=123
set.seed(seed)

library("TruncatedNormal"); library("mvtnorm"); library("Matrix")
library("ggplot2"); library("latex2exp")


load("Financial.RData")


p = ncol(matFF)
n = length(y)


mu0 = matrix(0, nrow = p, ncol = 1)
Sigma0 = diag(rep(3,p))

# Smoothing IID
nSim = 1e4

SigmaEps = diag(c(0.01,0.01))
GG = diag(p)

startTime = Sys.time()
smoot = smoothing(y=y,
                  mu0 = mu0,
                  Sigma0 = Sigma0,
                  SigmaEps = SigmaEps,
                  matFF = matFF,
                  GG = GG, seed = seed, nSim = nSim)

meanTheta1     = apply(X=smoot$smoothState$values[,,1], MARGIN=2, FUN = mean)
meanTheta2     = apply(X=smoot$smoothState$values[,,2], MARGIN=2, FUN = mean)

sdTheta1     = apply(X=smoot$smoothState$values[,,1], MARGIN=2, FUN = sd)
sdTheta2     = apply(X=smoot$smoothState$values[,,2], MARGIN=2, FUN = sd)

timeMC = difftime(Sys.time(), startTime, units=("secs"))[[1]]


smoothParams = getSmooothParams(y,
                                mu0 = mu0,
                                Sigma0 = Sigma0,
                                SigmaEps = SigmaEps,
                                matFF = matFF,
                                GG =GG)
Omega = smoothParams$Omega

### Smoothing PMF-VB
startTime = Sys.time()
paramsPFM = getParamsPFM(X = X,
                         y = y,
                         Omega = Omega,
                         moments = T,
                         tolerance = 1e-3)
timePFM = difftime(Sys.time(), startTime, units=("secs"))[[1]]


indSeq1 = seq(from = 1, by = 2, length.out = n)
indSeq2 = indSeq1 + 1

meanTheta1_PFM = paramsPFM$postMoments.meanBeta[indSeq1]
meanTheta2_PFM = paramsPFM$postMoments.meanBeta[indSeq2]

sdTheta1_PFM = sqrt(paramsPFM$postMoments.varBeta[indSeq1])
sdTheta2_PFM = sqrt(paramsPFM$postMoments.varBeta[indSeq2])



### Smoothing MF-VB
startTime = Sys.time()
paramsMF = getParamsMF(X = X,
                        y = y,
                        Omega = Omega,
                        tolerance = 1e-3)
timeMF = difftime(Sys.time(), startTime, units=("secs"))[[1]]

meanTheta1_MF = paramsMF$meanBeta[indSeq1]
meanTheta2_MF = paramsMF$meanBeta[indSeq2]

sdTheta1_MF = sqrt(paramsMF$diagV[indSeq1])
sdTheta2_MF = sqrt(paramsMF$diagV[indSeq2])



### Smoothing EP
startTime = Sys.time()
paramsEP = getParamsEP(X = X,
                         y = y,
                         Omega = Omega,
                         tolerance = 1e-3)
timeEP = difftime(Sys.time(), startTime, units=("secs"))[[1]]

meanTheta1_EP = paramsEP$meanBeta[indSeq1]
meanTheta2_EP = paramsEP$meanBeta[indSeq2]

sdTheta1_EP = sqrt(paramsEP$diagOmegaEP[indSeq1])
sdTheta2_EP = sqrt(paramsEP$diagOmegaEP[indSeq2])



# Comparison different methods smoothing trajectories
Mean = c(meanTheta1,meanTheta2,
         meanTheta1_PFM,meanTheta2_PFM,
         meanTheta1_MF,meanTheta2_MF,
         meanTheta1_EP,meanTheta2_EP) 
Low  = c(meanTheta1 - sdTheta1,meanTheta2 - sdTheta2,
         meanTheta1_PFM - sdTheta1_PFM, meanTheta2_PFM - sdTheta2_PFM,
         meanTheta1_MF - sdTheta1_MF, meanTheta2_MF - sdTheta2_MF,
         meanTheta1_EP - sdTheta1_EP, meanTheta2_EP - sdTheta2_EP)
Upp  = c(meanTheta1 + sdTheta1,meanTheta2 + sdTheta2,
         meanTheta1_PFM + sdTheta1_PFM, meanTheta2_PFM + sdTheta2_PFM,
         meanTheta1_MF + sdTheta1_MF, meanTheta2_MF + sdTheta2_MF,
         meanTheta1_EP + sdTheta1_EP, meanTheta2_EP + sdTheta2_EP)
Par        = as.factor(rep(c(rep(1,n),rep(2,n)),4))
Method  = as.factor(rep(c("IID","PMF-VB","MF-VB","EP"),each=2*n))
Time = rep(1:n,8)

Data_plot = data.frame(Mean,Low,Upp,Par,Method,Time)

# New label names for plot
levels(Data_plot$Par) = c("1"=TeX("$\\theta_{1}$"), "2"=TeX("$\\theta_{2}$"))
Data_plot$Method = factor(Data_plot$Method,
                          levels = c("IID", "EP", "PMF-VB", "MF-VB"))

Plot_smooth = ggplot(Data_plot,aes(x=Time,col=Method))+
  geom_line(aes(y=Mean),size=1.2)+
  geom_line(aes(y=Low),linetype = "dashed",size=1.2)+
  geom_line(aes(y=Upp),linetype = "dashed",size=1.2)+
  scale_colour_manual(values = alpha(c("black", "red","blue","green"),
                                     c(1,1,1,1)))+theme_bw()+
  theme(axis.title=element_blank(),legend.title = element_blank(),
        axis.text=element_text(size=20),strip.text = element_text(size=30),
        legend.text = element_text(size=30))+
  facet_wrap(~Par,nrow=2,scales = "free_y",labeller=label_parsed,
             strip.position = "right")

Plot_smooth
ggsave(filename='Plot_smooth.png', plot=Plot_smooth, device="png",
       width = 30, height = 12, units="cm")



### Comparison different methods via box-plot
PFM_Diff_mean_Theta1 = meanTheta1 - meanTheta1_PFM
PFM_Diff_mean_Theta2 = meanTheta2 - meanTheta2_PFM

MF_Diff_mean_Theta1 = meanTheta1 - meanTheta1_MF
MF_Diff_mean_Theta2 = meanTheta2 - meanTheta2_MF

EP_Diff_mean_Theta1 = meanTheta1 - meanTheta1_EP
EP_Diff_mean_Theta2 = meanTheta2 - meanTheta2_EP

PFM_Diff_logsd_Theta1 = log(sdTheta1) - log(sdTheta1_PFM)
PFM_Diff_logsd_Theta2 = log(sdTheta2) - log(sdTheta2_PFM)

MF_Diff_logsd_Theta1 = log(sdTheta1) - log(sdTheta1_MF)
MF_Diff_logsd_Theta2 = log(sdTheta2) - log(sdTheta2_MF)

EP_Diff_logsd_Theta1 = log(sdTheta1) - log(sdTheta1_EP)
EP_Diff_logsd_Theta2 = log(sdTheta2) - log(sdTheta2_EP)

Diff_Mean  = c(PFM_Diff_mean_Theta1, PFM_Diff_mean_Theta2, MF_Diff_mean_Theta1, 
               MF_Diff_mean_Theta2,EP_Diff_mean_Theta1, EP_Diff_mean_Theta2)
Diff_logsd = c(PFM_Diff_logsd_Theta1,PFM_Diff_logsd_Theta2,MF_Diff_logsd_Theta1, 
               MF_Diff_logsd_Theta2,EP_Diff_logsd_Theta1, EP_Diff_logsd_Theta2)


Par        = as.factor(rep(rep(c(rep(1,n),rep(2,n)),3),2))
Method     = as.factor(rep(rep(c("PFM-VB","MF-VB","EP"),each=2*n),2))
Functional = as.factor(rep(1:2,each=6*n))

Diff       = c(Diff_Mean,Diff_logsd)

Data_boxplot = data.frame(Diff,Method,Par,Functional)
col_meth   = c("red","blue","green")

# New label names for Plot
levels(Data_boxplot$Par)= c("1"=TeX("$\\theta_{1}$"), "2"=TeX("$\\theta_{2}$"))
levels(Data_boxplot$Functional) = 
  c("1"=TeX("$\\Delta E(\\theta_{jt} | y_{1 : n})$"), 
    "2"=TeX("$\\Delta \\log(\\sqrt{var(\\theta_{jt} | y_{1 : n})})$"))
Data_boxplot$Method = factor(Data_boxplot$Method, c("EP","PFM-VB","MF-VB"))

Plot_Diff = ggplot(Data_boxplot, aes(y=Diff, x=Method,col=Method)) + 
  geom_boxplot()+
  scale_colour_manual(values = alpha(col_meth,c(1,1,1)))+
  theme_bw()+ geom_hline(yintercept=0, linetype="dashed", size=2)+
  theme(axis.title=element_blank(),axis.text=element_text(size=20) , 
        axis.text.x = element_text(colour = col_meth), 
        strip.text = element_text(size=25),legend.position = "none")+
  facet_grid(Functional~Par,scales = "free", labeller=label_parsed)
Plot_Diff
ggsave(filename='Plot_Diff.png', plot=Plot_Diff, device="png",
       width = 30, height = 20, units="cm")

### Comparison different methods: time
timeMC
timeEP
timeMF
timePFM

