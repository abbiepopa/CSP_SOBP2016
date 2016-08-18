#Set WD#
setwd("~/Documents/Lab/SOBP2016/Data")

#Import Freesurfer Volume#
vol <- read.csv("freesurfer_volume_20150505-NJH.csv")

#Separate into two files: left and right
###note, for some reason 2RO5NormBR only has left
vol_r<-vol[which(vol$hemisphere == "right"),]
vol_l<-vol[which(vol$hemisphere == "left"),]

#Import CSP
CSP <- read.csv("CSP_Cycle2_MASTER_2014_AMP.csv", na.strings="!")
CSP_s <- CSP[,c("CABILID","STUDYID","VOLUME","LENGTH")]

#examine data to see if we need to log transform
library(psych)
describe(vol_r)[,c("skew","kurtosis")]
hist(vol_r$insula)
hist(vol_r$entorhinal)
#looks pretty okay, insula kurtosis = 1.59, entorhinal = 1.40

vol_l$parahippocampal[which(vol_l$parahippocampal == 0)]<-NA
describe(vol_l)[,c("skew","kurtosis")]
#caudateanteriorcingulate skew = 1.06, kurtosis = 2.33
#fusiform skew = 2.18, kurtosis = 15.72
#isthmuscingulate skew = 0.45, kurtosis = 2.95
#lateralorbitofrontal skew = -0.83, kurtosis = 1.28
#superiortemporal skew = -0.56, kurtosis = 1.61
#insula skew = 0.12, kurtosis =1.31

hist(vol_l$caudalanteriorcingulate)
hist(vol_l$fusiform) #def pointy
hist(vol_l$isthmuscingulate)
hist(vol_l$lateralorbitofrontal)
hist(vol_l$superiortemporal)
hist(vol_l$insula)

#based on examination of histograms a cut-off of 2 seems reasonable, though really only fusiform looks terrible
vol_l$log_caudalanteriorcingulate<-log(vol_l$caudalanteriorcingulate)
vol_l$log_fusiform<-log(vol_l$fusiform)
vol_l$log_isthmuscingulate<-log(vol_l$isthmuscingulate)

describe(CSP_s)[,c("skew","kurtosis")]
#Volume skew = 3.77, kurtosis = 15.43
#Length skew = 5.96, kurtosis = 50.26

hist(CSP_s$VOLUME)
hist(CSP_s$LENGTH)

CSP_s$log_volume<-log(CSP_s$VOLUME)
CSP_s$log_length<-log(CSP_s$LENGTH)

#Merge Files
colnames(CSP_s)[2]<-"studyid"
r<-merge(vol_r, CSP_s)
l<-merge(vol_l, CSP_s)
l_c<-na.omit(l)
r_c<-na.omit(r)

#for cor.test, statistic is t, p.value is p, estimate is r
#bankssts<-c(cor.test(l_c$log_volume, l_c$bankssts)$p.value, 
#	cor.test(l_c$log_volume, l_c$bankssts)$estimate)

l_cor <- data.frame()
for (i in 4:40)
{
	n <- cor.test(l_c$log_volume, l_c[,i])
	l_cor<-rbind(l_cor,c(n$p.value, n$estimate))
	 }
	 
row.names(l_cor)<-colnames(l_c)[4:40]
colnames(l_cor)<-c("p","r")

sig_l<-l_cor[which(l_cor$p < 0.05),]

r_cor<-data.frame()
for (i in 4:37)
{
	n <- cor.test(r_c$log_volume, r_c[,i])
	r_cor<-rbind(r_cor,c(n$p.value, n$estimate))
}

row.names(r_cor)<-colnames(r_c)[4:37]
colnames(r_cor)<-c("p","r")

sig_r<-r_cor[which(r_cor$p < 0.05),]

sig_l$hemi <- "left"
sig_r$hemi <- "right"

sig <- rbind(sig_l, sig_r)

sig <- sig[order(sig$p),]

sig$p_fdr<-p.adjust(sig$p, method = "fdr", n = (length(r_cor$p) + length(l_cor$p)))

#Dx#
l_c$Dx <- "22q"
l_c$Dx[grep("Norm",l_c$studyid, fixed = T)]<-"TD"
l_c$Dx[grep("SCA", l_c$studyid, fixed = T)]<-"SCA"


r_c$Dx <- "22q"
r_c$Dx[grep("Norm",r_c$studyid, fixed = T)]<-"TD"
r_c$Dx[grep("SCA", r_c$studyid, fixed = T)]<-"SCA"

####Do correlations for brain volume normalized versions of l_c and r_c####
ICV <- read.csv("Brain_ICV.csv")
r_c_norm<-r_c
ICVs<-ICV[,c("cabilid","subjectid","ICV_count","ICV_volume..mm.3.")]

#####################################################################
###Naomi's sheets are volume, so volume is the better normalizer!!###
#####################################################################

r_c_norm<-merge(r_c_norm, ICVs, by.x="studyid", by.y="subjectid", all.x=T)

for( i in 4:37){
	r_c_norm[,i]<-r_c_norm[,i]/r_c_norm[,"ICV_volume..mm.3."]
}

###################################
####Added to Norm CSP on 160501####
###################################

r_c_norm$VOLUME <- r_c_norm$VOLUME/r_c_norm$ICV_volume..mm.3.
r_c_norm$log_volume<-log(r_c_norm$VOLUME)

###################################
####Added to Norm CSP on 160501####
###################################

l_c_norm<-l_c

l_c_norm<-merge(l_c_norm, ICVs, by.x="studyid",by.y="subjectid")

for(i in 4:40){
	l_c_norm[,i]<-l_c_norm[,i]/l_c_norm[,"ICV_volume..mm.3."]
}

###################################
####Added to Norm CSP on 160501####
###################################

l_c_norm$VOLUME <- l_c_norm$VOLUME/l_c_norm$ICV_volume..mm.3.
l_c_norm$log_volume<-log(l_c_norm$VOLUME)

###################################
####Added to Norm CSP on 160501####
###################################

###################
####Normed Corrs###
###################

l_cor_norm <- data.frame()
for (i in 4:40)
{
	n <- cor.test(l_c_norm$log_volume, l_c_norm[,i])
	l_cor_norm<-rbind(l_cor_norm,c(n$p.value, n$estimate))
	 }
	 
row.names(l_cor_norm)<-colnames(l_c_norm)[4:40]
colnames(l_cor_norm)<-c("p","r")

sig_l_norm<-l_cor_norm[which(l_cor_norm$p < 0.05),]

r_cor_norm<-data.frame()
for (i in 4:37)
{
	n <- cor.test(r_c_norm$log_volume, r_c_norm[,i])
	r_cor_norm<-rbind(r_cor_norm,c(n$p.value, n$estimate))
}

row.names(r_cor_norm)<-colnames(r_c_norm)[4:37]
colnames(r_cor_norm)<-c("p","r")

sig_r_norm<-r_cor_norm[which(r_cor_norm$p < 0.05),]

sig_l_norm$hemi <- "left"
sig_r_norm$hemi <- "right"

sig_norm <- rbind(sig_l_norm, sig_r_norm)

sig_norm <- sig_norm[order(sig_norm$p),]

sig_norm$p_fdr<-p.adjust(sig_norm$p, method = "fdr", n = (length(r_cor_norm$p) + length(l_cor_norm$p)))

#######################################
###Violin Plots of CSP Size by Group###
#######################################
library (car)
Anova(lm(log_volume~Dx, data=l_c))
#Response: log_volume
#          Sum Sq  Df F value    Pr(>F)    
#Dx        108.75   2  15.139 1.078e-06 ***
#Residuals 517.20 144  

Anova(lm(log_volume~Dx, data=l_c_norm))
#Response: log_volume
#          Sum Sq  Df F value    Pr(>F)    
#Dx        115.57   2  16.538 3.423e-07 ***
#Residuals 503.12 144  

Anova(lm(log_volume~Dx, data=r_c))
#Response: log_volume
#          Sum Sq  Df F value    Pr(>F)    
#Dx        109.17   2  15.104 1.119e-06 ***
#Residuals 516.78 143  

Anova(lm(log_volume~Dx, data=r_c_norm))
#Response: log_volume
#          Sum Sq  Df F value    Pr(>F)    
#Dx        115.98   2  16.496 3.578e-07 ***
#Residuals 502.69 143    

pairwise.t.test(l_c$log_volume, l_c$Dx, p.adjust.method="fdr")
#    22q     SCA 
#SCA 2.2e-05 -   
#TD  2.5e-05 0.62

pairwise.t.test(l_c_norm$log_volume, l_c_norm$Dx, p.adjust.method="fdr")
#    22q     SCA 
#SCA 8.5e-06 -   
#TD  8.5e-06 0.68

pairwise.t.test(r_c$log_volume, r_c$Dx, p.adjust.method="fdr")
#    22q     SCA 
#SCA 2.3e-05 -   
#TD  2.6e-05 0.65

pairwise.t.test(r_c_norm$log_volume, r_c_norm$Dx, p.adjust.method="fdr")
#    22q     SCA
#SCA 8.7e-06 -  
#TD  8.7e-06 0.7

library(ggplot2)
ggplot(l_c, aes(Dx, log_volume, fill=Dx)) + geom_violin() + scale_fill_manual(values=c("hotpink","forestgreen","blue"))+geom_boxplot(width=0.05, col="gray", fill="black")+labs(title="Not Normalized for ICV", x = "Diagnostic Group", y="CSP Volume (log)")

quartz()
ggplot(l_c_norm, aes(Dx, log_volume, fill=Dx)) + geom_violin() + scale_fill_manual(values=c("hotpink","forestgreen","blue"))+geom_boxplot(width=0.05, col="gray", fill="black")+labs(title="Normalized for ICV", x = "Diagnostic Group", y="CSP Volume (log)")

############################################################
###identify who has extreme CSPs using different criteria###
############################################################
l_c$CSP5<-"Less Than 5mm"
l_c[which(l_c$LENGTH>5),"CSP5"]<-"Greater Than 5mm"

l_c$CSP6<-"Less Than 6mm"
l_c[which(l_c$LENGTH>6),"CSP6"]<-"Greater Than 6mm"

l_c$CSPml<-"Less Than Median Length"
l_c[which(l_c$LENGTH>4), "CSPml"]<-"Greater Than Median Length"

##################################################
###Correlation in The >6 group and the <6 group###
##################################################

###Example Plot###
pdf("rACC.pdf")
plot(r_c[which(r_c$Dx == "22q"),"log_volume"], r_c[which(r_c$Dx=="22q"),"rostralanteriorcingulate"], pch=16, col="hotpink", main = "CSP Volume Negatively Correlated with Right rACC Volume", xlab = "Volume (log) CSP", ylab = "Volume Right rACC", pin=c(5,5))
points(r_c[which(r_c$Dx == "TD"),"log_volume"], r_c[which(r_c$Dx=="TD"),"rostralanteriorcingulate"], pch=16, col="blue")
points(r_c[which(r_c$Dx == "SCA"),"log_volume"], r_c[which(r_c$Dx=="SCA"),"rostralanteriorcingulate"], pch=16, col="green")
fit <- lm(rostralanteriorcingulate ~ log_volume, data = r_c)
abline(fit, col="darkgrey", lwd = 4)
legend(6, 3500, legend=c("22q n=67","TD n=44","SCA n=35", "r=-0.22,", "p_fdr=0.07,","p_uncor=0.006"), pch = c(16,16,16,NA, NA, NA), col = c("hotpink","blue","green","darkgrey", NA, NA), lwd=c(NA,NA,NA,4, NA, NA))
dev.off()

pdf("cuneus.pdf")
#quartz()
plot(l_c[which(l_c$Dx == "22q"),"log_volume"], l_c[which(l_c$Dx=="22q"),"cuneus"], pch=16, col="hotpink", main = "CSP Volume Negatively Correlated with Left Cuneus Volume", xlab = "Volume (log) CSP", ylab = "Volume Left Cuneus", pin=c(5,5))
points(l_c[which(l_c$Dx == "TD"),"log_volume"], l_c[which(l_c$Dx=="TD"),"cuneus"], pch=16, col="blue")
points(l_c[which(l_c$Dx == "SCA"),"log_volume"], l_c[which(l_c$Dx=="SCA"),"cuneus"], pch=16, col="green")
fit <- lm(cuneus ~ log_volume, data = l_c)
abline(fit, col="darkgrey", lwd = 4)
legend(4.5, 1900, legend=c("22q n=67","TD n=45","SCA n=35", "r=-0.32,", "p_fdr=0.007,","p_uncor=0.00009"), pch = c(16,16,16,NA, NA, NA), col = c("hotpink","blue","green","darkgrey", NA, NA), lwd=c(NA,NA,NA,4, NA, NA))
dev.off()

write.csv(sig, "sig.csv")
write.csv(sig_norm, "sig_norm.csv")