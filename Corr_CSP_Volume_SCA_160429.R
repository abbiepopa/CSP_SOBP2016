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

####for TD grep("Norm",r_c$studyid)####
####for SCA grep("SCA", r_c$studyid)####
###for 22q... probably have to remove the others###

l_c$Dx<-"22q"
l_c[grep("Norm",l_c$studyid),"Dx"]<-"TD"
l_c[grep("SCA",l_c$studyid),"Dx"]<-"SCA"

r_c$Dx<-"22q"
r_c[grep("Norm",r_c$studyid),"Dx"]<-"TD"
r_c[grep("SCA",r_c$studyid),"Dx"]<-"SCA"

###################################################
####Change this part to change diagnostic group####
###################################################

l_c<-l_c[which(l_c$Dx == "SCA"),]
r_c<-r_c[which(r_c$Dx == "SCA"),]

#####################
####Then Run This####
#####################


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

########################
###Add Dx to Filename###
########################

write.csv(sig, "sig_SCA.csv")
write.csv(sig_norm, "sig_norm_SCA.csv")