#***************INPUT file name
## provide input file for the calculation of mean Divergence, X/A ratio and Alpha values
input <- "hg18.farFromGenes.20kb.nosep.schweinfurthii.bothsexes.reheader"   # input file name

## provide input directory where the bootstrap files are located (without the 'Afiles' and 'Xfiles' extension)
## Following output files are written: $input_Divergence (with mean and CI of Divergence in each bin)
##                                    $input_XA (with mean and CI of X/A Divergence ratio in each bin)
##                                    $input_Alpha (with mean and CI of Alpha (MMB) in each bin)

writefile <- "Schweinfurthii"

#*****************************MEAN CALCULATION*************************
#----Read input file---------------------
data <- read.table(input,header=TRUE,sep="\t")

#----Filter data based on Numbases and Dist2Root(divergence) values----------------
selectedX <- subset(subset(subset(data,Dist2Root > 0), NumBases > 5000),chrom == "chrX")
selectedA <- subset(subset(subset(data,Dist2Root > 0), NumBases > 5000),chrom != "chrX")

#----Divide filtered data into bins(cm2genes)---------------------
Xdatabin1 <- subset(selectedX,cm2genes >= 0 & cm2genes < 0.05)
Xdatabin2 <- subset(selectedX,cm2genes >= 0.05 & cm2genes < 0.1)
Xdatabin3 <- subset(selectedX,cm2genes >= 0.1 & cm2genes < 0.2)
Xdatabin4 <- subset(selectedX,cm2genes >= 0.2 & cm2genes < 0.4)
Xdatabin5 <- subset(selectedX,cm2genes >= 0.4 & cm2genes < 0.8)
Xdatabin6 <- subset(selectedX,cm2genes >= 0.8 & cm2genes < 2.0)
Adatabin1 <- subset(selectedA,cm2genes >= 0 & cm2genes < 0.05)
Adatabin2 <- subset(selectedA,cm2genes >= 0.05 & cm2genes < 0.1)
Adatabin3 <- subset(selectedA,cm2genes >= 0.1 & cm2genes < 0.2)
Adatabin4 <- subset(selectedA,cm2genes >= 0.2 & cm2genes < 0.4)
Adatabin5 <- subset(selectedA,cm2genes >= 0.4 & cm2genes < 0.8)
Adatabin6 <- subset(selectedA,cm2genes >= 0.8 & cm2genes < 2.0)


meanXAbin1 <- mean(Xdatabin1$Dist2Root)/mean(Adatabin1$Dist2Root)
meanXAbin2 <- mean(Xdatabin2$Dist2Root)/mean(Adatabin2$Dist2Root)
meanXAbin3 <- mean(Xdatabin3$Dist2Root)/mean(Adatabin3$Dist2Root)
meanXAbin4 <- mean(Xdatabin4$Dist2Root)/mean(Adatabin4$Dist2Root)
meanXAbin5 <- mean(Xdatabin5$Dist2Root)/mean(Adatabin5$Dist2Root)
meanXAbin6 <- mean(Xdatabin6$Dist2Root)/mean(Adatabin6$Dist2Root)

meanXAcorrbin1 <- (mean(Xdatabin1$Dist2Root)-mean(Xdatabin1$pi))/(mean(Adatabin1$Dist2Root)-mean(Adatabin1$pi))
meanXAcorrbin2 <- (mean(Xdatabin2$Dist2Root)-mean(Xdatabin2$pi))/(mean(Adatabin2$Dist2Root)-mean(Adatabin2$pi))
meanXAcorrbin3 <- (mean(Xdatabin3$Dist2Root)-mean(Xdatabin3$pi))/(mean(Adatabin3$Dist2Root)-mean(Adatabin3$pi))
meanXAcorrbin4 <- (mean(Xdatabin4$Dist2Root)-mean(Xdatabin4$pi))/(mean(Adatabin4$Dist2Root)-mean(Adatabin4$pi))
meanXAcorrbin5 <- (mean(Xdatabin5$Dist2Root)-mean(Xdatabin5$pi))/(mean(Adatabin5$Dist2Root)-mean(Adatabin5$pi))
meanXAcorrbin6 <- (mean(Xdatabin6$Dist2Root)-mean(Xdatabin6$pi))/(mean(Adatabin6$Dist2Root)-mean(Adatabin6$pi))

meanAlphabin1 <- ((4 - 3*meanXAbin1)/(3*meanXAbin1 - 2))
meanAlphabin2 <- ((4 - 3*meanXAbin2)/(3*meanXAbin2 - 2))
meanAlphabin3 <- ((4 - 3*meanXAbin3)/(3*meanXAbin3 - 2))
meanAlphabin4 <- ((4 - 3*meanXAbin4)/(3*meanXAbin4 - 2))
meanAlphabin5 <- ((4 - 3*meanXAbin5)/(3*meanXAbin5 - 2))
meanAlphabin6 <- ((4 - 3*meanXAbin6)/(3*meanXAbin6 - 2))


meanAlphacorrbin1 <- ((4 - 3*meanXAcorrbin1)/(3*meanXAcorrbin1 - 2))
meanAlphacorrbin2 <- ((4 - 3*meanXAcorrbin2)/(3*meanXAcorrbin2 - 2))
meanAlphacorrbin3 <- ((4 - 3*meanXAcorrbin3)/(3*meanXAcorrbin3 - 2))
meanAlphacorrbin4 <- ((4 - 3*meanXAcorrbin4)/(3*meanXAcorrbin4 - 2))
meanAlphacorrbin5 <- ((4 - 3*meanXAcorrbin5)/(3*meanXAcorrbin5 - 2))
meanAlphacorrbin6 <- ((4 - 3*meanXAcorrbin6)/(3*meanXAcorrbin6 - 2))

#------overall mean values---------------
meanXDiverTotal <- mean(selectedX$Dist2Root)
meanADiverTotal <- mean(selectedA$Dist2Root)
meanXATotal <- meanXDiverTotal/meanADiverTotal
meanXAcorrTotal <- (meanXDiverTotal-mean(selectedX$pi))/(meanADiverTotal-mean(selectedA$pi))
meanAlphaTotal <- ((4 - 3*meanXATotal)/(3*meanXATotal - 2))
meanAlphacorrTotal <- ((4 - 3*meanXAcorrTotal)/(3*meanXAcorrTotal - 2))

##----Mean Diversity and Divergence calculation of each bin---------------------
meanXDiver <- c(format(mean(Xdatabin1$Dist2Root)), format(mean(Xdatabin2$Dist2Root)), format(mean(Xdatabin3$Dist2Root)), format(mean(Xdatabin4$Dist2Root)), format(mean(Xdatabin5$Dist2Root)), format(mean(Xdatabin6$Dist2Root)), format(meanXDiverTotal))
meanADiver <- c(format(mean(Adatabin1$Dist2Root)), format(mean(Adatabin2$Dist2Root)), format(mean(Adatabin3$Dist2Root)), format(mean(Adatabin4$Dist2Root)), format(mean(Adatabin5$Dist2Root)), format(mean(Adatabin6$Dist2Root)), format(meanADiverTotal))

#*****************CONFIDENCE INTERVALS CALCULATION******************************
#----read the files in directories and read the data in each din
Xdatabin1 <- NULL
Xdatabin2 <- NULL
Xdatabin3 <- NULL
Xdatabin4 <- NULL
Xdatabin5 <- NULL
Xdatabin6 <- NULL
Adatabin1 <- NULL
Adatabin2 <- NULL
Adatabin3 <- NULL
Adatabin4 <- NULL
Adatabin5 <- NULL
Adatabin6 <- NULL

XAdatabin1 <- NULL
XAdatabin2 <- NULL
XAdatabin3 <- NULL
XAdatabin4 <- NULL
XAdatabin5 <- NULL
XAdatabin6 <- NULL
XATotal <- NULL

XAcorrbin1 <- NULL
XAcorrbin2 <- NULL
XAcorrbin3 <- NULL
XAcorrbin4 <- NULL
XAcorrbin5 <- NULL
XAcorrbin6 <- NULL
XA_Total_corr <- NULL

Xbin1D <- NULL
Xbin2D <- NULL
Xbin3D <- NULL
Xbin4D <- NULL
Xbin5D <- NULL
Xbin6D <- NULL
XTotalDiv <- NULL

Abin1D <- NULL
Abin2D <- NULL
Abin3D <- NULL
Abin4D <- NULL
Abin5D <- NULL
Abin6D <- NULL
ATotalDiv <- NULL

for (i in 1:1000){
    Xfilename <- paste(writefile,"Xfiles/Xfile",i,sep='')
    
    Afilename <- paste(writefile,"Afiles/Afile",i,sep='')
    dataX <- read.table(Xfilename,header=TRUE,sep="\t")
    dataA <- read.table(Afilename,header=TRUE,sep="\t")
    
    TotalX <- subset(subset(dataX, Dist2Root > 0), NumBases > 5000)
    TotalA <- subset(subset(dataA, Dist2Root > 0), NumBases > 5000)
    
    Xdatabin1 <- subset(subset(subset(dataX,Dist2Root > 0), NumBases > 5000), cm2genes >= 0 & cm2genes < 0.05)
    Xdatabin2 <- subset(subset(subset(dataX,Dist2Root > 0), NumBases > 5000), cm2genes >= 0.05 & cm2genes < 0.1)
    Xdatabin3 <- subset(subset(subset(dataX,Dist2Root > 0), NumBases > 5000), cm2genes >= 0.1 & cm2genes < 0.2)
    Xdatabin4 <- subset(subset(subset(dataX,Dist2Root > 0), NumBases > 5000), cm2genes >= 0.2 & cm2genes < 0.4)
    Xdatabin5 <- subset(subset(subset(dataX,Dist2Root > 0), NumBases > 5000), cm2genes >= 0.4 & cm2genes < 0.8)
    Xdatabin6 <- subset(subset(subset(dataX,Dist2Root > 0), NumBases > 5000), cm2genes >= 0.8 & cm2genes < 2.0)
    Adatabin1 <- subset(subset(subset(dataA,Dist2Root > 0), NumBases > 5000), cm2genes >= 0 & cm2genes < 0.05)
    Adatabin2 <- subset(subset(subset(dataA,Dist2Root > 0), NumBases > 5000), cm2genes >= 0.05 & cm2genes < 0.1)
    Adatabin3 <- subset(subset(subset(dataA,Dist2Root > 0), NumBases > 5000), cm2genes >= 0.1 & cm2genes < 0.2)
    Adatabin4 <- subset(subset(subset(dataA,Dist2Root > 0), NumBases > 5000), cm2genes >= 0.2 & cm2genes < 0.4)
    Adatabin5 <- subset(subset(subset(dataA,Dist2Root > 0), NumBases > 5000), cm2genes >= 0.4 & cm2genes < 0.8)
    Adatabin6 <- subset(subset(subset(dataA,Dist2Root > 0), NumBases > 5000), cm2genes >= 0.8 & cm2genes < 2.0)
    
    XAdatabin1 <- rbind(XAdatabin1,mean(Xdatabin1$Dist2Root)/mean(Adatabin1$Dist2Root))
    XAdatabin2 <- rbind(XAdatabin2,mean(Xdatabin2$Dist2Root)/mean(Adatabin2$Dist2Root))
    XAdatabin3 <- rbind(XAdatabin3,mean(Xdatabin3$Dist2Root)/mean(Adatabin3$Dist2Root))
    XAdatabin4 <- rbind(XAdatabin4,mean(Xdatabin4$Dist2Root)/mean(Adatabin4$Dist2Root))
    XAdatabin5 <- rbind(XAdatabin5,mean(Xdatabin5$Dist2Root)/mean(Adatabin5$Dist2Root))
    XAdatabin6 <- rbind(XAdatabin6,mean(Xdatabin6$Dist2Root)/mean(Adatabin6$Dist2Root))
    XATotal <- rbind(XATotal,mean(TotalX$Dist2Root)/mean(TotalA$Dist2Root))
    
    XAcorrbin1 <- rbind(XAcorrbin1,(mean(Xdatabin1$Dist2Root)-mean(Xdatabin1$pi))/(mean(Adatabin1$Dist2Root)-mean(Adatabin1$pi)))
    XAcorrbin2 <- rbind(XAcorrbin2,(mean(Xdatabin2$Dist2Root)-mean(Xdatabin2$pi))/(mean(Adatabin2$Dist2Root)-mean(Adatabin2$pi)))
    XAcorrbin3 <- rbind(XAcorrbin3,(mean(Xdatabin3$Dist2Root)-mean(Xdatabin3$pi))/(mean(Adatabin3$Dist2Root)-mean(Adatabin3$pi)))
    XAcorrbin4 <- rbind(XAcorrbin4,(mean(Xdatabin4$Dist2Root)-mean(Xdatabin4$pi))/(mean(Adatabin4$Dist2Root)-mean(Adatabin4$pi)))
    XAcorrbin5 <- rbind(XAcorrbin5,(mean(Xdatabin5$Dist2Root)-mean(Xdatabin5$pi))/(mean(Adatabin5$Dist2Root)-mean(Adatabin5$pi)))
    XAcorrbin6 <- rbind(XAcorrbin6,(mean(Xdatabin6$Dist2Root)-mean(Xdatabin6$pi))/(mean(Adatabin6$Dist2Root)-mean(Adatabin6$pi)))
    XA_Total_corr <- rbind(XA_Total_corr,(mean(TotalX$Dist2Root)-mean(TotalX$pi))/(mean(TotalA$Dist2Root)-mean(TotalA$pi)))
    
    Xbin1D <- rbind(Xbin1D,mean(Xdatabin1$Dist2Root))
    Xbin2D <- rbind(Xbin2D,mean(Xdatabin2$Dist2Root))
    Xbin3D <- rbind(Xbin3D,mean(Xdatabin3$Dist2Root))
    Xbin4D <- rbind(Xbin4D,mean(Xdatabin4$Dist2Root))
    Xbin5D <- rbind(Xbin5D,mean(Xdatabin5$Dist2Root))
    Xbin6D <- rbind(Xbin6D,mean(Xdatabin6$Dist2Root))
    XTotalDiv <- rbind(XTotalDiv, mean(TotalX$Dist2Root))
    
    Abin1D <- rbind(Abin1D,mean(Adatabin1$Dist2Root))
    Abin2D <- rbind(Abin2D,mean(Adatabin2$Dist2Root))
    Abin3D <- rbind(Abin3D,mean(Adatabin3$Dist2Root))
    Abin4D <- rbind(Abin4D,mean(Adatabin4$Dist2Root))
    Abin5D <- rbind(Abin5D,mean(Adatabin5$Dist2Root))
    Abin6D <- rbind(Abin6D,mean(Adatabin6$Dist2Root))
    ATotalDiv <- rbind(ATotalDiv, mean(TotalA$Dist2Root))
}

# Calculate CI of Divergence using quantile function
bin1_X_D <- quantile(Xbin1D, c(0.025,0.975))
bin2_X_D <- quantile(Xbin2D, c(0.025,0.975))
bin3_X_D <- quantile(Xbin3D, c(0.025,0.975))
bin4_X_D <- quantile(Xbin4D, c(0.025,0.975))
bin5_X_D <- quantile(Xbin5D, c(0.025,0.975))
bin6_X_D <- quantile(Xbin6D, c(0.025,0.975))
bin1_A_D <- quantile(Abin1D, c(0.025,0.975))
bin2_A_D <- quantile(Abin2D, c(0.025,0.975))
bin3_A_D <- quantile(Abin3D, c(0.025,0.975))
bin4_A_D <- quantile(Abin4D, c(0.025,0.975))
bin5_A_D <- quantile(Abin5D, c(0.025,0.975))
bin6_A_D <- quantile(Abin6D, c(0.025,0.975))
TotalDivX_CI <- quantile(XTotalDiv, c(0.025,0.975))
TotalDivA_CI <- quantile(ATotalDiv, c(0.025,0.975))

# Calculate CI of X/A using quantile function
bin1_XA_CI <- quantile(XAdatabin1, c(0.025,0.975))
bin2_XA_CI <- quantile(XAdatabin2, c(0.025,0.975))
bin3_XA_CI <- quantile(XAdatabin3, c(0.025,0.975))
bin4_XA_CI <- quantile(XAdatabin4, c(0.025,0.975))
bin5_XA_CI <- quantile(XAdatabin5, c(0.025,0.975))
bin6_XA_CI <- quantile(XAdatabin6, c(0.025,0.975))
XA_Total_CI <- quantile(XATotal, c(0.025,0.975))

#Calculate CI of Corrected X/A using quantile function
bin1_XA_corr_CI <- quantile(XAcorrbin1, c(0.025,0.975))
bin2_XA_corr_CI <- quantile(XAcorrbin2, c(0.025,0.975))
bin3_XA_corr_CI <- quantile(XAcorrbin3, c(0.025,0.975))
bin4_XA_corr_CI <- quantile(XAcorrbin4, c(0.025,0.975))
bin5_XA_corr_CI <- quantile(XAcorrbin5, c(0.025,0.975))
bin6_XA_corr_CI <- quantile(XAcorrbin6, c(0.025,0.975))
XA_Total_corr_CI <- quantile(XA_Total_corr, c(0.025,0.975))

#--------convert the X/A CI to Alpha CI ratio to alpha
bin1_Alpha_lower <- ((4 - 3*bin1_XA_CI[2])/(3*bin1_XA_CI[2] - 2))
bin1_Alpha_upper <- ((4 - 3*bin1_XA_CI[1])/(3*bin1_XA_CI[1] - 2))
bin2_Alpha_lower <- ((4 - 3*bin2_XA_CI[2])/(3*bin2_XA_CI[2] - 2))
bin2_Alpha_upper <- ((4 - 3*bin2_XA_CI[1])/(3*bin2_XA_CI[1] - 2))
bin3_Alpha_lower <- ((4 - 3*bin3_XA_CI[2])/(3*bin3_XA_CI[2] - 2))
bin3_Alpha_upper <- ((4 - 3*bin3_XA_CI[1])/(3*bin3_XA_CI[1] - 2))
bin4_Alpha_lower <- ((4 - 3*bin4_XA_CI[2])/(3*bin4_XA_CI[2] - 2))
bin4_Alpha_upper <- ((4 - 3*bin4_XA_CI[1])/(3*bin4_XA_CI[1] - 2))
bin5_Alpha_lower <- ((4 - 3*bin5_XA_CI[2])/(3*bin5_XA_CI[2] - 2))
bin5_Alpha_upper <- ((4 - 3*bin5_XA_CI[1])/(3*bin5_XA_CI[1] - 2))
bin6_Alpha_lower <- ((4 - 3*bin6_XA_CI[2])/(3*bin6_XA_CI[2] - 2))
bin6_Alpha_upper <- ((4 - 3*bin6_XA_CI[1])/(3*bin6_XA_CI[1] - 2))
Total_Alpha_lower <- ((4 - 3*XA_Total_CI[2])/(3*XA_Total_CI[2] - 2))
Total_Alpha_upper <- ((4 - 3*XA_Total_CI[1])/(3*XA_Total_CI[1] - 2))

#--------convert the Corrected X/A CI to Corrected Alpha CI ratio
bin1_Alpha_corr_lower <- ((4 - 3*bin1_XA_corr_CI[2])/(3*bin1_XA_corr_CI[2] - 2))
bin1_Alpha_corr_upper <- ((4 - 3*bin1_XA_corr_CI[1])/(3*bin1_XA_corr_CI[1] - 2))
bin2_Alpha_corr_lower <- ((4 - 3*bin2_XA_corr_CI[2])/(3*bin2_XA_corr_CI[2] - 2))
bin2_Alpha_corr_upper <- ((4 - 3*bin2_XA_corr_CI[1])/(3*bin2_XA_corr_CI[1] - 2))
bin3_Alpha_corr_lower <- ((4 - 3*bin3_XA_corr_CI[2])/(3*bin3_XA_corr_CI[2] - 2))
bin3_Alpha_corr_upper <- ((4 - 3*bin3_XA_corr_CI[1])/(3*bin3_XA_corr_CI[1] - 2))
bin4_Alpha_corr_lower <- ((4 - 3*bin4_XA_corr_CI[2])/(3*bin4_XA_corr_CI[2] - 2))
bin4_Alpha_corr_upper <- ((4 - 3*bin4_XA_corr_CI[1])/(3*bin4_XA_corr_CI[1] - 2))
bin5_Alpha_corr_lower <- ((4 - 3*bin5_XA_corr_CI[2])/(3*bin5_XA_corr_CI[2] - 2))
bin5_Alpha_corr_upper <- ((4 - 3*bin5_XA_corr_CI[1])/(3*bin5_XA_corr_CI[1] - 2))
bin6_Alpha_corr_lower <- ((4 - 3*bin6_XA_corr_CI[2])/(3*bin6_XA_corr_CI[2] - 2))
bin6_Alpha_corr_upper <- ((4 - 3*bin6_XA_corr_CI[1])/(3*bin6_XA_corr_CI[1] - 2))
Total_Alpha_corr_lower <- ((4 - 3*XA_Total_corr_CI[2])/(3*XA_Total_corr_CI[2] - 2))
Total_Alpha_corr_upper <- ((4 - 3*XA_Total_corr_CI[1])/(3*XA_Total_corr_CI[1] - 2))

#**************** WRITE OUTPUT FILES ******************************
#------output file for Divergence (X and Autosomes) and their CIs-------------
Diver_outfile <- paste(writefile,"_Divergence",sep='')
Bins <- c("Bin1", "Bin2", "Bin3", "Bin4", "Bin5", "Bin6", "Overall")
XDiverCI_lower <- c(format(bin1_X_D[1]), format(bin2_X_D[1]), format(bin3_X_D[1]), format(bin4_X_D[1]), format(bin5_X_D[1]), format(bin6_X_D[1]), format(TotalDivX_CI[1]))
XDiverCI_upper <- c(format(bin1_X_D[2]), format(bin2_X_D[2]), format(bin3_X_D[2]), format(bin4_X_D[2]), format(bin5_X_D[2]), format(bin6_X_D[2]), format(TotalDivX_CI[2]))
ADiverCI_lower <- c(format(bin1_A_D[1]), format(bin2_A_D[1]), format(bin3_A_D[1]), format(bin4_A_D[1]), format(bin5_A_D[1]), format(bin6_A_D[1]), format(TotalDivA_CI[1]))
ADiverCI_upper <- c(format(bin1_A_D[2]), format(bin2_A_D[2]), format(bin3_A_D[2]), format(bin4_A_D[2]), format(bin5_A_D[2]), format(bin6_A_D[2]), format(TotalDivA_CI[2]))
XA <- factor(c("XDiv", "XDiv", "XDiv", "XDiv", "XDiv", "XDiv", "XDiv" ))
File1 <- data.frame(Bins, XA, meanXDiver, XDiverCI_lower, XDiverCI_upper)
write.table(File1,Diver_outfile,quote = FALSE, row.names = FALSE,sep='\t')
XA <- factor(c("ADiv", "ADiv", "ADiv", "ADiv", "ADiv" ,"ADiv", "ADiv"))
File2 <- data.frame(Bins, XA, meanADiver, ADiverCI_lower, ADiverCI_upper)
write.table(File2,Diver_outfile,quote = FALSE, append =TRUE, row.names = FALSE,sep='\t', col.names= FALSE)

#------output file for X/A ratio (uncorrected and corrected) and their CIs------
XA_outfile <- paste(writefile,"_XA",sep='')
meanXA <- c(format(meanXAbin1), format(meanXAbin2), format(meanXAbin3), format(meanXAbin4), format(meanXAbin5), format(meanXAbin6), format(meanXATotal))
XA_CI_lower <- c(format(bin1_XA_CI[1]), format(bin2_XA_CI[1]), format(bin3_XA_CI[1]), format(bin4_XA_CI[1]), format(bin5_XA_CI[1]), format(bin6_XA_CI[1]), format(XA_Total_CI[1]))
XA_CI_upper <- c(format(bin1_XA_CI[2]), format(bin2_XA_CI[2]), format(bin3_XA_CI[2]), format(bin4_XA_CI[2]), format(bin5_XA_CI[2]), format(bin6_XA_CI[2]), format(XA_Total_CI[2]))
meanXAcorr <- c(format(meanXAcorrbin1), format(meanXAcorrbin2), format(meanXAcorrbin3), format(meanXAcorrbin4), format(meanXAcorrbin5), format(meanXAcorrbin6), format(meanXAcorrTotal))
XAcorr_CI_lower <- c(format(bin1_XA_corr_CI[1]), format(bin2_XA_corr_CI[1]), format(bin3_XA_corr_CI[1]), format(bin4_XA_corr_CI[1]), format(bin5_XA_corr_CI[1]), format(bin6_XA_corr_CI[1]), format(XA_Total_corr_CI[1]))
XAcorr_CI_upper <- c(format(bin1_XA_corr_CI[2]), format(bin2_XA_corr_CI[2]), format(bin3_XA_corr_CI[2]), format(bin4_XA_corr_CI[2]), format(bin5_XA_corr_CI[2]), format(bin6_XA_corr_CI[2]), format(XA_Total_corr_CI[2]))
Type <- factor(c("uncorr", "uncorr", "uncorr", "uncorr", "uncorr", "uncorr", "uncorr" ))
File3 <- data.frame(Bins, Type, meanXA, XA_CI_lower, XA_CI_upper)
write.table(File3,XA_outfile,quote = FALSE, row.names = FALSE,sep='\t')
Type <- factor(c("corr", "corr", "corr", "corr", "corr", "corr", "corr"))
File4 <- data.frame(Bins, Type, meanXAcorr, XAcorr_CI_lower, XAcorr_CI_upper)
write.table(File4,XA_outfile,quote = FALSE, append =TRUE, row.names = FALSE,sep='\t', col.names= FALSE)

#------output file for Alpha values (uncorrected and corrected) and their CIs
Alpha_outfile <- paste(writefile,"_Alpha",sep='')
meanAlpha <- c(format(meanAlphabin1), format(meanAlphabin2), format(meanAlphabin3), format(meanAlphabin4), format(meanAlphabin5), format(meanAlphabin6), format(meanAlphaTotal))
Alpha_lower <- c(format(bin1_Alpha_lower), format(bin2_Alpha_lower), format(bin3_Alpha_lower), format(bin4_Alpha_lower), format(bin5_Alpha_lower), format(bin6_Alpha_lower), format(Total_Alpha_lower))
Alpha_upper <- c(format(bin1_Alpha_upper), format(bin2_Alpha_upper), format(bin3_Alpha_upper), format(bin4_Alpha_upper), format(bin5_Alpha_upper), format(bin6_Alpha_upper), format(Total_Alpha_upper))
meanAlphacorr <- c(format(meanAlphacorrbin1), format(meanAlphacorrbin2), format(meanAlphacorrbin3), format(meanAlphacorrbin4), format(meanAlphacorrbin5), format(meanAlphacorrbin6), format(meanAlphacorrTotal))
Alpha_corr_lower <- c(format(bin1_Alpha_corr_lower), format(bin2_Alpha_corr_lower), format(bin3_Alpha_corr_lower), format(bin4_Alpha_corr_lower), format(bin5_Alpha_corr_lower), format(bin6_Alpha_corr_lower), format(Total_Alpha_corr_lower))
Alpha_corr_upper <- c(format(bin1_Alpha_corr_upper), format(bin2_Alpha_corr_upper), format(bin3_Alpha_corr_upper), format(bin4_Alpha_corr_upper), format(bin5_Alpha_corr_upper), format(bin6_Alpha_corr_upper), format(Total_Alpha_corr_upper))
Type <- factor(c("uncorr", "uncorr", "uncorr", "uncorr", "uncorr", "uncorr", "uncorr"))
File5 <- data.frame(Bins, Type, meanAlpha, Alpha_lower, Alpha_upper )
write.table(File5,Alpha_outfile,quote = FALSE, row.names = FALSE,sep='\t')
Type <- factor(c("corr", "corr", "corr", "corr", "corr", "corr","corr" ))
File6 <- data.frame(Bins, Type, meanAlphacorr, Alpha_corr_lower, Alpha_corr_upper)
write.table(File6,Alpha_outfile,quote = FALSE, append =TRUE, row.names = FALSE,sep='\t', col.names= FALSE)
