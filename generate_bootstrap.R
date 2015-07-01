#------------Read input file---------------------
## provide the input file to read
## If input file is not present the same directory where R is running, the provide the full path of the file
data <- read.table("hg18.farFromGenes.20kb.nosep.schweinfurthii.bothsexes.reheader",header=TRUE,sep="\t")


## provide the initial name of the directory to write the bootstrap files, the program will append Afiles (for autosomes) and Xfiles (for X-chromosomes)
## Two directories are created: (1) 1000 bootstrap files for autosomes in directory $DirectoryAfiles (2) 1000 bootstrap files for X chromosomes in directory $DirectoryXfiles
## the program will check if the directories already exist, if not it creates the new directories

Directory <- "Schweinfurthii"


#------------Filter data based on Number of bases(NumBases) and Divergence(Dist2Root) values----------------
selectedX <- subset(data,chrom == "chrX")
selectedA <- subset(data,chrom != "chrX")

#-----------Divide filtered data into windows on the bases of cM distance from genes(cm2genes)---------
##------- This step divide the data into six windows (0-0.05, 0.05-0.1, 0.1-0.2, 0.2-0.4, 0.4-0.8, 0.8-2.0)
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

#Randomly select specified number of values(locii) from each bin and write to file
for (i in 1:1000){
  XfileContent <- NULL
  AfileContent <- NULL

  selctedXbin1 <- Xdatabin1[sample(1:nrow(Xdatabin1),replace=TRUE),]
  selctedXbin2 <- Xdatabin2[sample(1:nrow(Xdatabin2),replace=TRUE),]
  selctedXbin3 <- Xdatabin3[sample(1:nrow(Xdatabin3),replace=TRUE),]
  selctedXbin4 <- Xdatabin4[sample(1:nrow(Xdatabin4),replace=TRUE),]
  selctedXbin5 <- Xdatabin5[sample(1:nrow(Xdatabin5),replace=TRUE),]
  selctedXbin6 <- Xdatabin6[sample(1:nrow(Xdatabin6),replace=TRUE),]
  XfileContent <- rbind(XfileContent,selctedXbin1, selctedXbin2, selctedXbin3, selctedXbin4, selctedXbin5, selctedXbin6)

  selctedAbin1 <- Adatabin1[sample(1:nrow(Adatabin1),replace=TRUE),]
  selctedAbin2 <- Adatabin2[sample(1:nrow(Adatabin2),replace=TRUE),]
  selctedAbin3 <- Adatabin3[sample(1:nrow(Adatabin3),replace=TRUE),]
  selctedAbin4 <- Adatabin4[sample(1:nrow(Adatabin4),replace=TRUE),]
  selctedAbin5 <- Adatabin5[sample(1:nrow(Adatabin5),replace=TRUE),]
  selctedAbin6 <- Adatabin6[sample(1:nrow(Adatabin6),replace=TRUE),]
  AfileContent <- rbind(AfileContent,selctedAbin1, selctedAbin2, selctedAbin3, selctedAbin4, selctedAbin5, selctedAbin6)

  DirX <- paste(Directory,"Xfiles",sep='')
  DirA <- paste(Directory,"Afiles",sep='')
  if (file.exists(DirX)){
    print("directory exists")
  }
  else {
    dir.create(file.path(DirX))
    print("creating directory")
  }
  if (file.exists(DirA)){
    print("directory exists")
  }
  else {
    dir.create(file.path(DirA))
    print("creating directory")
  }
  
  Xfilename <- paste(DirX,"/Xfile",i,sep='')
  Afilename <- paste(DirA,"/Afile",i,sep='')
  write.table(XfileContent,Xfilename,quote = FALSE, row.names = FALSE,sep='\t')
  write.table(AfileContent,Afilename,quote = FALSE, row.names = FALSE,sep='\t')
}
