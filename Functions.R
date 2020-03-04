 
  #Cohorts:
  IMT=grep("IMTCDIF\\-P",Integrated$Sample.ID) #c(1:21)
  Nodonor=-grep("D00",Integrated$Sample.ID)# -c(22:26)
  BioFire=c(grep("INV",Integrated$Sample.ID),grep("BF",Integrated$Sample.ID)) #c(27:49)
  Normal=c(grep("AUT",Integrated$Sample.ID,value=F),grep("PSYCH",Integrated$Sample.ID,value=F)) #c(50:57)
  DidNotResolve=which((Integrated$Resolution.of.symptoms=="no")=="TRUE")
  Resolved=which((Integrated$Resolution.of.symptoms=="yes")=="TRUE")
  Everything=c(1:58)
  AAD=grep("AAD",Integrated$Sample.ID)
  DonorThree<-which(Integrated$Donor=="D003")
  DonorTwo<-which(Integrated$Donor=="D002")
  AUT<-grep("AUTGI",Integrated$Sample.ID)
  PSYCH<-grep("PSYCH",Integrated$Sample.ID)
  Attribs<-list()
Attribs[[1]]<-Nodonor
Attribs[[2]]<-IMT
Attribs[[3]]<-BioFire
Attribs[[4]]<-Normal

#SANITY CHECK: NO FALSES, ONLY TRUES
which(unlist(sapply(seq(1,58), function(x) (SUPERFOCUS_outputbinning2[,x]$Read.Name)==(SUPERFOCUS_outputbinning2[,x]$Read.Name)))=="FALSE")

#This code involves what internal notes call WGS_BGI_2.
#Analysis involves metagenome data from which individual de novo assembly via megahit was carried out, and contig clustering performed via CD-HIT; taxonomy via Kraken using the MaxiKraken database; and annotation using SUPER-FOCUS.


#Here the main input is the cluster number, which when combined with CD-HIT input data and SUPER-FOCUS, outputs all the annotations for all contigs across all samples that were in this cluster.


CDHitToSuperFocus_original<-function(ClusterNum)
{
  return((((sapply(CDHitClusterReps3$ContigName[which(CDHitClusterReps3$ClusterNum==as.numeric(ClusterNum))], function(Temp) {
    A=gsub("^.*_k","k",Temp)
    B=gsub("_k.*$","",Temp)
    Col=grep(B,Integrated$Sample.ID)
    D=grep(A,SUPERFOCUS_outputbinning[,Col]$Read.Name)
    SUPERFOCUS_outputbinning[,Col]$Function[D]
  })))))
  
  
}

#Here the main input is the cluster number, which when combined with CD-HIT input data and SUPER-FOCUS, outputs all the annotations for all contigs across all samples that were in this cluster.
#Alternative version of above. Both are buggy and still need fixes. 

CDHitToSuperFocus<-function(CDHITInput,SUPERFOCUS_Input,ClusterNum)
{
  Vector=CDHITInput$ContigName[which(CDHITInput$ClusterNum==as.numeric(ClusterNum))]
  
  return(sapply(Vector, function(w)
  {
    #print(w)
    #print(gsub(".*_k","k",w))
    #print(gsub("_k.*","",w))
    #print(which(gsub("_k.*","",w)==Integrated$Sample.ID))
    #print(SUPERFOCUS_Input[,which(gsub("_k.*","",w)==Integrated$Sample.ID)]$Read.Name[1:10])
    #print(SUPERFOCUS_Input[,which(gsub("_k.*","",w)==Integrated$Sample.ID)]$Function[which(gsub(".*_k","k",w)==SUPERFOCUS_Input[,which(gsub("_k.*","",w)==Integrated$Sample.ID)]$Read.Name)])
    SUPERFOCUS_Input[,which(gsub("_k.*","",w)==Integrated$Sample.ID)]$Function[which(gsub(".*_k","k",w)==SUPERFOCUS_Input[,which(gsub("_k.*","",w)==Integrated$Sample.ID)]$Read.Name)]
  }
  ))
  
}




# tail(subset(QAZ[,59][(order(QAZ[,59]))],is.na(QAZ[,59][(order(QAZ[,59]))])=="FALSE"))


#This returns taxonomy of a given CDHIT cluster, by taking a CDHIT cluster number, selecting representative contigs and interrogating a pre-initiated Kraken matrix.


CDHitToKraken<-function(CDHitCluster3Input,ClusterNumber,KrakenInput)
{
  return(KrakenInput$TaxID[match(CDHitCluster3Input$ContigName[which(CDHitCluster3Input$ClusterNum==ClusterNumber)],KrakenInput$Header)])
}


#This function looks SUPERFOCUS annotations of a given L3 pathway, example, "Glycerol,... etc", infers all the child L4 unique members of this L3 pathway, and displays only those rows in a heatmap of all L4's. 

SUPERFOCUS_L3LookupLall<-function(ColOfInterest,InputMatrix,Integrated,Name)
{
  library(pheatmap)
  SUPERFOCUS_Lall<-sapply(seq(1,58), function(x) { read.csv(paste("/Users/cseto/Documents/CDiffMeta/Superfocus_diamond/",Integrated$Sample.ID[x],"_f13_megahit_superfocus/output_all_levels_and_function.xls",sep=''),sep='\t',header=T,skip=4,stringsAsFactors = F,strip.white=T)})
  
  Row=match(unique(unlist(sapply(seq(1,58), function(x) as.character(SUPERFOCUS_Lall[,x]$Function[which(SUPERFOCUS_Lall[,x]$Subsystem.Level.3==Name)])))),row.names(SUPERFOCUS_Lall_matrix_prime_prime..))
  
  pheatmap(InputMatrix[Row,],labels_col = paste(Integrated$CDI.risk.index_QW,Integrated$Sample.ID),angle_col = 315,labels_row=paste(sapply(Row,function(w) cor(as.numeric(matrix(Integrated$CDI.risk.index_QW))[ColOfInterest],InputMatrix[w,ColOfInterest])), row.names(SUPERFOCUS_Lall_matrix_prime_prime..)[Row]))
  
  
}

#The generalized version of above, displaying the children of the L level above. 

SUPERFOCUS_LevelSubLookup<-function(ColOfInterest,InputMatrix,InputMatrix2,Name,UpperLevel,Integrated)
{
  library(pheatmap)
  
  LowerLevel=UpperLevel+1;
  #Row=match(unique(unlist(sapply(seq(1,58), function(x) as.character(SUPERFOCUS_Lall[,x]$Function[which(SUPERFOCUS_Lall[,x]$Subsystem.Level.3==Name)])))),row.names(SUPERFOCUS_Lall_matrix_prime_prime..))
  
  SUPERFOCUS_Lall<-sapply(seq(1,58), function(x) { read.csv(paste("/Users/cseto/Documents/CDiffMeta/Superfocus_diamond/",Integrated$Sample.ID[x],"_f13_megahit_superfocus/output_all_levels_and_function.xls",sep=''),sep='\t',header=T,skip=4)})
  #Name<-row.names(SUPERFOCUS_L3_matrix_prime_prime..)[order((CorrVec))[is.na((CorrVec[order((CorrVec))]))=="FALSE"][(WhereNABegins-1)]];
  #Name<-row.names(SUPERFOCUS_L3_matrix_prime_prime..)[order((CorrVec))[is.na((CorrVec[order((CorrVec))]))=="FALSE"][(WhereNABegins-1)]];
#pheatmap(SUPERFOCUS_Lall_matrix_prime_prime..[Row,],labels_col = paste(Integrated$CDI.risk.index_QW,Integrated$Sample.ID),angle_col = 315,labels_row=)
  #Row=match(unique(unlist(sapply(seq(1,58), function(x) as.character(SUPERFOCUS_Lall[,x]$Function[which(SUPERFOCUS_Lall[,x]$Subsystem.Level.3==Name)])))),row.names(SUPERFOCUS_Lall_matrix_prime_prime..))
  Row=match(unique(unlist(sapply(seq(1,58), function(x) as.character(SUPERFOCUS_Lall[,x][[LowerLevel]][which(SUPERFOCUS_Lall[,x][[UpperLevel]]==Name)])))),row.names(InputMatrix2))
  #Row=match(unique(unlist(sapply(seq(1,58), function(x) as.character(SUPERFOCUS_Lall[,x]$Function[which(SUPERFOCUS_Lall[,x]$Subsystem.Level.3==Name)])))),row.names(SUPERFOCUS_Lall_matrix_prime_prime..))
  RowLabel=paste(gsub("^1","1.000",round(sapply(Row,function(w) cor(as.numeric(Integrated$CDI.risk.index_QW[ColOfInterest]),InputMatrix2[w,ColOfInterest])),3)), row.names(InputMatrix2)[Row],sep='   ')
  #RowLabel=paste(sapply(Row,function(w) cor(as.numeric(Integrated$CDI.risk.index_QW[-c(22:26)]),SUPERFOCUS_Lall_matrix_prime_prime..[w,-c(22:26)])), row.names(SUPERFOCUS_Lall_matrix_prime_prime..)[Row])
  #RowLabel=paste(sapply(Row,function(w) cor(as.numeric(Integrated$CDI.risk.index_QW[-c(22:26)]),SUPERFOCUS_Lall_matrix_prime_prime..[w,-c(22:26)])), row.names(SUPERFOCUS_Lall_matrix_prime_prime..)[Row]))
  CorrVec=sapply(Row,function(w) cor(as.numeric(Integrated$CDI.risk.index_QW[-c(22:26)]),SUPERFOCUS_Lall_matrix_prime_prime..[w,-c(22:26)]))
  PheatmapPlotter(InputMatrix2[Row,],c(1:10),Integrated$CDI.risk.index_QW,CorrVec,10,-c(22:26))
  #pheatmap(InputMatrix2[Row,],labels_col = paste(round(Integrated$CDI.risk.index_QW,3),Integrated$Group_QW),angle_col = 315,labels_row=RowLabel)
  #pheatmap(SUPERFOCUS_Lall_matrix_prime_prime..[Row,],labels_col = paste(Integrated$CDI.risk.index_QW,Integrated$Sample.ID),angle_col = 315,labels_row=paste(sapply(Row,function(w) cor(as.numeric(Integrated$CDI.risk.index_QW[-c(22:26)]),SUPERFOCUS_Lall_matrix_prime_prime..[w,-c(22:26)])), row.names(SUPERFOCUS_Lall_matrix_prime_prime..)[Row]))
  Returnables<-list()
  Returnables$Matrix<-InputMatrix2[Row,]
  Returnables$Row<-Row
  Returnables$RowLabel<-RowLabel
  return(Returnables)
}

#library(pheatmap)

#Name<-row.names(SUPERFOCUS_L3_matrix_prime_prime..)[order((CorrVec))[is.na((CorrVec[order((CorrVec))]))=="FALSE"][(WhereNABegins-1)]];
#pheatmap(SUPERFOCUS_Lall_matrix_prime_prime..[match(unique(unlist(sapply(seq(1,58), function(x) as.character(SUPERFOCUS_Lall[,x]$Function[which(SUPERFOCUS_Lall[,x]$Subsystem.Level.3==Name)])))),row.names(SUPERFOCUS_Lall_matrix_prime_prime..)),],labels_col = paste(Integrated$CDI.risk.index_QW,Integrated$Sample.ID),angle_col = 315,labels_row=paste(sapply(match(unique(unlist(sapply(seq(1,58), function(x) as.character(SUPERFOCUS_Lall[,x]$Function[which(SUPERFOCUS_Lall[,x]$Subsystem.Level.3==Name)])))),row.names(SUPERFOCUS_Lall_matrix_prime_prime..)),function(w) cor(as.numeric(Integrated$CDI.risk.index_QW[-c(22:26)]),SUPERFOCUS_Lall_matrix_prime_prime..[w,-c(22:26)])), row.names(SUPERFOCUS_Lall_matrix_prime_prime..)[match(unique(unlist(sapply(seq(1,58), function(x) as.character(SUPERFOCUS_Lall[,x]$Function[which(SUPERFOCUS_Lall[,x]$Subsystem.Level.3==Name)])))),row.names(SUPERFOCUS_Lall_matrix_prime_prime..))]))


  Calculon<-function(InputCols,InputMatrix,CorrVec,RiskIndex,OrderFlag)
  {
    RiskIndex<-as.numeric(matrix(RiskIndex))
    if(OrderFlag=="TopFour")
    {#print("T4");
      OrderVec<-order((CorrVec))}
    else if(OrderFlag=="BotFour")
    {#print("B4");
      OrderVec<-order((CorrVec))}
    else #(OrderFlag=="Abs")
    {#print("Abs");
      OrderVec<-order(abs(CorrVec))}
    WhereNABegins<-grep("TRUE",is.na(CorrVec[OrderVec]))[1]
    if(is.na(WhereNABegins)=="TRUE")
    {#print("No NAs?")
      WhereNABegins=nrow(InputMatrix)} else
      {#print("We see a NA");
        WhereNABegins=WhereNABegins-1} #This is more like, "last real number"
    if(OrderFlag=="Abs")
    {Range<-(WhereNABegins-3):(WhereNABegins)}
    if(OrderFlag=="TopFour")
    {Range<-(WhereNABegins-3):(WhereNABegins)}
    if(OrderFlag=="BotFour")
    {Range<-c(1:4)}
    else
    {Range<-(WhereNABegins-3):(WhereNABegins)}
    SubNameVec<-row.names(InputMatrix)[OrderVec]
    SubCorrVec<-CorrVec[OrderVec]
    #print(cbind(SubNameVec[Range],SubCorrVec[Range]))
    Returnable<-list()
    Returnable$OrderVec<-OrderVec
    Returnable$SubNameVec<-SubNameVec
    Returnable$SubCorrVec<-SubCorrVec
    return(Returnable)
  }
  
  
  QuadPlotter<-function(InputCols,Input4,RiskIndex,InputCorVec4)
  {
    Youngest=which(Integrated$Age==2)
    RiskIndex<-round(as.numeric(matrix(RiskIndex)),3)
    par(mfrow=c(2,2));sapply(seq(1,4), function(x) {
      plot(as.numeric(Input4[x,]),RiskIndex, axes=TRUE,xaxt='n',yaxt='n',lwd.axis=2,ylab="",main="",xlab="");
      box(lty=1,col='black',lwd=2)
      title(main = list(row.names(Input4)[x], cex = 1,col='black',font=2))
      points(as.numeric(Input4[x,c(1:21)]),RiskIndex[c(1:21)],pch=16,col='red',cex=2)
      points(as.numeric(Input4[x,c(27:49)]),RiskIndex[c(27:49)],pch=17,col='green',cex=2)
      points(as.numeric(Input4[x,c(50:57)]),RiskIndex[c(50:57)],pch=18,col='blue',cex=2)
      #points(as.numeric(Input4[x,InputCols]),RiskIndex[InputCols],cex=2,col='purple',cex=2)
      points(as.numeric(Input4[x,Youngest]),RiskIndex[Youngest],cex=4,col='black')
      legend("right",box.lwd=2,inset=0.02,legend=c("FMT","CDI/AAD","Healthy"),col=c("red","green","blue"),pch=c(16,17,18),cex=1.5)
     # plot(rnorm(99), bty="n", axes=FALSE, xlab="", ylab="")
      axis(1, col="black", col.ticks="black", col.axis="black", cex.axis=1.5)
      axis(2, col="black", col.ticks="black", col.axis="black", cex.axis=1.5)
      mtext(paste("Correlation: ",round(InputCorVec4[x],3)), side=1, line=3, col="black", cex=1.3)
      mtext("Risk Index", side=2, line=3, col="black", cex=1.3)
    })}
  
  OrderCaller<-function(CorrVec,OrderFlag)
  {
    if(OrderFlag=="TopFour")
    {#print("T4");
      OrderVec<-order((CorrVec))}
    else if(OrderFlag=="BotFour")
    {#print("B4");
      OrderVec<-order((CorrVec))}
    else #(OrderFlag=="Abs")
    {#print("Abs");
      OrderVec<-order(abs(CorrVec))}
    return(OrderVec)
  }
  
  RangeCaller<-function(WhereNABegins,OrderFlag)
  {
  
    if(OrderFlag=="Abs")
    {Range<-(WhereNABegins-3):(WhereNABegins)}
    if(OrderFlag=="TopFour")
    {Range<-(WhereNABegins-3):(WhereNABegins)}
    if(OrderFlag=="BotFour")
    {Range<-c(1:4)}
    else
    {Range<-(WhereNABegins-3):(WhereNABegins)}  
  return(Range)  
  }
  
  QuadPlot<-function(InputCols,InputMatrix,RiskIndex,CorMeth,OrderFlag)
  {
    RiskIndex<-as.numeric(matrix(RiskIndex))
    CorrVec<-CorCalc(InputCols,InputMatrix,RiskIndex,CorMeth)
    OrderVec<-OrderCaller(CorrVec,OrderFlag)

    SubNameVec<-row.names(InputMatrix)[OrderVec]
    SubCorrVec<-CorrVec[OrderVec]
    
    WhereNABegins<-grep("TRUE",is.na(SubCorrVec))[1]
    if(is.na(WhereNABegins)=="TRUE") 
    { 
      WhereNABegins=nrow(InputMatrix)} else
      { ; 
        WhereNABegins=WhereNABegins-1} #This is more like, "last real number"
    
    
    Range<-RangeCaller(WhereNABegins,OrderFlag)
    print(cbind(SubNameVec[Range],SubCorrVec[Range]))
    Input4<-InputMatrix[OrderVec,][Range,]
    CorVec4<-CorrVec[OrderVec][Range]
    QuadPlotter(InputCols,Input4,RiskIndex,CorVec4)
    #return(InputMatrix,CorrVec,OrderVec)
  }
    

  
  NO_NAME_Fixer<-function(InputMatrix,Source)
  {
    TestRows<-row.names(InputMatrix)
    for (I in seq(1,length(grep("NO_NAME",row.names(InputMatrix))))) TestRows[grep("NO_NAME",row.names(InputMatrix))[I]]<-paste(gsub(" NO_NAME","",row.names(InputMatrix)[grep("NO_NAME",row.names(InputMatrix))[I]]),unlist(Source[grep(unlist(strsplit(row.names(InputMatrix)[grep("NO_NAME",row.names(InputMatrix))[I]],":"))[1],Source[,1]),2]))
    row.names(InputMatrix)<-TestRows
    return(InputMatrix)
  }
  
  ZeroTrimmer<-function(NumZero,InputMatrix)
  {
    RowsWithGivenNumberofZeros=which(sapply(seq(1,nrow(InputMatrix)), function(x) length(which(as.numeric(InputMatrix[x,])==0)))<=NumZero)
    return(InputMatrix[RowsWithGivenNumberofZeros,]
    )}
  
  CorTestCalc<-function(InputCols,InputMatrix,RiskIndex)
  {
    RiskIndex=as.numeric(matrix(RiskIndex))
    CorrVec<-sapply(seq(1,nrow(InputMatrix)), function(x) cor.test(as.numeric(InputMatrix[x,InputCols]),RiskIndex[InputCols])$p.value)
    #CorrVec[is.na(CorrVec)=="TRUE"]<-1
    return(CorrVec)
  }
  
  CorCalc<-function(InputCols,InputMatrix,RiskIndex,CorMeth)
  {
    RiskIndex=as.numeric(matrix(RiskIndex))
    CorrVec<-sapply(seq(1,nrow(InputMatrix)), function(x) cor(as.numeric(InputMatrix[x,InputCols]),RiskIndex[InputCols],method=CorMeth))
    #CorrVec[is.na(CorrVec)=="TRUE"]<-0
    
    return(CorrVec)
  }
  
  CorPlot<-function(Tee1,Tee2,WhichP)
  {
    NAVector<-which(is.na(Tee1)=="FALSE")
    Tee1<-Tee1[NAVector]
    Tee2<-Tee2[NAVector]
    
    plot(Tee1,(-log10(Tee2)),xlab='Pearson Correlation Value',ylab='-log10 cor.test',xlim=c(-1,1),ylim=c(0,10));
    abline(h=min(-log10(Tee2[WhichP])),lty=2)
  }
  
  
  
  CorPoints<-function(Tee1,Tee2,WhichP,Color)
  {
    NAVector<-which(is.na(Tee1)=="FALSE")
    Tee1<-Tee1[NAVector]
    Tee2<-Tee2[NAVector]
    
    points(Tee1,(-log10(Tee2)),xlab='Pearson Correlation Value',col=Color);
    
    abline(h=min(-log10(Tee2[WhichP])),lty=2)
    
  }
  
  
  CorCheck<-function(InputCols,InputMatrix,RiskIndex,CorMeth)
  {
    T1<-CorCalc(InputCols,InputMatrix,as.numeric(matrix(RiskIndex)),CorMeth);
    T2=CorTestCalc(InputCols,InputMatrix,as.numeric(matrix(RiskIndex))) 
    WhichP<-which(T2<0.05); 
    Blob<-list()
    Blob$T1<-T1
    Blob$T2<-T2
    Blob$WhichP<-WhichP
    #Blob$HVal
    return(Blob)
  }
  
  Range2Caller<-function(OrderFlag,WhereNABegins,NumRankedforHeatmap)
  {
    if(OrderFlag=="Abs")
    {
      Range2<-((WhereNABegins-NumRankedforHeatmap+1):(WhereNABegins))
    }
    if(OrderFlag=="TopFour")
    {
      Range2<-((WhereNABegins-NumRankedforHeatmap+1):(WhereNABegins))
    }
    if(OrderFlag=="BotFour")
    {
      Range2<-c(1:NumRankedforHeatmap)
    }
    else
    {
      Range2<-((WhereNABegins-NumRankedforHeatmap+1):(WhereNABegins))
    }
  return(Range2)  
  }
  
  PheatmapPlotter<-function(InputMatrix,Range2,RiskIndex,CorrVec,NumRankedforHeatmap,InputColsForHeatMap)
  {
    library(pheatmap)
    #RowLabel=round(CorrVec[Range2],3);
    RowLabel=row.names(InputMatrix)[Range2];
    RowLabel=paste(round(CorrVec[Range2],3),row.names(InputMatrix)[Range2]); #row.names(InputMatrix)[Range2] #round(CorrVec[Range2],3) #
   ColLabel=paste((gsub("^1","1.000",round(RiskIndex[InputColsForHeatMap],3))),Integrated$Group_QW[InputColsForHeatMap],sep='    ')
    #jpeg('rplot.jpg')
    pheatmap(InputMatrix[Range2,InputColsForHeatMap],labels_col=ColLabel,angle_col = 315,cellheight=15,cellwidth=15,fontsize=14,labels_row=RowLabel);
    #dev.off()
    print("END HEATMAP")
  }
  
  MatrixCorCalculator<-function(InputCols,InputMatrix,NumRankedforHeatmap,CorMeth,InputColsForHeatMap,OrderFlag)
  {
    library(pheatmap)
    RiskIndex<-as.numeric(matrix(Integrated$CDI.risk.index_QW))
    #CorrVec<-CorCalc(InputCols,InputMatrix,RiskIndex,CorMeth)
    print("Line 1")
    CorThings<-CorCheck(InputCols,InputMatrix,RiskIndex,CorMeth)
    CorrVec<-unlist(CorThings$T1)
    CorrP<-unlist(CorThings$T2)
    ReturnMe<-cbind(InputMatrix,CorrVec,CorrP)
    
    print("Line 2")
    
    OrderVec<-OrderCaller(CorrVec,OrderFlag)
    SubNameVec<-row.names(InputMatrix)[OrderVec]
    SubCorrVec<-CorrVec[OrderVec]
    
    WhereNABegins<-grep("TRUE",is.na(SubCorrVec))[1]
    if(is.na(WhereNABegins)=="TRUE") 
    { 
      WhereRealEnds=nrow(InputMatrix)} else
      { ; 
        WhereRealEnds=WhereNABegins-1} #This is more like, "last real number"
    Range<-RangeCaller(WhereRealEnds,OrderFlag)
    Input4<-InputMatrix[OrderVec,][Range,]
    CorVec4<-CorrVec[OrderVec][Range]
    print(cbind(SubNameVec[Range],SubCorrVec[Range]))
    
   return(ReturnMe) 
  }

  Plotter<-function(InputCols,InputMatrix,NumRankedforHeatmap,CorMeth,InputColsForHeatMap,OrderFlag)
  {
    library(pheatmap)
    print("Line 1")
    CorThings<-CorCheck(InputCols,InputMatrix,RiskIndex,CorMeth)
    
    CorrVec<-unlist(CorThings$T1)
    CorrP<-unlist(CorThings$T2)
    ReturnMe<-cbind(InputMatrix,CorrVec,unlist(CorThings$T2)) #MatrixCorCalculator(InputCols,InputMatrix,NumRankedforHeatmap,CorMeth,InputColsForHeatMap,OrderFlag)
      
      
    print("Line 2")
    
    OrderVec<-OrderCaller(CorrVec,OrderFlag)
    
    SubNameVec<-row.names(InputMatrix)[OrderVec]
    SubCorrVec<-CorrVec[OrderVec]
    
    WhereNABegins<-grep("TRUE",is.na(SubCorrVec))[1]
    if(is.na(WhereNABegins)=="TRUE") 
    { 
    WhereRealEnds=nrow(InputMatrix)} else
    { ; 
      WhereRealEnds=WhereNABegins-1} #This is more like, "last real number"
    
        Range<-RangeCaller(WhereRealEnds,OrderFlag)
    Input4<-InputMatrix[OrderVec,][Range,]
    CorVec4<-CorrVec[OrderVec][Range]
    
    QuadPlotter(InputCols,Input4,RiskIndex,CorVec4)
    print(cbind(SubNameVec[Range],SubCorrVec[Range]))
    Range2<-Range2Caller(OrderFlag,WhereRealEnds,NumRankedforHeatmap)
    PheatmapPlotter(InputMatrix[OrderVec,],Range2,RiskIndex,CorrVec[OrderVec],NumRankedforHeatmap,InputColsForHeatMap)
    return(ReturnMe)
    
    }
  
HeatmapGenerator<-function(InputCols,InputMatrix,OrderedCorrVec,NumRankedforHeatmap,InputColsForHeatMap,OrderFlag)  
{
  WhereNABegins<-grep("TRUE",is.na(OrderedCorrVec))[1]
  
  if(is.na(WhereNABegins)=="TRUE") 
  { #print("No NAs?")
    WhereNABegins=nrow(InputMatrix)} else
    {#print("We see a NA"); 
      WhereNABegins=WhereNABegins-1} #This is more like, "last real number"
  if(OrderFlag=="Abs")
  {Range<-(WhereNABegins-3):(WhereNABegins);
  Range2<-((WhereNABegins-NumRankedforHeatmap+1):(WhereNABegins))
  }
  if(OrderFlag=="TopFour")
  {Range<-(WhereNABegins-3):(WhereNABegins);
  Range2<-((WhereNABegins-NumRankedforHeatmap+1):(WhereNABegins))
  }
  if(OrderFlag=="BotFour")
  {Range<-c(1:4);
  Range2<-c(1:NumRankedforHeatmap)
  }
  else
  {Range<-(WhereNABegins-3):(WhereNABegins);
  Range2<-((WhereNABegins-NumRankedforHeatmap+1):(WhereNABegins))
  }
  if(OrderFlag=="Abs")
  {
    pheatmap(InputMatrix[OrderVec,][Range2,][,InputColsForHeatMap],labels_col=paste(RiskIndex[InputColsForHeatMap],Integrated$Sample.ID[InputColsForHeatMap],sep='    '),angle_col = 315,cellheight=12,labels_row=paste(CorrVec[OrderVec[Range2]],row.names(InputMatrix)[OrderVec[Range2]]))
    
  }
  if(OrderFlag=="TopFour")
  {
    pheatmap(InputMatrix[OrderVec,][Range2,][,InputColsForHeatMap],labels_col=paste(RiskIndex[InputColsForHeatMap],Integrated$Sample.ID[InputColsForHeatMap],sep='    '),angle_col = 315,cellheight=12,labels_row=paste(CorrVec[OrderVec[Range2]],row.names(InputMatrix)[OrderVec[Range2]]))
    #Range<-(WhereNABegins-3):(WhereNABegins);
    #Range2<-((WhereNABegins-NumRankedforHeatmap+1):(WhereNABegins))
  }
  if(OrderFlag=="BotFour")
  {
    pheatmap(InputMatrix[OrderVec,][1:NumRankedforHeatmap,][,InputColsForHeatMap],labels_col=paste(RiskIndex[InputColsForHeatMap],Integrated$Sample.ID[InputColsForHeatMap],sep='    '),angle_col = 315,cellheight=12,labels_row=paste(CorrVec[OrderVec[Range2]],row.names(InputMatrix)[OrderVec[Range2]]))
    #Range<-c(1:4);
    # Range2<-c(1:NumRankedforHeatmap)
  }
  else
  {
    pheatmap(InputMatrix[OrderVec,][Range2,][,InputColsForHeatMap],labels_col=paste(RiskIndex[InputColsForHeatMap],Integrated$Sample.ID[InputColsForHeatMap],sep='    '),angle_col = 315,cellheight=12,labels_row=paste(CorrVec[OrderVec[Range2]],row.names(InputMatrix)[OrderVec[Range2]]))
    
    # Range<-(WhereNABegins-3):(WhereNABegins);
    # Range2<-((WhereNABegins-NumRankedforHeatmap+1):(WhereNABegins))
  }  
}
  
  
  PlotterNoNA<-function(InputCols,InputMatrix,NumRankedforHeatmap,CorMeth,InputColsForHeatMap)
  {
    RiskIndex=as.numeric(matrix(Integrated$CDI.risk.index_QW))
    #ColsOfInterest=c(1:21);
    #Input<-100*Humann_GeneFamilies_combined_cpm_names_relabund
    CorrVec<-sapply(seq(1,nrow(InputMatrix)), function(x) cor(RiskIndex[InputCols],as.numeric(InputMatrix[x,InputCols]),method=CorMeth))
    OrderVec<-order(abs(CorrVec))
    
    #WhereNABegins<-grep("TRUE",is.na(CorrVec[OrderVec]))[1]
    
    Range<-(nrow(InputMatrix)-3):nrow(InputMatrix)
    SubNameVec<-row.names(InputMatrix)[OrderVec]
    print(SubNameVec[Range])
    Range2<-(nrow(InputMatrix)-NumRankedforHeatmap):(nrow(InputMatrix))
    SubCorrVec<-CorrVec[OrderVec]
    print(SubCorrVec[Range])
    print(cbind(SubNameVec[Range],SubCorrVec[Range]))
    
    par(mfrow=c(2,2));sapply(seq(1,4), function(x) {plot(as.numeric(InputMatrix[OrderVec,][(nrow(InputMatrix)-x),]),matrix(Integrated$CDI.risk.index_QW),ylab='Risk Index',xlab=row.names(InputMatrix)[OrderVec][(nrow(InputMatrix)-x)]);
      points(as.numeric(InputMatrix[OrderVec,][(nrow(InputMatrix)-x),c(1:21)]),RiskIndex[c(1:21)],pch=17,col='red')
      points(as.numeric(InputMatrix[OrderVec,][(nrow(InputMatrix)-x),c(27:49)]),RiskIndex[c(27:49)],pch=17,col='green')
      points(as.numeric(InputMatrix[OrderVec,][(nrow(InputMatrix)-x),c(50:57)]),RiskIndex[c(50:57)],pch=17,col='blue')
      points(as.numeric(InputMatrix[OrderVec,][(nrow(InputMatrix)-x),InputCols]),RiskIndex[InputCols],cex=2,col='purple')
      
    })
    
    pheatmap(InputMatrix[OrderVec,][Range2,][,InputColsForHeatMap],labels_col=paste(RiskIndex[InputColsForHeatMap],Integrated$Sample.ID[InputColsForHeatMap],sep='    '),angle_col = 315,cellheight=12,labels_row=paste(CorrVec[OrderVec[Range2]],row.names(InputMatrix)[OrderVec[Range2]]))
  }

  
  MetaphlanParser<-function(Level,MetaphlanMatrix,TaxonomyTable)
  {
    if(Level==1)
    {MainText=("Abunds of Kingdom in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if (Level==2)
    {MainText=("Abunds of Phylum in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if (Level==3)
    {MainText=("Abunds of Class in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if(Level==4)
    {MainText=("Abunds of Order in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if (Level==5)
    {MainText=("Abunds of Family in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if (Level==6)
    {MainText=("Abunds of Genus in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if (Level==7)
    {MainText=("Abunds of Species in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if (Level==8)
    {MainText=("Abunds of Strains in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,8]!="0"),]
    }
    else
    {
      MainText=("Abunds of Something in Metaphlan2")
      TempehFunction=0 #Deliberate function explosion
    }
    return(TempehFunction)  
  }
  
  AbundanceRiskPlots<-function(Level,MetaphlanMatrix,TaxonomyTable)
  {
    if(Level==1)
    {MainText=("Abunds of Kingdom in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if (Level==2)
    {MainText=("Abunds of Phylum in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if (Level==3)
    {MainText=("Abunds of Class in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if(Level==4)
    {MainText=("Abunds of Order in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if (Level==5)
    {MainText=("Abunds of Family in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if (Level==6)
    {MainText=("Abunds of Genus in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if (Level==7)
    {MainText=("Abunds of Species in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,Level]!="0" & TaxonomyTable[,(Level+1)]=="0"),]
    }
    else if (Level==8)
    {MainText=("Abunds of Strains in Metaphlan2")
    TempehFunction<-MetaphlanMatrix[which(TaxonomyTable[,8]!="0"),]
    }
    else
    {
      MainText=("Abunds of Something in Metaphlan2")
      TempehFunction=0 #Deliberate function explosion
    }
    
    plot(-10000,-10000, xlim=c(0,round(max(TempehFunction))),ylim=c(0,1), xlab='Metaphlan Abundance',ylab='Risk Index',main=MainText)
    sapply(seq(1,21), function(x) points(TempehFunction[,x],rep(as.numeric(matrix(Integrated$CDI.risk.index_QW[x])),nrow(TempehFunction)),col='red',pch=16))
    sapply(seq(27,49), function(x) points(TempehFunction[,x],rep(as.numeric(matrix(Integrated$CDI.risk.index_QW[x])),nrow(TempehFunction)),col='green',pch=16))
    sapply(seq(50,56), function(x) points(TempehFunction[,x],rep(as.numeric(matrix(Integrated$CDI.risk.index_QW[x])),nrow(TempehFunction)),col='blue',pch=16))
    sapply(which(Integrated$Resolution.of.symptoms=="no"), function(x) points(TempehFunction[,x],rep(as.numeric(matrix(Integrated$CDI.risk.index_QW[x])),nrow(TempehFunction)),col='black',pch=16))
    return(TempehFunction)
  }
  
  
  
