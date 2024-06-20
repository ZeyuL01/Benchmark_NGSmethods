require(data.table)
library(cowplot)
require(gridExtra)
library(dplyr)
library(ggthemes)
require(GGally)
library(org.Hs.eg.db)
library(limma)
library(reactable)
library(reactablefmtr)

############################################################
#Set work directory
############################################################

work_dir<-"/Users/zeyulu/Dropbox/datasets/clean_results_new/"
knockTF_folder<-c("knockTF_200/","knockTF_600/","knockTF_1000/")

############################################################
#Define metrics calculation function
############################################################

MRR_cal<-function(li,K,term){
  values=which(li[1:K]==term)
  if(length(values)>=1){
    return(1/values[1])
  }else if(length(values)==0){
    return(0)
  }
}

AP_cal<-function(li,K,term){
  values<-sum((1:length(which(li[1:K]==term)))/which(li[1:K]==term))/length(which(li[1:K]==term))
  if(is.na(values)){
    values=0
  }
  return(values)
}

Hit_cal<-function(li,K,term){
  if(term%in%li[1:K]){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

NDCG_cal<-function(li,K,term){
  values<-sum(1/log2(1+which(li[1:K]==term)))/sum(1/log2(1+1:length(which(li[1:K]==term))))
  if(is.na(values)){
    values=0
  }
  return(values)
}

############################################################
Method_names<-c("BART","ChEA3","ChIP-Atlas","Cscan","Enrichr","HOMER","i-cisTarget","Lisa","MAGIC","Pscan","RcisTarget","RegulatorTrail","TFEA.ChIP")
#MRR_Table
MRR_table_200<-data.frame(matrix(nrow=13,ncol=4))
MRR_table_600<-data.frame(matrix(nrow=13,ncol=4))
MRR_table_1000<-data.frame(matrix(nrow=13,ncol=4))
colnames(MRR_table_200)<-c("Method","Top10","Top50","Top100")
colnames(MRR_table_600)<-c("Method","Top10","Top50","Top100")
colnames(MRR_table_1000)<-c("Method","Top10","Top50","Top100")
MRR_table_200$Method<-Method_names
MRR_table_600$Method<-Method_names
MRR_table_1000$Method<-Method_names
#MAP Table
MAP_table_200<-data.frame(matrix(nrow=13,ncol=4))
MAP_table_600<-data.frame(matrix(nrow=13,ncol=4))
MAP_table_1000<-data.frame(matrix(nrow=13,ncol=4))
colnames(MAP_table_200)<-c("Method","Top10","Top50","Top100")
colnames(MAP_table_600)<-c("Method","Top10","Top50","Top100")
colnames(MAP_table_1000)<-c("Method","Top10","Top50","Top100")
MAP_table_200$Method<-Method_names
MAP_table_600$Method<-Method_names
MAP_table_1000$Method<-Method_names
#Hit Rate Table
Hit_rate_table_200<-data.frame(matrix(nrow=13,ncol=4))
Hit_rate_table_600<-data.frame(matrix(nrow=13,ncol=4))
Hit_rate_table_1000<-data.frame(matrix(nrow=13,ncol=4))
colnames(Hit_rate_table_200)<-c("Method","Top10","Top50","Top100")
colnames(Hit_rate_table_600)<-c("Method","Top10","Top50","Top100")
colnames(Hit_rate_table_1000)<-c("Method","Top10","Top50","Top100")
Hit_rate_table_200$Method<-Method_names
Hit_rate_table_600$Method<-Method_names
Hit_rate_table_1000$Method<-Method_names
#NDCG
NDCG_table_200<-data.frame(matrix(nrow=13,ncol=4))
colnames(NDCG_table_200)<-c("Method","Top10","Top50","Top100")
NDCG_table_200$Method<-Method_names
NDCG_table_600<-data.frame(matrix(nrow=13,ncol=4))
colnames(NDCG_table_600)<-c("Method","Top10","Top50","Top100")
NDCG_table_600$Method<-Method_names
NDCG_table_1000<-data.frame(matrix(nrow=13,ncol=4))
colnames(NDCG_table_1000)<-c("Method","Top10","Top50","Top100")
NDCG_table_1000$Method<-Method_names

############################################################
#Fill the tables (1000), 200, 600 can be similarly modified
############################################################
work_dir<-"/Users/zeyulu/Dropbox/datasets/clean_results_new/"
knockTF_folder<-c("knockTF_200/","knockTF_600/","knockTF_1000/")
work_files<-list.files(paste0(work_dir,knockTF_folder[1]))
Method_names<-c("BART","ChEA3","ChIP-Atlas","Cscan","Enrichr","HOMER","i-cisTarget","Lisa","MAGIC","Pscan","RcisTarget","RegulatorTrail","TFEA.ChIP")

work_files
Top<-c(10,50,100)
for(k in 1:3){
  for(i in 1:length(work_files)){
    table<-read.csv(paste0(work_dir,knockTF_folder[1],work_files[i]))
    TF_names<-sapply(strsplit(colnames(table),"_",fixed=TRUE),function(x){return(x[[1]])})
    AP_list<-c()
    NDCG_list<-c()
    Hit_list<-c()
    MRR_list<-c()
    for(j in 1:ncol(table)){
      print(paste0(i,"_",j,"_",k))
      MRR_list<-c(MRR_list,MRR_cal(table[,j],Top[k],TF_names[j]))
      AP_list<-c(AP_list,AP_cal(table[,j],Top[k],TF_names[j]))
      NDCG_list<-c(NDCG_list,NDCG_cal(table[,j],Top[k],TF_names[j]))
      Hit_list<-c(Hit_list,Hit_cal(table[,j],Top[k],TF_names[j]))
    }
    Hit_rate_table_200[i,k+1]<-sum(Hit_list)/length(Hit_list)
    NDCG_table_200[i,k+1]<-mean(NDCG_list)
    MAP_table_200[i,k+1]<-mean(AP_list)
    MRR_table_200[i,k+1]<-sum(MRR_list,na.rm=TRUE)/length(MRR_list)
  }
}

############################################################
#Figure 6 bottomleft
############################################################
#par(mfrow=c(3,3),oma = c(2,2,2,0) + 0.1,mar = c(0.5,1,1,1) + 0.1,xpd=TRUE)
pdf("/Users/zeyulu/Desktop/Review_Manuscript/BiB/Figures/FIgure 5_bottom_left.pdf",width=8.9,height=8.29)
layout_matrix<-matrix(c(1:2),2,1,byrow=FALSE)
par(mfcol=c(2,1),oma = c(3,3,2,0) + 0.1,mgp=c(3,0.5,0),mar = c(5,4,4,4) + 0.1,xpd=TRUE,bty="o")
layout(layout_matrix,heights=c(1,1.32),widths=c(1))

for(k in 1:1){
  work_files<-list.files(paste0(work_dir,knockTF_folder[k]))

  method_name<-c("BART","ChEA3","ChIP-Atlas","Cscan","Enrichr","HOMER","i-cisTarget",
                 "Lisa","MAGIC","Pscan","RcisTarget","RegulatorTrail","TFEA.ChIP")

  ratio_table<-data.frame(matrix(nrow=13,ncol=3))
  colnames(ratio_table)<-c("TOP.1%","TOP.5%","TOP.10%")
  rownames(ratio_table)<-method_name

  abs_table<-data.frame(matrix(nrow=13,ncol=3))
  colnames(abs_table)<-c("TOP.10","TOP.50","TOP.100")
  rownames(abs_table)<-method_name

  na_table<-data.frame(matrix(nrow=13,ncol=1))
  colnames(na_table)<-c("NA_NUM")
  rownames(na_table)<-method_name

  for(i in 1:length(work_files)){
    table_res<-read.csv(paste0(work_dir,knockTF_folder[k],work_files[i]),header=TRUE)
    TF_index<-colnames(table_res)
    TF_names<-sapply(strsplit(TF_index,"_",fixed=TRUE),function(x){return(x[[1]])})
    rank_list<-c()
    rank_num<-c()
    for(j in 1:length(TF_index)){
      print(paste0(i,"_",j))
      TF_rank_list<-table_res[,TF_index[j]]
      rank_num<-c(rank_num,sum(!is.na(TF_rank_list)))

      rank_list<-c(rank_list,which(TF_rank_list==TF_names[j])[1])
      ratio_list<-rank_list/rank_num
    }
    top1p<-sum(ratio_list<=0.01,na.rm=TRUE)
    top5p<-sum(ratio_list<=0.05,na.rm=TRUE)
    top10p<-sum(ratio_list<=0.1,na.rm=TRUE)

    top10<-sum(rank_list<=10,na.rm=TRUE)
    top50<-sum(rank_list<=50,na.rm=TRUE)
    top100<-sum(rank_list<=100,na.rm=TRUE)

    na_val<-sum(is.na(rank_list))

    ratio_table[method_name[i],1]<-top1p
    ratio_table[method_name[i],2]<-top5p
    ratio_table[method_name[i],3]<-top10p

    abs_table[method_name[i],1]<-top10
    abs_table[method_name[i],2]<-top50
    abs_table[method_name[i],3]<-top100

    na_table[method_name[i],1]<-na_val
  }
  if(k==1){
    abs_table_200=abs_table
  }
  col_index_1<-c("BART","Lisa","ChIP-Atlas","MAGIC","TFEA.ChIP","Enrichr","RcisTarget","i-cisTarget","HOMER","Cscan","RegulatorTrail",
                 "ChEA3","Pscan")

  abs_table<-abs_table[col_index_1,]

  col_abs<-c("black","black","black","black","black","black","red","black",
             "red","black","black","black","red")

  ratio_table<-ratio_table[col_index_1,]

  col_ratio<-col_abs


  if(k==1){
    par(mar=c(0.5,2.6,1,1)+0.1)
    b<-barplot(abs_table[,"TOP.100"],col="#FFDF34",ylim=c(0,140),xaxt = "n", yaxt = "n",cex.main=1.2,tcl=-0.2)
    complement<-c(0,0,0,0,0,0,0,9,0,0,0,0,0)
    axis(2,at=c(0,30,60,90,120),labels=c(0,30,60,90,120),cex.axis=1.3,tcl=-0.2)
    text(b,abs_table[,"TOP.100"]+3+complement,abs_table[,"TOP.100"],cex=0.8,font=2)
    par(new=T)
    barplot(abs_table[,"TOP.50"],col="#FC8C5A",ylim=c(0,140),xaxt = "n", yaxt = "n",tcl=-0.2,cex.axis=1.3)
    complement<-c(0,0,0,0,0,0,0,0,0,0,0,0,0)
    text(b,abs_table[,"TOP.50"]+5+complement,abs_table[,"TOP.50"],cex=0.8,font=2,yaxt="n")
    par(new=T)
    b<-barplot(abs_table[,"TOP.10"],col="#DB3124",ylim=c(0,140),xaxt = "n", yaxt = "n",tcl=-0.2,cex.axis=1.3)
    complement<-c(0,0,0,0,0,0,0,-1,0,0,0,0,0)
    text(b,abs_table[,"TOP.10"]+5+complement,abs_table[,"TOP.10"],cex=0.8,font=2,yaxt="n")
    box()
    title(ylab="Number of Perturbed TRs",cex.lab=1.2,font.lab=1,line=1.7)
    legend("topright",legend=c("TOP 10","TOP 50","TOP 100"),col="black",pch=22,pt.bg=c("#DB3124","#FC8C5A","#FFDF34"),text.font=2)
  }
  if(k==1){
    par(mar=c(5,2.6,1,1)+0.1)
    b<-barplot(ratio_table[,"TOP.10%"],col="#98c7df",ylim=c(0,330),tcl=-0.2,xaxt = "n", yaxt = "n",cex.main=1.1)
    labs<-col_index_1
    axis(2,at=c(0,100,200,300),labels=c(0,100,200,300),cex.axis=1.3,tcl=-0.2)
    complement1<-c(0,0,0,0,0,0,0,7,0,0,0,0,0)
    text(b,ratio_table[,"TOP.10%"]+10+complement1,ratio_table[,"TOP.10%"],cex=0.7,font=2)
    par(new=T)
    barplot(ratio_table[,"TOP.5%"],col="#2c81be",ylim=c(0,330),xaxt = "n", yaxt = "n",tcl=-0.2)
    complement2<-c(0,0,0,0,0,0,0,0,0,0,0,0,0)
    text(b,ratio_table[,"TOP.5%"]+10+complement2,ratio_table[,"TOP.5%"],cex=0.7,font=2,yaxt="n")
    par(new=T)
    b<-barplot(ratio_table[,"TOP.1%"],col="#1c3e71",ylim=c(0,330),xaxt = "n", yaxt = "n",tcl=-0.2)
    complement3<-c(0,0,0,0,0,0,0,-3,0,0,0,0,0)
    text(b,ratio_table[,"TOP.1%"]+10+complement3,ratio_table[,"TOP.1%"],cex=0.7,font=2,yaxt="n")
    box()
    text(x=b,y=-60,labs,xpd=TRUE,cex=0.8,srt=90,font=2,col=col_abs)
    title(ylab="Number of Perturbed TRs",cex.lab=1.2,font.lab=1,line=1.7)
    legend("topright",legend=c("TOP 1%","TOP 5%","TOP 10%"),col="black",pch=22,pt.bg=c("#1c3e71","#2c81be","#98c7df"),text.font=2)
  }
}
title(xlab="Computational method",font.lab=1,line=-0.5,outer=TRUE,cex.main=1.2,cex.lab=1.2)
dev.off()


############################################################
#Figure 6 bottomright
############################################################

#Sort the order based on high to low
Hit_rate_table_200

col_MRR<-c("black","black","black","black","black","black","red","black",
           "red","black","black","black","red")

Method_names[c(1,8,3,9,13,5,11,7,6,4,12,2,10)]
method_name<-Method_names[c(1,8,3,9,13,5,11,7,6,4,12,2,10)]
Hit_rate_table_200<-Hit_rate_table_200[c(1,8,3,9,13,5,11,7,6,4,12,2,10),]
MRR_table_200<-MRR_table_200[c(1,8,3,9,13,5,11,7,6,4,12,2,10),]
MAP_table_200<-MAP_table_200[c(1,8,3,9,13,5,11,7,6,4,12,2,10),]
NDCG_table_200<-NDCG_table_200[c(1,8,3,9,13,5,11,7,6,4,12,2,10),]

labs=method_name
layout_matrix<-matrix(c(1:4),4,1,byrow=FALSE)
par(mfcol=c(4,1),oma = c(3,2,2,0) + 0.1,mar = c(1,5,1,1) + 0.1,xpd=TRUE,bty="o")
layout(layout_matrix,heights=c(1,1,1,1.5),widths=c(1))

plot(1:13,Hit_rate_table_200$Top10,type="b",ylim=c(0,0.2),lwd=1.5,axes=FALSE,col="#DB3124",ylab="")
complement<-c(0.02,0.02,0.02,0.02,0.02,0.02,0.02,-0.02,0.02,0.02,0.02,0.02,0.02)
text(1:13,Hit_rate_table_200$Top10+complement,round(Hit_rate_table_200$Top10,digits=3),cex=1.1,col="#DB3124",font=2)
title(ylab="Hit Rate",cex.lab=1.6)
axis(2,at=c(0,0.05,0.1,0.15,0.2),label=c("0",".05",".10",".15",".20"),las=1,cex.axis=1.4)
box()
lines(1:13,type="b",Hit_rate_table_200$Top50,col="#FC8C5A",lwd=1.5)
complement<-c(0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.01,0.02,0.02,0.02,0.02,0.02)
text(1:13,Hit_rate_table_200$Top50+complement-0.0013,round(Hit_rate_table_200$Top50,digits=3),cex=1.1,col="#DB3124",font=2)
text(1:13,Hit_rate_table_200$Top50+complement,round(Hit_rate_table_200$Top50,digits=3),cex=1.1,col="#FC8C5A",font=2)
lines(1:13,type="b",Hit_rate_table_200$Top100,col="#FFDF34",lwd=1.5)
complement<-c(0.02,0.02,0.02,0.02,-0.02,-0.02,0.02,0.04,0.02,0.02,0.02,0.02,0.02)
text(1:13,Hit_rate_table_200$Top100+complement-0.0013,round(Hit_rate_table_200$Top100,digits=3),cex=1.1,col="#DB3124",font=2)
text(1:13,Hit_rate_table_200$Top100+complement,round(Hit_rate_table_200$Top100,digits=3),cex=1.1,col="#FFDF34",font=2)


plot(1:13,MRR_table_200$Top10,type="b",ylim=c(0,0.04),axes=FALSE,col="#DB3124",ylab="")
complement<--c(0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.004,0.0034,0.004,0.004,0.004,0.004)
text(1:13,MRR_table_200$Top10+complement,round(MRR_table_200$Top10,digits=3),cex=1.1,col="#DB3124",font=2)
box()
axis(2,at=c(0,0.01,0.02,0.03,0.04),label=c("0",".01",".02",".03",".04"),lwd=1.5,las=1,cex.axis=1.4)
lines(1:13,type="b",MRR_table_200$Top50,col="#FC8C5A",lwd=1.5)
complement<-c(rep(0.002,13))
text(1:13,MRR_table_200$Top50+complement,round(MRR_table_200$Top50,digits=3),cex=1.1,col="#FC8C5A",font=2)
lines(1:13,type="b",MRR_table_200$Top100,col="#FFDF34",lwd=1.5)
complement<-c(rep(0.005,13))
text(1:13,MRR_table_200$Top100+complement-0.0003,round(MRR_table_200$Top100,digits=3),cex=1.1,col="#DB3124",font=2)
text(1:13,MRR_table_200$Top100+complement,round(MRR_table_200$Top100,digits=3),cex=1.1,col="#FFDF34",font=2)
legend("topright",legend = c("Top K=10","Top K=50","Top K=100"),lwd=1.5,pch=1,col=c("#DB3124","#FC8C5A","#FFDF34"))
title(ylab="MRR",cex.lab=1.6)

plot(1:13,MAP_table_200$Top10,type="b",ylim=c(0,0.04),lwd=1.5,axes=FALSE,col="#DB3124",ylab="")
complement<-c(-0.004,0.004,0.004,-0.004,-0.004,-0.004,-0.004,-0.004,-0.002,-0.003,-0.003,-0.003,-0.003)
text(1:13,MAP_table_200$Top10+complement,round(MAP_table_200$Top10,digits=3),cex=1.1,col="#DB3124",font=2)
title(ylab="MAP",cex.lab=1.6)
axis(2,at=c(0,0.01,0.02,0.03,0.04),label=c("0",".01",".02",".03",".04"),las=1,cex.axis=1.4)
box()
lines(1:13,type="b",MAP_table_200$Top50,col="#FC8C5A",lwd=1.5)
complement<-c(0.002,0.002,0.002,0.002,0.002,0.002,0.002,0,0.002,0.002,0.002,0.002,0.002)
text(1:13,MAP_table_200$Top50+complement,round(MAP_table_200$Top50,digits=3),cex=1.1,col="#FC8C5A",font=2)
lines(1:13,type="b",MAP_table_200$Top100,col="#FFDF34",lwd=1.5)
complement<-c(0.0055,-0.002,-0.002,0.0055,0.0055,0.0055,0.0055,0.0055,0.0055,0.0055,0.0055,0.0055,0.0055)
text(1:13,MAP_table_200$Top100+complement-0.0003,round(MAP_table_200$Top100,digits=3),cex=1.1,col="#DB3124",font=2)
text(1:13,MAP_table_200$Top100+complement,round(MAP_table_200$Top100,digits=3),cex=1.1,col="#FFDF34",font=2)

par(mar=c(6.5,5,1,1)+0.1)
plot(1:13,NDCG_table_200$Top10,type="b",ylim=c(0,0.065),lwd=1.5,axes=FALSE,col="#DB3124",ylab="",xlab="")
complement<--c(rep(0.005,13))
text(1:13,NDCG_table_200$Top10+complement,round(NDCG_table_200$Top10,digits=3),cex=1.1,col="#DB3124",font=2)
title(ylab="Mean NDCG",cex.lab=1.6)
axis(2,at=c(0,0.02,0.04,0.06),label=c("0",".02",".04",".06"),las=1,cex.axis=1.4)
box()
lines(1:13,type="b",NDCG_table_200$Top50,col="#FC8C5A")
complement<-c(-0.005,-0.005,-0.005,-0.005,-0.005,-0.005,-0.005,0.005,-0.005,-0.005,-0.005,-0.005,-0.005)
text(1:13,NDCG_table_200$Top50+complement,round(NDCG_table_200$Top50,digits=3),cex=1.1,col="#FC8C5A",font=2)
lines(1:13,type="b",NDCG_table_200$Top100,col="#FFDF34")
complement<-c(-0.005,-0.005,-0.005,-0.005,-0.005,0.004,0.004,0.009,0.004,0.004,0.004,0.004,0.004)
text(1:13,NDCG_table_200$Top100+complement-0.0005,round(NDCG_table_200$Top100,digits=3),cex=1.1,col="#DB3124",font=2)
text(1:13,NDCG_table_200$Top100+complement,round(NDCG_table_200$Top100,digits=3),cex=1.1,col="#FFDF34",font=2)
text(x=1:13,y=-0.028,labs,xpd=TRUE,cex=1.1,srt=90,font=2,col=col_MRR)



############################################################
#Figure 6 Upper
############################################################

#remember to reset order of table 200 if you run code above to draw figure 6 bottomright part.
#Hit_rate_table_200<-Hit_rate_table_200[c(1,12,3,10,7,9,8,2,4,13,6,11,5),]
#MRR_table_200<-MRR_table_200[c(1,12,3,10,7,9,8,2,4,13,6,11,5),]
#MAP_table_200<-MAP_table_200[c(1,12,3,10,7,9,8,2,4,13,6,11,5),]
#NDCG_table_200<-NDCG_table_200[c(1,12,3,10,7,9,8,2,4,13,6,11,5),]

rank_table_200<-data.frame(matrix(nrow=13,ncol=17))
colnames(rank_table_200)<-c("Method","Hit_rate_10","MRR_10","MAP_10","NDCG_10",
                            "Hit_rate_50","MRR_50","MAP_50","NDCG_50",
                            "Hit_rate_100","MRR_100","MAP_100","NDCG_100",
                            "Top1p","Top5p","Top10p","Mean_Rank")
rank_table_200$Method<-Method_names
rank_table_200$Hit_rate_10<-Hit_rate_table_200$Top10
rank_table_200$Hit_rate_50<-Hit_rate_table_200$Top50
rank_table_200$Hit_rate_100<-Hit_rate_table_200$Top100
rank_table_200$MRR_10<-MRR_table_200$Top10
rank_table_200$MRR_50<-MRR_table_200$Top50
rank_table_200$MRR_100<-MRR_table_200$Top100
rank_table_200$MAP_10<-MAP_table_200$Top10
rank_table_200$MAP_50<-MAP_table_200$Top50
rank_table_200$MAP_100<-MAP_table_200$Top100
rank_table_200$NDCG_10<-NDCG_table_200$Top10
rank_table_200$NDCG_50<-NDCG_table_200$Top50
rank_table_200$NDCG_100<-NDCG_table_200$Top100
Method_names
#fill the top percentages ranks
rank_table_200$`Top1p`<-c(38,17,117,8,25,1,8,78,22,8,26,9,15)
rank_table_200$`Top5p`<-c(67,82,234,31,108,19,19,158,53,37,63,51,46)
rank_table_200$`Top10p`<-c(99,128,313,50,175,33,22,212,86,55,116,71,70)

#create a rank table to calculate mean rank
rank_table_rank_200<-data.frame(matrix(nrow=nrow(rank_table_200),ncol=ncol(rank_table_200)))
colnames(rank_table_rank_200)<-colnames(rank_table_200)
rank_table_rank_200$Method<-rank_table_200$Method
for(i in 2:(ncol(rank_table_200)-1)){
  rank_table_rank_200[,i]<-rank(-rank_table_200[,i])
}
rank_table_rank_200$`Mean Rank`<-rowMeans(rank_table_rank_200[,2:(ncol(rank_table_rank_200)-1)])
min_rank<-min(rank_table_rank_200$`Mean Rank`)
max_rank<-max(rank_table_rank_200$`Mean Rank`)
rank_table_200$mean_rank_method<-paste0(rank_table_200$Method,": ",rank_table_200$formated_number)
rank_table_200$formated_number<-round(rank_table_rank_200$`Mean Rank`,digits = 2)
rank_table_200$Mean_Rank<-max(rank_table_rank_200$`Mean Rank`)+1-rank_table_rank_200$`Mean Rank`
#reset order
#rank_table_200=rank_table_200[c(1,3,8,13,5,9,11,4,7,12,2,10,6),]
rank_table_200$color_ramp<-colorRampPalette(c("green","red"))(13)

reactable(rank_table_200, defaultPageSize = 20,defaultColDef = colDef(
  minWidth = 35, headerStyle = list(fontWeight="normal",fontFamily="sans-serif",textAlign="right",fontSize='12px',borderBottom = "1px solid black",borderRight = "1px solid black")
),
columns=list(
  Method=colDef(name="Method",minWidth=35,align="center",style=function(value){
    if(value%in%c("RcisTarget","Pscan","HOMER")){
      return(list(fontWeight="normal",fontFamily="sans-serif",fontSize="12px",color="red"))
    }else{
      return(list(fontWeight="normal",fontFamily="sans-serif",fontSize="12px"))
    }
  },headerStyle = list(fontWeight="normal",fontFamily="sans-serif",textAlign="right",fontSize='12px',borderBottom = "1px solid black"),cell=data_bars(rank_table_200,text_position = "none")),
  Hit_rate_10=colDef(name="Hit Rate",minWidth=25,align="center",style=list(borderLeft="1px solid black"),headerStyle = list(fontWeight="normal",fontFamily="sans-serif",textAlign="right",fontSize='12px',borderBottom = "1px solid black",borderLeft = "1px solid black",borderRight = "1px solid black"),cell=data_bars(rank_table_200,fill_color = "#DB3124",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  Hit_rate_50=colDef(name="Hit Rate",minWidth=25,align="center",cell=data_bars(rank_table_200,fill_color = "#FC8C5A",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  Hit_rate_100=colDef(name="Hit Rate",minWidth=25,align="center",cell=data_bars(rank_table_200,fill_color = "#FFDF34",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  MRR_10=colDef(name="MRR",minWidth=25,align="center",cell=data_bars(rank_table_200,fill_color = "#DB3124",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  MRR_50=colDef(name="MRR",minWidth=25,align="center",cell=data_bars(rank_table_200,fill_color = "#FC8C5A",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  MRR_100=colDef(name="MRR",minWidth=25,align="center",cell=data_bars(rank_table_200,fill_color = "#FFDF34",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  MAP_10=colDef(name="MAP",minWidth=25,align="center",cell=data_bars(rank_table_200,fill_color = "#DB3124",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  MAP_50=colDef(name="MAP",minWidth=25,align="center",cell=data_bars(rank_table_200,fill_color = "#FC8C5A",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  MAP_100=colDef(name="MAP",minWidth=25,align="center",cell=data_bars(rank_table_200,fill_color = "#FFDF34",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  NDCG_10=colDef(name="Mean NDCG",minWidth=25,align="center",style=list(borderRight="1px solid black"),cell=data_bars(rank_table_200,fill_color = "#DB3124",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  NDCG_50=colDef(name="Mean NDCG",minWidth=25,align="center",style=list(borderRight="1px solid black"),cell=data_bars(rank_table_200,fill_color = "#FC8C5A",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  NDCG_100=colDef(name="Mean NDCG",minWidth=25,align="center",style=list(borderRight="1px solid black"),cell=data_bars(rank_table_200,fill_color = "#FFDF34",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  Top1p=colDef(name="Top 1%",minWidth=25,align="center",cell=data_bars(rank_table_200,fill_color = "#1c3e71",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  Top5p=colDef(name="Top 5%",minWidth=25,align="center",cell=data_bars(rank_table_200,fill_color = "#2c81be",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  Top10p=colDef(name="Top 10%",minWidth=25,align="center",style=list(borderRight="1px solid black"),cell=data_bars(rank_table_200,fill_color="#98C7DF",border_style="solid",border_width="thin",border_color = "black",text_position = "none")),
  mean_rank_method=colDef(name="Mean Rank",style=list(fontWeight="normal",fontFamily="sans-serif",fontSize="10px"),align="center",minWidth=50,headerStyle = list(fontWeight="normal",align="center",fontFamily="sans-serif",textAlign="right",borderBottom = "1px solid black"),cell=data_bars(rank_table_200,fill_by="Mean_Rank",text_color = "black",brighten_text = FALSE,fill_color_ref = "color_ramp",text_position = "inside-base",border_style="solid",border_width="thin",border_color = "black",)),
  Mean_Rank=colDef(show=FALSE),
  color_ramp=colDef(show=FALSE),
  formated_number=colDef(show=FALSE),
  mean_rank_method<-colDef(show=FALSE)
),columnGroups = list(
  colGroup(name="Top K = 10",columns=c("Hit_rate_10","MRR_10","MAP_10","NDCG_10"),headerStyle = list(textAlign="center",fontWeight="normal",fontFamily="sans-serif",textAlign="right",fontSize='12px',borderRight = "1px solid black",borderLeft="1px solid black")),
  colGroup(name="Top K = 50",columns=c("Hit_rate_50","MRR_50","MAP_50","NDCG_50"),headerStyle = list(textAlign="center",fontWeight="normal",fontFamily="sans-serif",textAlign="right",fontSize='12px',borderRight = "1px solid black")),
  colGroup(name="Top K = 100",columns=c("Hit_rate_100","MRR_100","MAP_100","NDCG_100"),headerStyle = list(textAlign="center",fontWeight="normal",fontFamily="sans-serif",textAlign="right",fontSize='12px',borderRight = "1px solid black")),
  colGroup(name="Top Percentiles",columns=c("Top1p","Top5p","Top10p"),headerStyle = list(textAlign="center",fontWeight="normal",fontFamily="sans-serif",textAlign="right",fontSize='12px',borderRight = "1px solid black"))
))

############################################################
#Figure 7
############################################################
work_dir<-"/Users/zeyulu/Dropbox/datasets/clean_results_new/"
knockTF_folder<-c("knockTF_200/","knockTF_600/","knockTF_1000/")
work_files<-list.files(paste0(work_dir,knockTF_folder[1]))
genes<-c(200,600,1000)

method_name<-c("BART","ChEA3","ChIP-Atlas","Cscan","Enrichr","HOMER","i-cisTarget",
               "Lisa","MAGIC","Pscan","RcisTarget","RegulatorTrail","TFEA.ChIP")

table_res<-read.csv(paste0(work_dir,knockTF_folder[1],work_files[1]),header=TRUE)
TF_index<-colnames(table_res)[2:ncol(table_res)]
TF_index2<-sapply(strsplit(TF_index,".",fixed=TRUE),function(x){return(x[[1]])})
TF_names<-sapply(strsplit(TF_index,"_",fixed=TRUE),function(x){return(x[[1]])})

TF_rank_table_200<-data.frame(matrix(nrow=13,ncol=length(TF_index2)))
colnames(TF_rank_table_200)<-TF_index2
rownames(TF_rank_table_200)<-method_name

TF_rank_table_600<-data.frame(matrix(nrow=13,ncol=length(TF_index2)))
colnames(TF_rank_table_600)<-TF_index2
rownames(TF_rank_table_600)<-method_name

TF_rank_table_1000<-data.frame(matrix(nrow=13,ncol=length(TF_index2)))
colnames(TF_rank_table_1000)<-TF_index2
rownames(TF_rank_table_1000)<-method_name

TF_rank_table_ratio_200<-data.frame(matrix(nrow=13,ncol=length(TF_index2)))
colnames(TF_rank_table_ratio_200)<-TF_index2
rownames(TF_rank_table_ratio_200)<-method_name

TF_rank_table_ratio_600<-data.frame(matrix(nrow=13,ncol=length(TF_index2)))
colnames(TF_rank_table_ratio_600)<-TF_index2
rownames(TF_rank_table_ratio_600)<-method_name

TF_rank_table_ratio_1000<-data.frame(matrix(nrow=13,ncol=length(TF_index2)))
colnames(TF_rank_table_ratio_1000)<-TF_index2
rownames(TF_rank_table_ratio_1000)<-method_name

for(k in 1:3){
  work_files<-list.files(paste0(work_dir,knockTF_folder[k]))

  method_name<-c("BART","ChEA3","ChIP-Atlas","Cscan","Enrichr","HOMER","i-cisTarget",
                 "Lisa","MAGIC","Pscan","RcisTarget","RegulatorTrail","TFEA.ChIP")

  for(i in 1:length(work_files)){
    print(i)
    table_res<-read.csv(paste0(work_dir,knockTF_folder[k],work_files[i]),header=TRUE)
    TF_index<-colnames(table_res)[1:ncol(table_res)]
    TF_index2<-sapply(strsplit(TF_index,".",fixed=TRUE),function(x){return(x[[1]])})
    TF_names<-sapply(strsplit(TF_index,"_",fixed=TRUE),function(x){return(x[[1]])})
    rank_num<-c()
    TF_in_list<-c()
    for(j in 1:length(TF_index)){
      TF_rank_list<-table_res[,TF_index[j]]
      rank_num<-c(rank_num,sum(!is.na(TF_rank_list)))
      if(k==1){
        TF_rank_table_200[i,TF_index2[j]]<-which(TF_rank_list==TF_names[j])[1]
        TF_rank_table_ratio_200[i,TF_index2[j]]<-(which(TF_rank_list==TF_names[j])[1])/sum(!is.na(TF_rank_list))
      }else if(k==2){
        TF_rank_table_600[i,TF_index2[j]]<-which(TF_rank_list==TF_names[j])[1]
        TF_rank_table_ratio_600[i,TF_index2[j]]<-(which(TF_rank_list==TF_names[j])[1])/sum(!is.na(TF_rank_list))
      }else if(k==3){
        TF_rank_table_1000[i,TF_index2[j]]<-which(TF_rank_list==TF_names[j])[1]
        TF_rank_table_ratio_1000[i,TF_index2[j]]<-(which(TF_rank_list==TF_names[j])[1])/sum(!is.na(TF_rank_list))
      }
    }
  }
}


max_rank_table_200<-list()
for(i in 1:length(method_name)){
  for(j in 1:length(method_name)){
    max_rank_table_200[[paste0(i,"_",j)]]=c(0)
  }
}

MRR_min_rank_table_200_10<-data.frame(matrix(nrow=13,ncol=13))
colnames(MRR_min_rank_table_200_10)<-method_name
colnames(MRR_min_rank_table_200_10)<-method_name

MRR_min_rank_table_200_50<-data.frame(matrix(nrow=13,ncol=13))
colnames(MRR_min_rank_table_200_50)<-method_name
colnames(MRR_min_rank_table_200_50)<-method_name

MRR_min_rank_table_200_100<-data.frame(matrix(nrow=13,ncol=13))
colnames(MRR_min_rank_table_200_100)<-method_name
colnames(MRR_min_rank_table_200_100)<-method_name

Jaccard_rank_table_200_10<-data.frame(matrix(nrow=13,ncol=13))
colnames(Jaccard_rank_table_200_10)<-method_name
colnames(Jaccard_rank_table_200_10)<-method_name

Jaccard_rank_table_200_50<-data.frame(matrix(nrow=13,ncol=13))
colnames(Jaccard_rank_table_200_50)<-method_name
colnames(Jaccard_rank_table_200_50)<-method_name

Jaccard_rank_table_200_100<-data.frame(matrix(nrow=13,ncol=13))
colnames(Jaccard_rank_table_200_100)<-method_name
colnames(Jaccard_rank_table_200_100)<-method_name

combo_table_200_10<-data.frame(matrix(nrow=13,ncol=13))
colnames(combo_table_200_10)<-method_name
rownames(combo_table_200_10)<-method_name

combo_table_200_1p<-data.frame(matrix(nrow=13,ncol=13))
colnames(combo_table_200_1p)<-method_name
rownames(combo_table_200_1p)<-method_name

combo_table_200_50<-data.frame(matrix(nrow=13,ncol=13))
colnames(combo_table_200_50)<-method_name
rownames(combo_table_200_50)<-method_name

combo_table_200_5p<-data.frame(matrix(nrow=13,ncol=13))
colnames(combo_table_200_5p)<-method_name
rownames(combo_table_200_5p)<-method_name

combo_table_200_100<-data.frame(matrix(nrow=13,ncol=13))
colnames(combo_table_200_100)<-method_name
rownames(combo_table_200_100)<-method_name

combo_table_200_10p<-data.frame(matrix(nrow=13,ncol=13))
colnames(combo_table_200_10p)<-method_name
rownames(combo_table_200_10p)<-method_name

par(mfrow=c(4,3))

Jac_cal<-function(li1,li2){
  li1<-li1[!is.na(li1)]
  li2<-li2[!is.na(li2)]
  Jac_val<-length(intersect(li1,li2))/length(union(li1,li2))
  if(is.na(Jac_val)){
    return(0)
  }else{
    return(Jac_val)
  }
}



for(k in 1:3){
  work_files<-list.files(paste0(work_dir,knockTF_folder[k]))

  method_name<-c("BART","ChEA3","ChIP-Atlas","Cscan","Enrichr","HOMER","i-cisTarget",
                 "Lisa","MAGIC","Pscan","RcisTarget","RegulatorTrail","TFEA.ChIP")

  for(i in 1:length(work_files)){
    for(m in 1:length(work_files)){
      print(i)
      table_res1<-read.csv(paste0(work_dir,knockTF_folder[k],work_files[i]),header=TRUE)
      table_res2<-read.csv(paste0(work_dir,knockTF_folder[k],work_files[m]),header=TRUE)
      TF_index<-colnames(table_res1)[1:ncol(table_res1)]
      TF_index3<-colnames(table_res2)[1:ncol(table_res2)]
      TF_index<-intersect(TF_index,TF_index3)
      TF_index2<-sapply(strsplit(TF_index,".",fixed=TRUE),function(x){return(x[[1]])})
      TF_names<-sapply(strsplit(TF_index,"_",fixed=TRUE),function(x){return(x[[1]])})
      rank_num<-c()
      TF_in_list<-c()

      for(j in 1:length(TF_index)){
        TF_rank_list1<-table_res1[,TF_index[j]]
        TF_rank_list2<-table_res2[,TF_index[j]]
        if(k==1){
          Jaccard_rank_table_200_10[i,m]<-Jac_cal(TF_rank_list1[1:10],TF_rank_list2[1:10])
          Jaccard_rank_table_200_50[i,m]<-Jac_cal(TF_rank_list1[1:50],TF_rank_list2[1:50])
          Jaccard_rank_table_200_100[i,m]<-Jac_cal(TF_rank_list1[1:100],TF_rank_list2[1:100])
        }else if(k==2){
        }else if(k==3){
        }
      }
    }
  }
}




for(i in 1:nrow(TF_rank_table_200)){
  for(k in 1:nrow(TF_rank_table_200)){
    rank_list=c()
    ratio_list=c()
    for(j in 1:ncol(TF_rank_table_200)){
      rank_list=c(rank_list,min(TF_rank_table_200[i,j],TF_rank_table_200[k,j],na.rm=TRUE))
      ratio_list=c(ratio_list,min(TF_rank_table_ratio_200[i,j],TF_rank_table_ratio_200[k,j],na.rm=TRUE))
    }
    rank_10_list<-rank_list
    rank_10_list[rank_10_list>10]=NA
    rank_50_list<-rank_list
    rank_50_list[rank_50_list>50]=NA
    rank_100_list<-rank_list
    rank_100_list[rank_100_list>100]=NA


    MRR_min_rank_table_200_10[i,k]<-sum(1/rank_10_list,na.rm=TRUE)/570
    MRR_min_rank_table_200_50[i,k]<-sum(1/rank_50_list,na.rm=TRUE)/570
    MRR_min_rank_table_200_100[i,k]<-sum(1/rank_100_list,na.rm=TRUE)/570

    combo_table_200_10[i,k]<-sum(rank_list<=10,na.rm=TRUE)
    combo_table_200_1p[i,k]<-sum(ratio_list<=0.01,na.rm=TRUE)

    combo_table_200_50[i,k]<-sum(rank_list<=50,na.rm=TRUE)
    combo_table_200_5p[i,k]<-sum(ratio_list<=0.05,na.rm=TRUE)

    combo_table_200_100[i,k]<-sum(rank_list<=100,na.rm=TRUE)
    combo_table_200_10p[i,k]<-sum(ratio_list<=0.1,na.rm=TRUE)
  }
}

Jac_cal<-function(li1,li2){
  return(length(intersect(li1,li2))/length(union(li1,li2)))
}

MRR_min_rank_table_200_10<-round(MRR_min_rank_table_200_10,3)
MRR_min_rank_table_200_50<-round(MRR_min_rank_table_200_50,3)
MRR_min_rank_table_200_100<-round(MRR_min_rank_table_200_100,3)

Hit_rate_table_200_10<-round(combo_table_200_10/570,3)
Hit_rate_table_200_50<-round(combo_table_200_50/570,3)
Hit_rate_table_200_100<-round(combo_table_200_100/570,3)

Jaccard_rank_table_200_10<-round(Jaccard_rank_table_200_10,3)
Jaccard_rank_table_200_50<-round(Jaccard_rank_table_200_50,3)
Jaccard_rank_table_200_100<-round(Jaccard_rank_table_200_100,3)

create_mask<-function(mat,adj=0){
  n=ncol(mat)
  mask_upper <- matrix(NA, ncol = n, nrow = n)
  for (i in 1:(n)) {
    for (j in 1:(n)) {
      if(adj==0){
        if (j >= i) {
          mask_upper[i, j] <- TRUE
        } else {
          mask_upper[i, j] <- NA
        }
      }else{
        if (j > i) {
          mask_upper[i, j] <- TRUE
        } else {
          mask_upper[i, j] <- NA
        }
      }
    }
  }
  return(mask_upper)
}

plot_heatmap<-function(table,adj=0){
  table_mat<-create_mask(as.matrix(table),adj)
  table_upper<-as.matrix(table)*table_mat
  return(table_upper)
}

draw_rect<-function(n,table,cex_text=0.8){
  for(i in 1:n) {
    for(j in 1:n) {
      if (j >= i) {
        rect(i-0.5, j-0.5, i+0.5, j+0.5)
        text(i,j,table[i,j],col="black",cex=cex_text,font=2)
      }
    }
  }
}
###########
collab<-c("black")

layout_matrix=matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),5,3,byrow=TRUE)
layout(layout_matrix,heights=c(1.01,0.9,0.9,0.9,0.9),widths=c(1,0.92,0.92))
par(mar = c(1,3,4,1) + 0.1)
h<-image(1:13,1:13,plot_heatmap(combo_table_200_10),ylab="",axes=FALSE,col=colorRampPalette(c("yellow","red"))(20))
draw_rect(13,combo_table_200_10)
axis(3, at = seq(1, 13, by = 1), labels = method_name, las=2,cex.axis=0.7,tcl=-0.2,font.axis=2,col.axis=collab)
axis(2, at = seq(1, 13, by = 1), labels = method_name, las = 1,cex.axis=0.7,tcl=-0.2,font.axis=2,col.axis=collab)
mtext("Top 10", side=1, line=0,cex=0.7,font=2)
segments(1,8,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)

par(mar=c(1,1,4,1))
h<-image(1:13,1:13,plot_heatmap(combo_table_200_50),axes=FALSE,col=colorRampPalette(c("yellow","red"))(20))
draw_rect(13,combo_table_200_50)
axis(3, at = seq(1, 13, by = 1), labels = method_name, las=2,cex.axis=0.7,tcl=-0.2,font.axis=2,col.axis=collab)
mtext("Top 50", side=1, line=0,cex=0.7,font=2)
segments(5,13,8.5,5.3,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)

par(mar=c(1,1,4,1))
h<-image(1:13,1:13,plot_heatmap(combo_table_200_100),axes=FALSE,col=colorRampPalette(c("yellow","red"))(20))
draw_rect(13,combo_table_200_100)
axis(3, at = seq(1, 13, by = 1), labels = method_name, las=2,cex.axis=0.7,tcl=-0.2,font.axis=2,col.axis=collab)
mtext("Top 100", side=1, line=0,cex=0.7,font=2)
legend("bottomright",legend=c("Low","High"),fill=c("yellow","red"),title="Top N")
segments(5,13,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)

par(mar=c(1,3,1,1))
h<-image(1:13,1:13,plot_heatmap(combo_table_200_1p),ylab="",axes=FALSE,col=colorRampPalette(c("lightblue","lightblue4"))(20))
draw_rect(13,combo_table_200_1p)
axis(2, at = seq(1, 13, by = 1), labels = method_name, las = 1,cex.axis=0.7,tcl=-0.2,font.axis=2,col.axis=collab)
mtext("Top 1%", side=1, line=0,cex=0.7,font=2)
segments(3,8,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)

par(mar=c(1,1,1,1))
h<-image(1:13,1:13,plot_heatmap(combo_table_200_5p),axes=FALSE,col=colorRampPalette(c("lightblue","lightblue4"))(20))
draw_rect(13,combo_table_200_5p)
mtext("Top 5%", side=1, line=0,cex=0.7,font=2)
segments(3,5,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)

par(mar=c(1,1,1,1))
h<-image(1:13,1:13,plot_heatmap(combo_table_200_10p),axes=FALSE,col=colorRampPalette(c("lightblue","lightblue4"))(20))
draw_rect(13,combo_table_200_10p)
mtext("Top 10%", side=1, line=0,cex=0.7,font=2)
segments(3,5,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)
legend("bottomright",legend=c("Low","High"),fill=c("lightblue","lightblue4"),title="Top %")
mtext("200 Input Genes",side=1,outer=TRUE,line=1,font=2,cex=0.8,col.axis=collab)


par(mar=c(1,3,1,1))
h<-image(1:13,1:13,plot_heatmap(Hit_rate_table_200_10),ylab="",axes=FALSE,col=colorRampPalette(c("#FAE2C5","#EB8E47"))(20))
draw_rect(13,Hit_rate_table_200_10,0.7)
axis(2, at = seq(1, 13, by = 1), labels = method_name, las = 1,cex.axis=0.7,tcl=-0.2,font.axis=2,col.axis=collab)
mtext("K = 10", side=1, line=0,cex=0.7,font=2)
segments(1,8,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)

par(mar=c(1,1,1,1))
h<-image(1:13,1:13,plot_heatmap(Hit_rate_table_200_50),axes=FALSE,col=colorRampPalette(c("#FAE2C5","#EB8E47"))(20))
draw_rect(13,Hit_rate_table_200_50,0.7)
mtext("K = 50", side=1, line=0,cex=0.7,font=2)
segments(5,13,8.5,5.3,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)

par(mar=c(1,1,1,1))
h<-image(1:13,1:13,plot_heatmap(Hit_rate_table_200_100),axes=FALSE,col=colorRampPalette(c("#FAE2C5","#EB8E47"))(20))
draw_rect(13,Hit_rate_table_200_100,0.7)
mtext("K = 100", side=1, line=0,cex=0.7,font=2)
segments(5,13,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)
legend("bottomright",legend=c("Low","High"),fill=c("#FAE2C5","#EB8E47"),title="Hit Rate @ K")

par(mar=c(1,3,1,1))
h<-image(1:13,1:13,plot_heatmap(MRR_min_rank_table_200_10),ylab="",axes=FALSE,col=colorRampPalette(c("#EFE2ED","#A977A6"))(20))
draw_rect(13,MRR_min_rank_table_200_10,0.7)
axis(2, at = seq(1, 13, by = 1), labels = method_name, las = 1,cex.axis=0.7,tcl=-0.2,font.axis=2,col.axis=collab)
mtext("K = 10", side=1, line=0,cex=0.7,font=2)
segments(1,3,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)

par(mar=c(1,1,1,1))
h<-image(1:13,1:13,plot_heatmap(MRR_min_rank_table_200_50),axes=FALSE,col=colorRampPalette(c("#EFE2ED","#A977A6"))(20))
draw_rect(13,MRR_min_rank_table_200_50,0.7)
mtext("K = 50", side=1, line=0,cex=0.7,font=2)
segments(1,3,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)

par(mar=c(1,1,1,1))
h<-image(1:13,1:13,plot_heatmap(MRR_min_rank_table_200_100),axes=FALSE,col=colorRampPalette(c("#EFE2ED","#A977A6"))(20))
draw_rect(13,MRR_min_rank_table_200_100,0.7)
mtext("K = 100", side=1, line=0,cex=0.7,font=2)
segments(1,3,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)
legend("bottomright",legend=c("Low","High"),fill=c("#EFE2ED","#A977A6"),title="MRR @ K")

par(mar=c(1,3,1,1))
h<-image(1:13,1:13,plot_heatmap(Jaccard_rank_table_200_10,-1),ylab="",axes=FALSE,col=colorRampPalette(c("#D3EEE2","#2AA371"))(20))
draw_rect(13,Jaccard_rank_table_200_10,0.7)
axis(2, at = seq(1, 13, by = 1), labels = method_name, las = 1,cex.axis=0.7,tcl=-0.2,font.axis=2,col.axis=collab)
mtext("K = 10", side=1, line=0,cex=0.7,font=2)
segments(7,9,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)

par(mar=c(1,1,1,1))
h<-image(1:13,1:13,plot_heatmap(Jaccard_rank_table_200_50,-1),axes=FALSE,col=colorRampPalette(c("#D3EEE2","#2AA371"))(20))
draw_rect(13,Jaccard_rank_table_200_50,0.7)
mtext("K = 50", side=1, line=0,cex=0.7,font=2)
segments(2,9,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)

par(mar=c(1,1,1,1))
h<-image(1:13,1:13,plot_heatmap(Jaccard_rank_table_200_100,-1),axes=FALSE,col=colorRampPalette(c("#D3EEE2","#2AA371"))(20))
draw_rect(13,Jaccard_rank_table_200_100,0.7)
mtext("K = 100", side=1, line=0,cex=0.7,font=2)
segments(2,9,8,5,lwd=0.7)
text(8.5,5,"Max",font=2,cex=0.8)
legend("bottomright",legend=c("Low","High"),fill=c("#D3EEE2","#2AA371"),title="Jaccard @ K")
title(main="Performance Evaluation by Combining Two Methods",outer=TRUE,line=1,font=2,font.lab=2,col.axis=collab)
mtext("200 Input Genes",side=1,outer=TRUE,line=1,font=2,cex=0.8,col.axis=collab)






############################################################
#Figure 8
############################################################
work_dir<-"/Users/zeyulu/Dropbox/datasets/clean_results_new/"
knockTF_folder<-c("knockTF_200/","knockTF_600/","knockTF_1000/")
genes<-c(200,600,1000)

par(mfrow=c(2,2),oma = c(2,2,2,0) + 0.1,mar = c(2,2,1,1) + 0.1)
method_name<-c("BART","ChEA3","ChIP-Atlas","Cscan","Enrichr","HOMER","i-cisTarget",
               "Lisa","MAGIC","Pscan","RcisTarget","RegulatorTrail","TFEA.ChIP")
mappable_table<-data.frame(matrix(nrow=13,ncol=3))
colnames(mappable_table)<-c("top200","top600","top1000")
rownames(mappable_table)<-method_name

na_table<-data.frame(matrix(nrow=13,ncol=3))
colnames(na_table)<-c("top200","top600","top1000")
rownames(na_table)<-method_name

ratio_table<-data.frame(matrix(nrow=13,ncol=3))
colnames(ratio_table)<-c("200","600","1000")
rownames(ratio_table)<-method_name

abs_table<-data.frame(matrix(nrow=13,ncol=3))
colnames(abs_table)<-c("200","600","1000")
rownames(abs_table)<-method_name

TF_table<-read.csv("/Users/zeyulu/Dropbox/datasets/clean_results_new copy/TF.csv")
TF_names_all<-TF_table$HGNC.symbol

top10_rank_list_200<-list()
top50_rank_list_600<-list()
top100_rank_list_1000<-list()

for(k in 1:3){
  work_files<-list.files(paste0(work_dir,knockTF_folder[k]))

  method_name<-c("BART","ChEA3","ChIP-Atlas","Cscan","Enrichr","HOMER","i-cisTarget",
                 "Lisa","MAGIC","Pscan","RcisTarget","RegulatorTrail","TFEA.ChIP")

  top10_rank_list<-list()
  top50_rank_list<-list()
  top100_rank_list<-list()
  for(i in 1:length(work_files)){
    print(i)
    table_res<-read.csv(paste0(work_dir,knockTF_folder[k],work_files[i]),header=TRUE)
    TF_index<-colnames(table_res)[2:ncol(table_res)]
    TF_names<-sapply(strsplit(TF_index,"_",fixed=TRUE),function(x){return(x[[1]])})
    rank_list<-c()
    rank_num<-c()
    TF_in_list<-c()
    top_10_list<-c()
    top_50_list<-c()
    top_100_list<-c()
    for(j in 1:length(TF_index)){
      TF_rank_list<-table_res[,TF_index[j]]
      top_10_list<-c(top_10_list,unique(TF_rank_list[1:10]))
      top_50_list<-c(top_50_list,unique(TF_rank_list[1:50]))
      top_100_list<-c(top_100_list,unique(TF_rank_list[1:100]))

      TF_in_list<-c(TF_in_list,length(which(unique(TF_rank_list)%in%TF_names_all)))
      rank_num<-c(rank_num,sum(!is.na(TF_rank_list)))
      rank_list<-c(rank_list,which(TF_rank_list==TF_names[j])[1])
      ratio_list<-rank_list/rank_num
    }

    top10_rank_list[[method_name[i]]]<-table(top_10_list)
    top50_rank_list[[method_name[i]]]<-table(top_50_list)
    top100_rank_list[[method_name[i]]]<-table(top_100_list)

    mappable_table[method_name[i],k]<-max(TF_in_list,rm.na=TRUE)
    top1p<-sum(ratio_list<=0.01,na.rm=TRUE)
    top5p<-sum(ratio_list<=0.05,na.rm=TRUE)
    top10p<-sum(ratio_list<=0.1,na.rm=TRUE)



    top10<-sum(rank_list<=10,na.rm=TRUE)
    top50<-sum(rank_list<=50,na.rm=TRUE)
    top100<-sum(rank_list<=100,na.rm=TRUE)

    na_val<-sum(is.na(rank_list))

    ratio_table[method_name[i],k]<-top1p


    abs_table[method_name[i],k]<-top10

    na_table[method_name[i],k]<-na_val
  }
  top10_rank_list_200[[k]]<-top10_rank_list
  top50_rank_list_600[[k]]<-top50_rank_list
  top100_rank_list_1000[[k]]<-top100_rank_list
}

method_rowm<-c()
method_nam<-c()

for(i in 1:13){
  method_nam<-c(method_nam,min(na_table[i,1:3]))
}
work_files<-list.files("/Users/zeyulu/Dropbox/datasets/clean_results_new/knockTF_200/")
for(i in 1:13){
  input_df<-as.data.frame(fread(paste0("/Users/zeyulu/Dropbox/datasets/clean_results_new/knockTF_200/",work_files[i])))
  method_rowm<-c(method_rowm,length(unique(input_df[,10])))
}


method_rowm
names(method_rowm)<-rownames(mappable_table)
years<-c(2018,2019,2018,2013,2013,2010,2015,2020,2020,2009,2016,2017,2019)

method_name<-c("BART","ChEA3","ChIP-Atlas","Cscan","Enrichr","HOMER","i-cisTarget",
               "Lisa","MAGIC","Pscan","RcisTarget","RegulatorTrail","TFEA.ChIP")

coverage_table<-data.frame(years=years,nas=method_nam,avail=method_rowm,top10=abs_table$'200',top1p=ratio_table$'200')
rownames(coverage_table)<-method_name
coverage_table

library(basicPlotteR)

abs_table
ratio_table
coverage_table

par(mfrow=c(2,2),oma = c(2,2,2,0) + 0.1,mar = c(0.5,1,1,1) + 0.1,mgp=c(2,0.5,0))
plot(top10~years,data=coverage_table,main="Top 10 TRs (200 input genes)",
     font.lab=2,cex=0.8,pch=21,xaxt="n",col="black",bg=c("blue","blue","blue","blue","blue","red",
                                                         "blue","blue","blue","red","red","blue","blue"))
addTextLabels(coverage_table$years,coverage_table$top10,method_name,cex.label=1,col.label="black",lwd=0.75)
arrows(x0=2014,y0=20,x1=2017,y1=24,col="black",length=0.05,lwd=1.5)
text(2015.2,23.1,"Increase",srt=25,cex=0.8)
text(2009,36.7,"A",font=2)

plot(top1p~years,data=coverage_table,main="Top 1% TRs (200 input genes)",
     font.lab=2,cex=0.8,pch=21,xaxt="n",col="black",bg=c("blue","blue","blue","blue","blue","red",
                                                         "blue","blue","blue","red","red","blue","blue"))
addTextLabels(coverage_table$years,coverage_table$top1p,method_name,cex.label=1,col.label="black",lwd=0.75)
arrows(x0=2012,y0=35,x1=2015,y1=55,col="black",length=0.05,lwd=1.5)
text(2013,48,"Increase",srt=36,cex=0.8)
text(2009,114,"B",font=2)
plot(nas~years,data=coverage_table,main="Perturbed TRs Not Found",
     font.lab=2,cex=0.8,pch=21,col="black",bg=c("blue","blue","blue","blue","blue","red",
                                                "blue","blue","blue","red","red","blue","blue"))
addTextLabels(coverage_table$years,coverage_table$nas,method_name,cex.label=1,col.label="black",lwd=0.75)
arrows(x0=2013,y0=370,x1=2016,y1=250,col="black",length=0.05,lwd=1.5)
text(2014.8,333,"Decrease",srt=-44,cex=0.8)
text(2009,545,"C",font=2)
legend("topright",pch=c(21,21),col="black",pt.bg=c("blue","red"),legend=c("NGS","Motif"),
       pt.cex=1,cex=1)
plot(avail~years,data=coverage_table,main="Number of Unique TRs",
     font.lab=2,cex=0.8,pch=21,col="black",bg=c("blue","blue","blue","blue","blue","red",
                                                "blue","blue","blue","red","red","blue","blue"))
addTextLabels(coverage_table$years,coverage_table$avail,method_name,cex.label=1,col.label="black",lwd=0.75)
arrows(x0=2012,y0=400,x1=2015,y1=700,col="black",length=0.05,lwd=1.5)
text(2013.2,630,"Increase",srt=34,cex=0.8)
text(2009,1650,"D",font=2)

coverage_table

title(xlab="Publication date",outer=TRUE,font.lab=2,line=1,cex.main=1.3,cex.lab=1.4)
title(ylab="Number of TRs",outer=TRUE, font.lab=2,line=0.6,cex.main=1.4,cex.lab=1.4)


############################################################
#Figure 9
############################################################
library(ggplot2)

MRR_diff_table_top10<-data.frame(matrix(nrow=13,ncol=3))
MRR_diff_table_top50<-data.frame(matrix(nrow=13,ncol=3))
MRR_diff_table_top100<-data.frame(matrix(nrow=13,ncol=3))
colnames(MRR_diff_table_top10)<-c("g200","g600","g1000")
colnames(MRR_diff_table_top50)<-c("g200","g600","g1000")
colnames(MRR_diff_table_top100)<-c("g200","g600","g1000")

MRR_diff_table_top10$g200<-0
MRR_diff_table_top10$g600<-MRR_table_600$Top10-MRR_table_200$Top10
MRR_diff_table_top10$g1000<-MRR_table_1000$Top10-MRR_table_200$Top10

MRR_diff_table_top50$g200<-0
MRR_diff_table_top50$g600<-MRR_table_600$Top50-MRR_table_200$Top50
MRR_diff_table_top50$g1000<-MRR_table_1000$Top50-MRR_table_200$Top50

MRR_diff_table_top100$g200<-0
MRR_diff_table_top100$g600<-MRR_table_600$Top100-MRR_table_200$Top100
MRR_diff_table_top100$g1000<-MRR_table_1000$Top100-MRR_table_200$Top100

MAP_diff_table_top10<-data.frame(matrix(nrow=13,ncol=3))
MAP_diff_table_top50<-data.frame(matrix(nrow=13,ncol=3))
MAP_diff_table_top100<-data.frame(matrix(nrow=13,ncol=3))
colnames(MAP_diff_table_top10)<-c("g200","g600","g1000")
colnames(MAP_diff_table_top50)<-c("g200","g600","g1000")
colnames(MAP_diff_table_top100)<-c("g200","g600","g1000")

MAP_diff_table_top10$g200<-0
MAP_diff_table_top10$g600<-MAP_table_600$Top10-MAP_table_200$Top10
MAP_diff_table_top10$g1000<-MAP_table_1000$Top10-MAP_table_200$Top10

MAP_diff_table_top50$g200<-0
MAP_diff_table_top50$g600<-MAP_table_600$Top50-MAP_table_200$Top50
MAP_diff_table_top50$g1000<-MAP_table_1000$Top50-MAP_table_200$Top50

MAP_diff_table_top100$g200<-0
MAP_diff_table_top100$g600<-MAP_table_600$Top100-MAP_table_200$Top100
MAP_diff_table_top100$g1000<-MAP_table_1000$Top100-MAP_table_200$Top100

Hit_rate_diff_table_top10<-data.frame(matrix(nrow=13,ncol=3))
Hit_rate_diff_table_top50<-data.frame(matrix(nrow=13,ncol=3))
Hit_rate_diff_table_top100<-data.frame(matrix(nrow=13,ncol=3))
colnames(Hit_rate_diff_table_top10)<-c("g200","g600","g1000")
colnames(Hit_rate_diff_table_top50)<-c("g200","g600","g1000")
colnames(Hit_rate_diff_table_top100)<-c("g200","g600","g1000")

Hit_rate_diff_table_top10$g200<-0
Hit_rate_diff_table_top10$g600<-Hit_rate_table_600$Top10-Hit_rate_table_200$Top10
Hit_rate_diff_table_top10$g1000<-Hit_rate_table_1000$Top10-Hit_rate_table_200$Top10

Hit_rate_diff_table_top50$g200<-0
Hit_rate_diff_table_top50$g600<-Hit_rate_table_600$Top50-Hit_rate_table_200$Top50
Hit_rate_diff_table_top50$g1000<-Hit_rate_table_1000$Top50-Hit_rate_table_200$Top50

Hit_rate_diff_table_top100$g200<-0
Hit_rate_diff_table_top100$g600<-Hit_rate_table_600$Top100-Hit_rate_table_200$Top100
Hit_rate_diff_table_top100$g1000<-Hit_rate_table_1000$Top100-Hit_rate_table_200$Top100

NDCG_diff_table_top10<-data.frame(matrix(nrow=13,ncol=3))
NDCG_diff_table_top50<-data.frame(matrix(nrow=13,ncol=3))
NDCG_diff_table_top100<-data.frame(matrix(nrow=13,ncol=3))
colnames(NDCG_diff_table_top10)<-c("g200","g600","g1000")
colnames(NDCG_diff_table_top50)<-c("g200","g600","g1000")
colnames(NDCG_diff_table_top100)<-c("g200","g600","g1000")

NDCG_diff_table_top10$g200<-0
NDCG_diff_table_top10$g600<-NDCG_table_600$Top10-NDCG_table_200$Top10
NDCG_diff_table_top10$g1000<-NDCG_table_1000$Top10-NDCG_table_200$Top10

NDCG_diff_table_top50$g200<-0
NDCG_diff_table_top50$g600<-NDCG_table_600$Top50-NDCG_table_200$Top50
NDCG_diff_table_top50$g1000<-NDCG_table_1000$Top50-NDCG_table_200$Top50

NDCG_diff_table_top100$g200<-0
NDCG_diff_table_top100$g600<-NDCG_table_600$Top100-NDCG_table_200$Top100
NDCG_diff_table_top100$g1000<-NDCG_table_1000$Top100-NDCG_table_200$Top100

Hit_rate_table_200
Hit_rate_table_600

par(mfrow=c(4,3),oma = c(3,3,2,0) + 0.1,mar = c(5,4,4,4) + 0.1,xpd=TRUE,bty="o")

data <- mutate(MRR_diff_table_top10, RowID = row_number())
data <- mutate(data, RowID = row_number())

# Reshape the data to a long format
data_long <- tidyr::pivot_longer(data, -RowID, names_to = "Column", values_to = "Value")

# Plotting
data_long$Column <- factor(data_long$Column, levels = names(data[-ncol(data)]))

# Plotting
ggplot(data_long, aes(x = Column, y = Value, group = RowID, color = as.factor(RowID))) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Line Plot of Columns by Rows", x = "Column", y = "Value", color = "Row ID") +
  scale_color_discrete(name = "Row ID")


tables<-list("Hit Rate K=10"=Hit_rate_diff_table_top10,"Hit Rate K=50"=Hit_rate_diff_table_top50,"Hit Rate K=100"=Hit_rate_diff_table_top100,
             "MRR K=10"=MRR_diff_table_top10,"MRR K=50"=MRR_diff_table_top50,"MRR K=100"=MRR_diff_table_top100,
             "MAP K=10"=MAP_diff_table_top10,"MAP K=50"=MAP_diff_table_top50,"MAP K=100"=MAP_diff_table_top100,
             "NDCG K=10"=NDCG_diff_table_top10,"NDCG K=50"=NDCG_diff_table_top50,"NDCG K=100"=NDCG_diff_table_top100)

tables <- lapply(tables, function(x) mutate_all(x, as.numeric))

combined_data <- bind_rows(lapply(names(tables), function(t_id) {
  mutate(tables[[t_id]], rowName = Hit_rate_table_200$Method, TableID = t_id)
}), .id = "id")  # '.id' is used to create a new column from the names of the list elements

# Reshape the combined data to a long format
combined_long <- tidyr::pivot_longer(combined_data, c(g200,g600,g1000), names_to = "Column", values_to = "Value")

# Convert the 'Column' to a factor with the levels in the desired order (assuming all tables have the same column names)
combined_long$Column <- factor(combined_long$Column, levels = names(tables[[1]]))
combined_long$TableID <- factor(combined_long$TableID, levels = c("Hit Rate K=10", "Hit Rate K=50", "Hit Rate K=100", "MAP K=10", "MAP K=50", "MAP K=100", "MRR K=10", "MRR K=50", "MRR K=100", "NDCG K=10", "NDCG K=50", "NDCG K=100"))
combined_long$rowName <- factor(combined_long$rowName, levels = c(Hit_rate_table_200$Method))

Hit_rate_diff_table_top10

strip <- strip_themed(background_x = elem_list_rect(fill = c("#DB3124","#FC8C5A","#FFDF34","#DB3124","#FC8C5A","#FFDF34","#DB3124","#FC8C5A","#FFDF34","#DB3124","#FC8C5A","#FFDF34")))

color_more<-c("Lisa"="#FF0000",
              "Cscan"="#00FF00",
              "HOMER"="#0000FF",
              "i-cisTarget"="#008000",
              "ChIP-Atlas"="#FFA500",
              "TFEA.ChIP"="#800080",
              "Pscan"="#00FFFF",
              "Enrichr"="#FF00FF",
              "RcisTarget"="#FFC0C8",
              "ChEA3"="#008080",
              "MAGIC"="#808000",
              "BART"="#800000")


library(ggh4x)

combined_long[combined_long$rowName=="HOMER",]

Hit_rate_table_600
Hit_rate_table_200



gg <- ggplot(combined_long, aes(x = Column, y = Value, group = rowName, color = as.factor(rowName))) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(title = "", x = "Number of input genes", y = "Metric score deviation", color = "Methods") +
  scale_color_manual(name = "Methods", values=color_more) +
  scale_x_discrete(labels = c("g200" = "200", "g600" = "600", "g1000" = "1000")) +
  facet_wrap2(~ TableID, scales = "free_y", ncol = 3, strip=strip)
print(gg)
#############################
#Stouffer test
#############################
work_dir
knockTF_folder<-c("knockTF_200/","knockTF_600/","knockTF_1000/")
genes<-c(200,600,1000)

method_name<-c("BART","CHEA3","CHIPATLAS","CSCAN","ENRICHR","HOMER","ICIS",
               "LISA","MAGIC","PSCAN","RCIS","REG","TFEA")

list.files(paste0(work_dir,knockTF_folder[1]))
test1<-as.data.frame(fread(paste0(work_dir,knockTF_folder[3],method_name[1],"_RANK_",genes[3],".csv")))
test1<-test1[,3:ncol(test1)]
test_results<-c()
i=1
TR_cols<-sapply(strsplit(colnames(test1),"_",fixed=TRUE),function(x){return(x[[1]])})
TR_col_indexes<-which(!TR_cols%in%c(high_datasets_TR,low_datasets_TR))

for(i in TR_col_indexes){
  print(i)
  high_group_ranks<-which(test1[,i]%in%high_datasets_TR)
  low_group_ranks<-which(test1[,i]%in%low_datasets_TR)
  test_p<-wilcox.test(high_group_ranks,low_group_ranks,alternative = "less")
  test_results<-c(test_results,test_p$p.value)
}
length(test_results_adjust)
test_results_adjust<-p.adjust(test_results,"BH")
stouffer_results<-stouffer(test_results_adjust,side=1)
stouffer_results
stouffer_results$p


high_datasets_TR<-names(TR_freq_table)[1:20]
low_datasets_TR<-sample(names(TR_freq_table[which(TR_freq_table==1)]),20)
low_datasets_TR<-unique(low_list)[1:20]

#ChEA3
ChEA3<-as.data.frame(fread(paste0(work_dir,knockTF_folder[3],method_name[2],"_RANK_",genes[3],".csv")))
ChEA3<-ChEA3[,3:ncol(ChEA3)]
test_results<-c()
TR_cols<-sapply(strsplit(colnames(ChEA3),"_",fixed=TRUE),function(x){return(x[[1]])})
TR_col_indexes<-which(!TR_cols%in%c(high_datasets_TR,low_datasets_TR))

for(i in TR_col_indexes){
  print(i)
  high_group_ranks<-which(ChEA3[,i]%in%high_datasets_TR)
  low_group_ranks<-which(ChEA3[,i]%in%low_datasets_TR)
  test_p<-wilcox.test(high_group_ranks,low_group_ranks,alternative = "less")
  test_results<-c(test_results,test_p$p.value)
}

test_results_adjust<-p.adjust(test_results,"BH")
stouffer_results<-stouffer(test_results_adjust,side=1)
stouffer_results
print(stouffer_results$p)

#Lisa
Lisa<-as.data.frame(fread(paste0(work_dir,knockTF_folder[3],method_name[8],"_RANK_",genes[3],".csv")))
Lisa
Lisa<-Lisa[,3:ncol(Lisa)]
test_results<-c()
TR_cols<-sapply(strsplit(colnames(Lisa),"_",fixed=TRUE),function(x){return(x[[1]])})
TR_col_indexes<-which(!TR_cols%in%c(high_datasets_TR,low_datasets_TR))

for(i in TR_col_indexes){
  print(i)
  high_group_ranks<-match(high_datasets_TR,Lisa[,i])
  low_group_ranks<-match(low_datasets_TR,Lisa[,i])
  test_p<-wilcox.test(high_group_ranks,low_group_ranks,alternative = "less")
  test_results<-c(test_results,test_p$p.value)
}

test_results_adjust<-p.adjust(test_results,"BH")
stouffer_results<-stouffer(test_results_adjust,side=1)
stouffer_results
print(stouffer_results$p)

#CHIPATLAS
Chipatlas<-as.data.frame(fread(paste0(work_dir,knockTF_folder[3],method_name[3],"_RANK_",genes[3],".csv")))

Chipatlas<-Chipatlas[,3:ncol(Chipatlas)]
test_results<-c()
TR_cols<-sapply(strsplit(colnames(Chipatlas),"_",fixed=TRUE),function(x){return(x[[1]])})
TR_col_indexes<-which(!TR_cols%in%c(high_datasets_TR,low_datasets_TR))

for(i in TR_col_indexes){
  print(i)
  if(i%in%c(67,456)){
    next
  }
  high_group_ranks<-match(high_datasets_TR,Chipatlas[,i])
  low_group_ranks<-match(low_datasets_TR,Chipatlas[,i])
  test_p<-wilcox.test(high_group_ranks,low_group_ranks,alternative = "less")
  test_results<-c(test_results,test_p$p.value)
}

test_results_adjust<-p.adjust(test_results,"BH")
stouffer_results<-stouffer(test_results_adjust,side=1)
print(stouffer_results$p)
stouffer_results

#Cscan
Cscan<-as.data.frame(fread(paste0(work_dir,knockTF_folder[3],method_name[4],"_RANK_",genes[3],".csv")))
Cscan<-Cscan[,3:ncol(Cscan)]
test_results<-c()
TR_cols<-sapply(strsplit(colnames(Cscan),"_",fixed=TRUE),function(x){return(x[[1]])})
TR_col_indexes<-which(!TR_cols%in%c(high_datasets_TR,low_datasets_TR))

low_datasets_TR_Cscan<-names(TR_freq_table[which(TR_freq_table==1)])[names(TR_freq_table[which(TR_freq_table==1)])%in%Cscan[,1]]

for(i in TR_col_indexes){
  print(i)
  high_group_ranks<-match(high_datasets_TR,Cscan[,i])
  low_group_ranks<-match(low_datasets_TR_Cscan,Cscan[,i])
  test_p<-wilcox.test(high_group_ranks,low_group_ranks,alternative = "less")
  test_results<-c(test_results,test_p$p.value)
}

test_results_adjust<-p.adjust(test_results,"BH")
stouffer_results<-stouffer(test_results,side=1)
print(stouffer_results$p)
stouffer_results

#Enrichr
enrichr<-as.data.frame(fread(paste0(work_dir,knockTF_folder[3],method_name[5],"_RANK_",genes[3],".csv")))
enrichr<-enrichr[,3:ncol(enrichr)]
test_results<-c()
TR_cols<-sapply(strsplit(colnames(enrichr),"_",fixed=TRUE),function(x){return(x[[1]])})
TR_col_indexes<-which(!TR_cols%in%c(high_datasets_TR,low_datasets_TR))

which(high_datasets_TR%in%enrichr[,1])
which(low_datasets_TR%in%enrichr[,1])

for(i in TR_col_indexes){
  print(i)
  high_group_ranks<-match(high_datasets_TR,enrichr[,i])
  low_group_ranks<-match(low_datasets_TR,enrichr[,i])
  test_p<-wilcox.test(high_group_ranks,low_group_ranks,alternative = "less")
  test_results<-c(test_results,test_p$p.value)
}

test_results_adjust<-p.adjust(test_results,"BH")
stouffer_results<-stouffer(test_results,side=1)
print(stouffer_results$p)
stouffer_results

#MAGIC
magic<-as.data.frame(fread(paste0(work_dir,knockTF_folder[3],method_name[9],"_RANK_",genes[3],".csv")))
magic<-magic[,3:ncol(magic)]
test_results<-c()
TR_cols<-sapply(strsplit(colnames(magic),"_",fixed=TRUE),function(x){return(x[[1]])})
TR_col_indexes<-which(!TR_cols%in%c(high_datasets_TR,low_datasets_TR))

which(high_datasets_TR%in%magic[,1])
which(low_datasets_TR%in%magic[,1])

low_datasets_TR_magic<-names(TR_freq_table[which(TR_freq_table==1)])[names(TR_freq_table[which(TR_freq_table==1)])%in%magic[,1]][1:20]

for(i in TR_col_indexes){
  print(i)
  if(i%in%c(123,311)){
    next
  }
  high_group_ranks<-match(high_datasets_TR,magic[,i])
  low_group_ranks<-match(low_datasets_TR_magic,magic[,i])
  test_p<-wilcox.test(high_group_ranks,low_group_ranks,alternative = "less")
  test_results<-c(test_results,test_p$p.value)
}

test_results_adjust<-p.adjust(test_results,"BH")
# Extremely small p-values # Example of very small p-values
test_results_adjust
# Convert to Z-scores
z_scores <- qnorm(test_results_adjust, lower.tail = FALSE)
length(z_scores)

# Check for extremely large values and handle them
z_scores <- ifelse(abs(z_scores) == Inf, -37, z_scores)  # -37 is a safe lower bound for qnorm
length(z_scores)
# Calculate combined Z-score
combined_z <- sum(z_scores) / sqrt(length(z_scores))

# Convert combined Z-score back to a p-value
combined_p_value <- pnorm(combined_z, lower.tail = FALSE)

# Output
combined_z
combined_p_value

####TFEA
TFEA<-as.data.frame(fread(paste0(work_dir,knockTF_folder[3],method_name[13],"_RANK_",genes[3],".csv")))
TFEA<-TFEA[,3:ncol(TFEA)]
test_results<-c()
TR_cols<-sapply(strsplit(colnames(TFEA),"_",fixed=TRUE),function(x){return(x[[1]])})
TR_col_indexes<-which(!TR_cols%in%c(high_datasets_TR,low_datasets_TR))

for(i in TR_col_indexes){
  print(i)
  if(i%in%c(123,311)){
    next
  }
  high_group_ranks<-match(high_datasets_TR,TFEA[,i])
  low_group_ranks<-match(low_datasets_TR_TFEA,TFEA[,i])
  test_p<-wilcox.test(high_group_ranks,low_group_ranks,alternative = "less")
  test_results<-c(test_results,test_p$p.value)
}

test_results_adjust<-p.adjust(test_results,"BH")
stouffer_results<-stouffer(test_results,side=1)
print(stouffer_results$p)
stouffer_results

