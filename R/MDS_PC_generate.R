#' @title Calculate MDS of the Euclidean distance matrix
#' @description Input the reference panel and calculate the associated MDS of the Euclidean distance matrix.
#' @param data_name A file containing all GWAS file names and the corresponding locations.
#' @param batch_loc The location of the batches.
#' @param qt Whether the trait is quantitative or not. This method is only used for quantitative traits.
#' @param pcCount The number of axes of genetic variation. Default the number of PCs is 1 and the limitation of pcCount should be pcCount+envCount<n_cohort-2
#' @param ncores The number of cores which would be used for running in parallel.
#'
#' @author Siru Wang
#' @import parallel
#' @import doParallel
#' @import foreach
#' @importFrom Rdpack reprompt
#' @return null
#' @references \insertRef{magi2017trans}{env.MRmega}
#'@details The meta-regression model was proposed by \insertCite{magi2017trans;textual}{env.MRmega}.
#'@export

MR_mega_MDS<-function(data_name,batch_loc,qt=TRUE,pcCount=2,ncores){


 sum_stat_pop<-readInputfile(data_name,TRUE)
 sum_stat_pop<-check_eaf(sum_stat_pop)
 n_cohort=length(sum_stat_pop)
  marker_name_pop<-lapply(sum_stat_pop,function(ssp){
    return(unique(ssp$MARKERNAME))
  })

  snp_pop=unique(unlist(marker_name_pop))

  rm(marker_name_pop)

  num_snp=length(snp_pop)
  cat("There are ",num_snp,"different SNPs among",length(sum_stat_pop),"populations\n")


  #Creat a distance matrix
  sel_stat_snp<-sel_MafChrPos(sum_stat_pop,snp_pop)
  cat("Altogether ",dim(sel_stat_snp)[1]," good markers.\n")
  #sel_markernum=matrix(0,nrow=30,ncol=300)
  ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))

  if(ncores>1){
    cl<-makeCluster(ncores,type="FORK")#shared memory
    registerDoParallel(cl)
    sel_markernum=foreach(i_chr=1:23,.combine=rbind)%dopar%{
      sel_markernum=rep(0,300)
      sel_snp_pos=(sel_stat_snp[(sel_stat_snp$CHR==i_chr),]$POS)%/%(10^6)
      dup_last=(!duplicated(sel_snp_pos,fromLast=TRUE))
      sel_markernum[unique(sel_snp_pos,fromLast=TRUE)+1]=sel_stat_snp[(sel_stat_snp$CHR==i_chr),]$MARKERNUM[dup_last]
      return(sel_markernum)
    }
  }else{
    sel_markernum=foreach(i_chr=1:23,.combine=rbind)%do%{
      sel_markernum=rep(0,300)
      sel_snp_pos=(sel_stat_snp[(sel_stat_snp$CHR==i_chr),]$POS)%/%(10^6)
      dup_last=(!duplicated(sel_snp_pos,fromLast=TRUE))
      sel_markernum[unique(sel_snp_pos,fromLast=TRUE)+1]=sel_stat_snp[(sel_stat_snp$CHR==i_chr),]$MARKERNUM[dup_last]
      return(sel_markernum)
    }

  }


  cat("There are ",length(sel_markernum[sel_markernum!=0]),"independent variants for EAF correlation calculation.\n")

  print(dim(sel_markernum))

 for(row.idx in 1:dim(sel_markernum)[1]){

  cat("In cohort ",row.idx,"\n")
  sel_chr=(sel_markernum[row.idx,]!=0)
  cat("There are ",sum(sel_chr),"independent variants for EAF correlation calculation.\n")

 }

  exp3=paste("forDistance=cbind(",paste("sel_stat_snp$EAF",1:n_cohort,"[sel_stat_snp$MARKERNUM%in%(sel_markernum[sel_markernum!=0])]",sep="",collapse=","),")",sep="")
  eval(parse(text=exp3))

  #dim(forDistance)#25*4
  dist<-dist_pop(t(forDistance))
  dist<-double_center(dist)
  print("The dist matrix (MDS) is ")
  print(dist)

  eigen_d<-svd(dist)
  eigvalue<-ifelse(eigen_d$d>0,eigen_d$d,0)
  eigvector<-eigen_d$u
  pcs<-eigvector%*%diag(sqrt(eigvalue))[,1:pcCount]

  print("The pcs are ")
  print(pcs)

  rm(sel_stat_snp)
  rm(forDistance)
  #gc()
  #subpopulation without sex stratification
  write.table(pcs,file=paste(batch_loc,"/","pcs_spc.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(dist,file=paste(batch_loc,"/","dist_spc.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE)

  # with sex stratification
  #write.table(pcs,file=paste(batch_loc,"/","pcs_spc_sex_sim2.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE)
  #write.table(dist,file=paste(batch_loc,"/","dist_spc_sex_sim2.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE)

}




