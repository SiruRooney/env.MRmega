################################################################
#' @title read inputted file
#' @description Read the inputted of file including all names of summary statistic data across different cohorts and select useful columns which
#' would be used for meta analysis.
#' @param data_name A file involving all names of summary statistic data across different cohorts
#' @param qt Whether the trait is quantitative or not. This method is only used for quantitative traits.
#' @return  A list of of all study files.
#' @import data.table
#' @import dplyr
#' @import parallel
#' @import doParallel
#' @import foreach
#' @export

readInputfile<-function(data_name,qt=TRUE){
  sum_dat<-list()
  if(isTRUE(qt)){
    sum_stat_nd=c("MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","BETA","SE","N")
  }else{
    sum_stat_nd=c("MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","OR","OR_95L","OR_95U","N")
  }


  if (!is.null(data_name)){
    for (i in 1:length(data_name)){
      print(paste("i in filesize ",i,sep=""))

      print(paste("Reading file ",data_name[i],sep=""))

      sum_dat[[i]]=fread(data_name[i],select=sum_stat_nd)#[,pop_sum_stat%in%sum_stat_nd]
    }
  }else{
    print("data_name should be a file containing names of several summary statistics files")
  }
  names(sum_dat)=paste("Pop",1:length(data_name),sep="")
  return(sum_dat)
}

################################################################
#'@title Check EA and NEA.
#'@description Check whether EA and NEA of each genetic variant follow the same reference allele.
#'@param ssp A list of all GWAS
#'@return A list of all GWAS which would be aligned to the same reference allele.

check_eaf<-function(ssp){

  print("All GWAS are aligned to the same reference allele")
  print("Set the ancestry 1 as the reference allele")
  for(i_pop in 2:length(ssp)){
    #print(paste("Check ancestry",i_pop,sep=""))
    for(j_pop in (i_pop-1):1){
      cover=ssp[[i_pop-j_pop]]$MARKERNAME%in%ssp[[i_pop]]$MARKERNAME
      loc=match(ssp[[i_pop-j_pop]]$MARKERNAME,ssp[[i_pop]]$MARKERNAME)
      same_ea=(ssp[[i_pop]]$EA[loc[complete.cases(loc)]]!=ssp[[i_pop-j_pop]]$EA[cover])

      #print(paste(ssp[[i_pop]]$MARKERNAME[loc[complete.cases(loc)]][same_ea]," are different from the reference alleles in ancestry 1.",sep=""))
      if(any(same_ea)){
        ssp[[i_pop]]$EA[loc[complete.cases(loc)]][same_ea]=ssp[[i_pop-j_pop]]$EA[cover][same_ea]
        ssp[[i_pop]]$NEA[loc[complete.cases(loc)]][same_ea]=ssp[[i_pop-j_pop]]$NEA[cover][same_ea]
        ssp[[i_pop]]$EAF[loc[complete.cases(loc)]][same_ea]=1-ssp[[i_pop]]$EAF[loc[complete.cases(loc)]][same_ea]
        ssp[[i_pop]]$BETA[loc[complete.cases(loc)]][same_ea]=-ssp[[i_pop]]$BETA[loc[complete.cases(loc)]][same_ea]
      }
    }
  }
  return(ssp)
}

#'@title  Euclidean distance matrix
#'@description Calculate the Euclidean distance matrix.
#'@param X MAF of selected variants with MAF>5\% across all cohorts.
#'@return Output the Euclidean distance matrix.
#'
dist_pop<-function(X){

  ind_X=(!is.na(X))
  denomi<-ind_X%*%t(ind_X)
  nomi <- X %*% t(X)
  vec <- apply(X, 1, function(x) sum(x * x))
  d <- -2 * nomi + vec[col(nomi)] + vec[row(nomi)]
  diag(d) <- 0
  d <- d/denomi
  return(d)
}

#'@title Double-center distance matrix
#'@description Calculate multi-dimensional scaling of the distance matrix.
#'@param D The Euclidean distance matrix.
#'@return The multi-dimensional scaling of the Euclidean distance matrix.
double_center<-function(D){
  row_meanD=rowMeans(D)
  mean_dist=mean(row_meanD)

  D_center=(-0.5)*(D-row_meanD[col(D)]-row_meanD[row(D)]+mean_dist)

  return(D_center)
}


#'@title Select the quantified genetic variants with the specified threshold of MAF, chromosome and position.
#'@param ssp A list of all GWAS.
#'@param sp A vector including the names of the genetic variants across all cohorts.
#'@return A dataframe including the names of the genetic variants across all cohorts, the corresponding chromosomes, positions,ID and EAF across all cohorts.
#'@export
#'@author Siru Wang
sel_MafChrPos<-function(ssp,sp){

  exp=paste("sel_stat_snp=data.frame(MARKERNAME=sp,CHR=NA,POS=NA,MARKERNUM=1:length(sp),",paste("EAF",1:length(ssp),"=NA",sep ="",collapse=","),")",sep="")
  eval(parse(text=exp))
  #sel_stat_snp=data.frame(MARKERNAME=sp,CHR=NA,POS=NA,MARKERNUM=1:length(sp),EAF=NA)
  for(i_pop in 1:length(ssp)){
    ind_maf=as.matrix((ssp[[i_pop]]$EAF>0.05)&(ssp[[i_pop]]$EAF<0.95))#rownames(ind_maf)=ssp[[i_pop]]$MARKNERNAME
    ind_chr=as.matrix(ssp[[i_pop]]$CHROMOSOME>0)
    ind_pos=as.matrix(ssp[[i_pop]]$POSITION>0)
    sub_snp=sp[(sp%in%ssp[[i_pop]]$MARKERNAME)]
    loc=match(sub_snp,ssp[[i_pop]]$MARKERNAME)

    exp2=paste("sel_stat_snp$EAF",i_pop,"[(sp%in%ssp[[i_pop]]$MARKERNAME)][ind_maf[loc]]=(ssp[[i_pop]]$EAF[loc])[ind_maf[loc]]",sep="")
    eval(parse(text=exp2))
    sel_stat_snp$CHR[(sp%in%ssp[[i_pop]]$MARKERNAME)][ind_chr[loc]]=(ssp[[i_pop]]$CHROMOSOME[loc])[ind_chr[loc]]

    sel_stat_snp$POS[(sp%in%ssp[[i_pop]]$MARKERNAME)][ind_pos[loc]]=(ssp[[i_pop]]$POSITION[loc])[ind_pos[loc]]

  }

  sel=complete.cases(sel_stat_snp)
  return(sel_stat_snp[sel,])
}

#'@title Calculate the weighted MAF.
#'@param ssp A list of all GWAS.
#'@param ssf A vector including the names of the filtered genetic variants across all cohorts .
#'@param ncores The the number of cores which would be used for running in parallel.
#'@return The output is the weighted MAF
#'@author Siru Wang
getN_AverageEAF<-function(ssp,ssf,ncores){
  if(ncores>1){
    cl<-makeCluster(ncores)
    registerDoParallel(cl)
    n_eaf_pop=foreach(i_pop=1:length(ssp))%dopar%{
      ind_marker=(ssp[[i_pop]]$MARKERNAME%in%ssf)
      dtf=data.frame(MARKERNAME=ssp[[i_pop]][ind_marker,]$MARKERNAME,N=ssp[[i_pop]][ind_marker,]$N,sum_eaf=ssp[[i_pop]][ind_marker,]$N*ssp[[i_pop]][ind_marker,]$EAF)
      return(dtf)
    }
  }else{
    cl<-makeCluster(ncores)
    registerDoParallel(cl)
    n_eaf_pop=foreach(i_pop=1:length(ssp))%do%{
      ind_marker=(ssp[[i_pop]]$MARKERNAME%in%ssf)
      dtf=data.frame(MARKERNAME=ssp[[i_pop]][ind_marker,]$MARKERNAME,N=ssp[[i_pop]][ind_marker,]$N,sum_eaf=ssp[[i_pop]][ind_marker,]$N*ssp[[i_pop]][ind_marker,]$EAF)
      return(dtf)
    }
  }
  n_eaf<-data.frame(MARKERNAME=ssf,N=0,EAF_avg=0)
  for(i_pop in 1:length(ssp)){
    cover=n_eaf$MARKERNAME%in%n_eaf_pop[[i_pop]]$MARKERNAME
    loc=match(n_eaf$MARKERNAME,n_eaf_pop[[i_pop]]$MARKERNAME)
    n_eaf$N[cover]=n_eaf$N[cover]+n_eaf_pop[[i_pop]]$N[loc[complete.cases(loc)]]
    n_eaf$EAF_avg[cover]=n_eaf$EAF_avg[cover]+n_eaf_pop[[i_pop]]$sum_eaf[loc[complete.cases(loc)]]
  }

  n_eaf$EAF_avg=n_eaf$EAF_avg/n_eaf$N
  return(n_eaf)
}


#'@title Prepare data for MR-mega approach and environment-adjusted MR-MEGA approach (mr-mega)
#'
#'@param sum_stat_pop A list of all GWAS file data.
#'@param n_cohort The number of the cohorts.
#'@param usefor "mr-mega" or "env-mr-mega". The default option is "mr-mega".
#'@param qt Whether the trait is quantitative or not. This method is only used for quantitative traits.
#'@param pcCount The number of axes of genetic variation. Default the number of PCs is 1. The limitation of \code{pcCount} should be \code{pcCount}+\code{envCount}<\code{n_cohort}-2 for env-adjusted MR MEGA approach and
#'The limitation of \code{pcCount} should be \code{pcCount}<\code{n_cohort}-2 for MR MEGA approach.
#'@param calpc Logical. whether need to calculate PCs based on the inputted GWAS data.
#'@param envCount The number of environment covariates. Default the number of environment covariates is null.The limitation of envCount should be \code{pcCount}+\code{envCount}<\code{n_cohort}-2 for env-adjusted MR MEGA approach.
#'
#'@param ncores The the number of cores which would be used for running in parallel.
#' @seealso \code{\link{MR_mega_MDS}} is the function to calculate the MDS of Euclidean distance matrix and the associated axes of genetic variation.
#' @author Siru Wang
#' @importFrom data.table fwrite
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import dplyr
#' @return A list containing BETA, inverse of squared SE, PCs if required, marker information, cohort_count_filt and filtered-out genetic variants.
#'@export
MR_mega_prep<-function(sum_stat_pop,n_cohort=NULL,usefor="mr-mega",
                       qt=TRUE,pcCount=2,calpc=FALSE,envCount=NULL,ncores){

  if(is.null(sum_stat_pop)){
      print("Need to input GWAS data sets!")
  }

  if(is.null(n_cohort)) {n_cohort=length(sum_stat_pop)}
  if(ncores>1){print(paste("Parallel running and use",ncores," cores.",sep=""))}

  sum_stat_pop<-check_eaf(sum_stat_pop)

  marker_name_pop<-lapply(sum_stat_pop,function(ssp){
    return(unique(ssp$MARKERNAME))
  })

  snp_pop=unique(unlist(marker_name_pop))

  rm(marker_name_pop)

  num_snp=length(snp_pop)
  cat("There are ",num_snp,"different SNPs among",length(sum_stat_pop),"populations\n")

  if(isTRUE(calpc)){
    #Creat a distance matrix
    sel_stat_snp<-sel_MafChrPos(sum_stat_pop,snp_pop)
    cat("Altogether ",dim(sel_stat_snp)[1]," good markers.\n")
    #sel_markernum=matrix(0,nrow=30,ncol=300)

    ncores <- min(c(ncores, detectCores(logical = TRUE)))

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

    exp3=paste("forDistance=cbind(",paste("sel_stat_snp$EAF",1:n_cohort,"[sel_stat_snp$MARKERNUM%in%(sel_markernum[sel_markernum!=0])]",sep="",collapse=","),")",sep="")
    eval(parse(text=exp3))

    dist<-dist_pop(t(forDistance))
    dist<-double_center(dist)
    print("The dist matrix (MDS) is ")
    print(dist)

    eigen_d<-svd(dist)
    eigvalue<-ifelse(eigen_d$d>0,eigen_d$d,0)
    eigvector<-eigen_d$u
    pcs<-eigvector%*%diag(sqrt(eigvalue))[,1:pcCount]

    rm(sel_stat_snp)
    rm(forDistance)
    #gc()
    #write.table(pcs,file=paste0(batch_loc,"/pcs.txt"),sep="\t",row.names=FALSE,col.names=FALSE)
  }

  marker_cohort_count=array(0,dim=c(num_snp,1),dimnames=list(snp_pop,"Count_cohorts"))
  eaf_filt=list()

  for (i_pop in 1:n_cohort){
    sub_snp=snp_pop[(snp_pop%in%sum_stat_pop[[i_pop]]$MARKERNAME)]
    eaf_filt[[i_pop]]=(sum_stat_pop[[i_pop]]$EAF!=-1)
    marker_cohort_count[(snp_pop%in%(sum_stat_pop[[i_pop]]$MARKERNAME[eaf_filt[[i_pop]]]))]=marker_cohort_count[(snp_pop%in%(sum_stat_pop[[i_pop]]$MARKERNAME[eaf_filt[[i_pop]]]))]+1
  }
  rm(eaf_filt)
  rm(sub_snp)
  gc()

  #n_avgeaf<-getN_AverageEAF(sum_stat_pop,snp_pop)

  #filtering the SNPs
  if(is.null(envCount)){
    snp_pop_filt=rownames(marker_cohort_count)[marker_cohort_count-2>pcCount]
    cohort_count_filt=data.frame(MARKERNAME=snp_pop_filt,Count_cohorts=marker_cohort_count[marker_cohort_count-2>pcCount])
    if(all(marker_cohort_count-2>pcCount)){
      logout=NULL
    }else{
      exp=paste("logout=data.frame(MARKERNAME=rownames(marker_cohort_count)[(marker_cohort_count-2)<=pcCount],",paste(rep(c("beta","se"),(pcCount+1)),rep(0:(pcCount*2),each=2),"=NA",collapse=",",sep=""),",chisq_association=NA,ndf_association=NA,pvalue_association=NA,chisq_heter=NA,ndf_chisq_heter=NA,pvalue_heter=NA,chisq_residual=NA,ndf_residual=NA,pvalue_residual=NA,comments=\"SmallcohortCount\")",collapse=",",sep="")
      eval(parse(text=exp))
    }
  }else{
    snp_pop_filt=rownames(marker_cohort_count)[marker_cohort_count-2>(pcCount+envCount)]
    cohort_count_filt=data.frame(MARKERNAME=snp_pop_filt,Count_cohorts=marker_cohort_count[marker_cohort_count-2>(pcCount+envCount)])
    if(all(marker_cohort_count-2>(pcCount+envCount))){
      logout=NULL
    }else{
      exp=paste("logout=data.frame(MARKERNAME=rownames(marker_cohort_count)[(marker_cohort_count-2)<=(pcCount+envCount)],",paste(rep(c("beta","se"),((pcCount+envCount)+1)),rep(0:((pcCount+envCount)*2),each=2),"=NA",collapse=",",sep=""),",chisq_association=NA,ndf_association=NA,pvalue_association=NA,chisq_heter=NA,ndf_chisq_heter=NA,pvalue_heter=NA,chisq_residual=NA,ndf_residual=NA,pvalue_residual=NA,comments=\"SmallcohortCount\")",collapse=",",sep="")
      eval(parse(text=exp))
    }
  }

  for (i_pop in 1:n_cohort){
    filt_threshold=(sum_stat_pop[[i_pop]]$MARKERNAME%in%snp_pop_filt)
    sum_stat_pop[[i_pop]]=sum_stat_pop[[i_pop]][filt_threshold,]
  }

  rm(filt_threshold)
  n_avgeaf<-getN_AverageEAF(sum_stat_pop,snp_pop_filt,ncores)
  #########################################################################
  exp=paste("beta_pop=data.frame(MARKERNAME=snp_pop_filt,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
  eval(parse(text=exp))
  exp2=paste("invse2_pop=data.frame(MARKERNAME=snp_pop_filt,",paste("Pop",1:n_cohort,"=NA",sep ="",collapse=","),")",sep="")
  eval(parse(text=exp2))
  exp3=paste("marker_inf_pop=data.frame(MARKERNAME=snp_pop_filt,",paste(c("CHROMOSOME","EA","NEA"),"=NA",sep="",collapse=","),")",sep="")
  eval(parse(text=exp3))
  Dir_pop=data.frame(MARKERNAME=snp_pop_filt,dir="?")
  #beta_pop=data.frame(MARKERNAME=snp_pop_filt)
  #invse2_pop=data.frame(MARKERNAME=snp_pop_filt)
  #marker_inf_pop=data.frame(MARKERNAME=snp_pop_filt)

  for(i_pop in 1:n_cohort){
    cover=beta_pop$MARKERNAME%in%sum_stat_pop[[i_pop]]$MARKERNAME
    loc=match(beta_pop$MARKERNAME,sum_stat_pop[[i_pop]]$MARKERNAME)
    beta_pop[,i_pop+1][cover]=sum_stat_pop[[i_pop]]$BETA[loc[complete.cases(loc)]]
    invse2_pop[,i_pop+1][cover]=sum_stat_pop[[i_pop]]$SE[loc[complete.cases(loc)]]
    marker_inf_pop[cover,-1]=sum_stat_pop[[i_pop]][,c("CHROMOSOME","EA","NEA")][loc[complete.cases(loc)],]
    if(any(is.numeric(sum_stat_pop[[i_pop]]$BETA[loc[complete.cases(loc)]]))){
      Dir_pop$dir[cover][which(sum_stat_pop[[i_pop]]$BETA[loc[complete.cases(loc)]]>0)]=paste0(Dir_pop$dir[cover][which(sum_stat_pop[[i_pop]]$BETA[loc[complete.cases(loc)]]>0)],"+")
      Dir_pop$dir[cover][which(sum_stat_pop[[i_pop]]$BETA[loc[complete.cases(loc)]]<=0)]=paste0(Dir_pop$dir[cover][which(sum_stat_pop[[i_pop]]$BETA[loc[complete.cases(loc)]]<=0)],"-")
    }
  }

  marker_inf_pop=left_join(marker_inf_pop,n_avgeaf,by="MARKERNAME")
  marker_inf_pop=left_join(marker_inf_pop,Dir_pop,by="MARKERNAME")


  #rm(cover)
  invse2_pop[-1]=round(invse2_pop[-1]^{-2})

  if(isTRUE(calpc)){
    gwas_input<-list(
      beta_pop=beta_pop,
      invse2_pop=invse2_pop,
      marker_inf_pop=marker_inf_pop,
      pcs=pcs,
      cohort_count_filt=cohort_count_filt,
      sel_gene=logout
    )
  }else{
    gwas_input<-list(
      beta_pop=beta_pop,
      invse2_pop=invse2_pop,
      marker_inf_pop=marker_inf_pop,
      cohort_count_filt=cohort_count_filt,
      sel_gene=logout
    )
  }

    return(gwas_input)
  #}
}









