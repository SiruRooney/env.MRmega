##########################################################
#' @title Check unique
#' @description Function: Check whether there are some genetic variants shown more than once in each cohort.
#' @param ssp A dataframe of one GWAS
#' @param ncores The number of cores which would be used for running in parallel.
#' @return A list of all study files in which genetic variants shown more than once would be filtered out.
#' @author Siru Wang
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import iterators
#' @import dplyr
check_uniq<-function(ssp,ncores){
  print("Check whether there are some markers shown more than once in each ancestry")

    rep_pop=duplicated(ssp$MARKERNAME)
    print(paste("There are ",length(unique(ssp$MARKERNAME[rep_pop]))," marker shown more than once in the ancestry",sep=""))
    rep_marker=as.matrix(unique(ssp$MARKERNAME[rep_pop]))
    if(length(unique(rep_pop))!=1){
      ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
      cl<-makeCluster(ncores,type="FORK")#shared memory
      registerDoParallel(cl)
      uniq_marker=foreach(i_rep=iter(rep_marker,by="row"),.combine=rbind)%dopar%{
        cat("i_rep",i_rep,"\n")

        filt=ssp%>%filter(MARKERNAME%in%i_rep)%>%filter(N==max(N))
        if(dim(filt)[1]>1){
          filt=filt[1,]
        }
        return(filt)
      }
      ssp=ssp%>%distinct(MARKERNAME,.keep_all=TRUE)
      ssp[ssp$MARKERNAME%in%rep_marker,]=uniq_marker
    }
    print(paste("There are ",dim(ssp)[1]," unique markers in the ancestry.",sep=""))

  return(ssp)
}


#############################################################################
#'@title MR_mega_normalize
#'@description Normalize the inputted file of summary statistic data across different cohorts, which includes
#'deleting all genetic variants with NA's and filtering out genetic variants shown more once.
#'@param data_name A vector of data name indicating the specific location of the summary statistic data of one GWAS.
#'@param qt Whether the trait is quantitative or not. This method is only used for quantitative traits.
#'@param ncores The number of cores which would be used for running in parallel.
#'@return Output the normalized summary statistic data across different cohorts
#'@export
#' @author Siru Wang
#' @import parallel
#' @import doParallel
#' @import foreach
MR_mega_normalize<-function(data_name,qt=TRUE,ncores){

  #sum_stat$Pop_i is data.frame
  if(isTRUE(qt)){
    sum_stat_nd=c("MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","BETA","SE","N")
  }else{
    sum_stat_nd=c("MARKERNAME","CHROMOSOME","POSITION","EA","NEA","EAF","OR","OR_95L","OR_95U","N")
  }
  print(paste("Reading file ",data_name,sep=""))
  sum_stat_pop=fread(data_name,select=sum_stat_nd)


  if(ncores>1){print(paste("Parallel running and use ",ncores," cores.",sep=""))}

  #Delete the genetic variants of which eaf is NA and MARKERNAME is NA
  #check_na<-complete.cases(sum_stat_pop)
  #sum_stat_pop<-sum_stat_pop[check_na,]

 check_eaf<-complete.cases(sum_stat_pop$EAF)
 check_marker<-complete.cases(sum_stat_pop$MARKERNAME)
 check_na<-(check_eaf)&(check_marker)
 sum_stat_pop<-sum_stat_pop[check_na,]

  print(dim(sum_stat_pop))
  sum_stat_pop<-check_uniq(sum_stat_pop,ncores)

 fwrite(sum_stat_pop,file=data_name,sep="\t")
}




