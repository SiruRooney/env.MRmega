#'@title The weighted environment-adjusted linear regression.
#'@description Fit the linear regression with weights using GWAS data across cohorts.
#'The contribution of each GWAS is weighted by the estimated inverse variance of the reference allele effect at the corresponding variant.
#'@param Y BETA values of the genetic variants across all cohorts.
#'@param X PCs plus environment covariates.
#'@param W The vector of weights, each component referring to the estimated inverse variance of the allele effect at the corresponding variants. The length of the vector is the number of cohorts.
#'@param n_pc the number of PCs
#'@return Output the data frame containing information of fitting the regression model: the estimated coefficients, standard errors, the deviance of the different constrained  models and the corresponding degrees of freedom.
#'@author Siru Wang
#'@export
lr_env_w<-function(Y,X,W,n_pc){
  X=as.matrix(X)
  if((length(Y)-dim(X)[2])<1){
    df_lr_w=NULL
 }else{
  m0<-lm(Y~1,weights=W)
  m1<-lm(Y~X,weights=W)

  m_env<-lm(Y~X[,-(((dim(X)[2])-n_pc+1):(dim(X)[2]))],weights=W)
  m_pc<-lm(Y~X[,(((dim(X)[2])-n_pc+1):(dim(X)[2]))],weights=W)

  RSS=sum(W*(m1$residuals^2))
  TSS0=sum(W*Y^2)
  TSS=sum(W*m0$residuals^2)

  RSS_env=sum(W*(m_env$residuals^2))
  RSS_pc=sum(W*(m_pc$residuals^2))

  RES_each<-apply(as.matrix(1:dim(X)[2]),1,function(idx_col){
    m_each<-lm(Y~X[,-idx_col],weights=W)
    return(m_each$residuals)
  })

  RSS_each=colSums(W*RES_each^2)
  names(RSS_each)=colnames(X)

  if(all(!is.na(m1$coefficients))){
    sum_m1<-summary(m1)
    coef_se<-c(t(sum_m1$coefficients[,1:2]))
  }else{
    coef_se<-rep(NA,2*(dim(X)[2]+1))
    sum_m1<-summary(m1)
    coef_se[seq(1,2*(dim(X)[2]+1),2)]<-m1$coefficients
    coef_se[seq(2,2*(dim(X)[2]+1),2)][!is.na(m1$coefficients)]<-sum_m1$coefficients[,2]
  }
  names(coef_se)=sprintf(rep(c("beta%s","se%s"),dim(X)[2]+1),rep(0:dim(X)[2],each=2))

  coef.se=as.data.frame(as.list(coef_se),col.names=names(coef_se))
  RSS.each=as.data.frame(as.list(RSS_each),col.names=names(RSS_each))
  df_lr_w=data.frame(coef.se,RSS=RSS,df_RSS=m1$df.residual,TSS0=TSS0,TSS=TSS,RSS_env=RSS_env,RSS_pc=RSS_pc,RSS.each)
 }
  return(df_lr_w)
}


#'@title Linear regression model with weights for environment-adjusted MR-MEGA approach
#'@description Fit the environment-adjusted meta-regression model using the normalized statistic data across all cohorts.
#'@param dta A list containing BETA, SE^{-2},PCs and envs, cohort_count_filt which may be used for Bayesian factor.
#'@param qt Whether the trait is quantitative or not. This method is only used for quantitative traits.
#'@param ncores The the number of cores which would be used for running in parallel.
#'@param pcCount The number of axes of genetic variation.
#'@returns Output a file containing names of genetic variants, estimated coefficients, standard errors, chisq value of the association, the number of degrees of freedom of the association, p-value of the association, chisq value of the heterogeneity due to different ancestry, ndf of the heterogeneity due to different ancestry, p-value of the heterogeneity due to different ancestry, chisq value of the residual heterogeneity, ndf of the residual heterogeneity,  p-value of the residual heterogeneity.
#'
#'@author Siru Wang
#'@export
#'@import data.table
#'@import parallel
#'@import doParallel
#'@import foreach
#'@import dplyr
MR_mega_env_lrw<-function(dta,qt=TRUE,ncores,pcCount){
     #print(names(dta))
      beta_pop<-dta$beta_pop
      invse2_pop<-dta$invse2_pop
      pcs<-data.matrix(dta$pcs)
      colnames(pcs)=paste0("pc",1:pcCount)
      envs<-data.matrix(dta$env)
      colnames(envs)=paste0("env",1:dim(envs)[2])


  covs=cbind(envs,pcs)
  #Run linear regression in parallel and return results of hypothesis testings
  ncores <- min(c(ncores,detectCores(logical = TRUE)))
  print("Linear regression model for each genetic varaiant")
  if(ncores>1){
    cl<-makeCluster(ncores,type="FORK")#shared memory
    registerDoParallel(cl)
    lr_out<-foreach(i_snp=1:dim(beta_pop)[1],.combine=rbind)%dopar%{
     check_w=!is.na(invse2_pop[i_snp,-1])
      lr_out<-lr_env_w(unlist(beta_pop[i_snp,-1])[check_w],as.matrix(covs)[check_w,],unlist(invse2_pop[i_snp,-1])[check_w],pcCount)
      return(lr_out)
    }
  }else{
    lr_out<-foreach(i_snp=1:dim(beta_pop)[1],.combine=rbind)%do%{
      cat("i_snp=",i_snp,"\n")
      check_w=!is.na(invse2_pop[i_snp,-1])
      lr_out<-lr_env_w(unlist(beta_pop[i_snp,-1])[check_w],as.matrix(covs)[check_w,],unlist(invse2_pop[i_snp,-1])[check_w],pcCount)
      return(lr_out)
    }
  }

  TSS_RSS=abs(lr_out$TSS-lr_out$RSS)
  RSS=abs(lr_out$RSS)
  TSS0_RSS=abs(lr_out$TSS0-lr_out$RSS)
  RSS_env_RSS=abs(lr_out$RSS_env-lr_out$RSS)
  RSS_pc_RSS=abs(lr_out$RSS_pc-lr_out$RSS)


  pModelHet=pchisq(TSS_RSS,dim(covs)[2],lower.tail = FALSE)
  #df_pModelHet=dim(covs)[2]
  pResidHet=pchisq(RSS,lr_out$df_RSS,lower.tail = FALSE)
  #df_pResidHet=lr_out$df_RSS
  pModelTest=pchisq(TSS0_RSS,dim(covs)[2]+1,lower.tail = FALSE)
  #df_pModelTest=dim(covs)[2]+1
  pResEnv=pchisq(RSS_env_RSS,pcCount,lower.tail = FALSE)
  pResPC=pchisq(RSS_pc_RSS,dim(envs)[2],lower.tail = FALSE)

  RSS_each_RSS=abs(as.matrix(lr_out[,-(1:((dim(covs)[2]+1)*2+6))])-lr_out$RSS)
  #RSS_each_RSS=abs(as.matrix(lr_out[,-(1:((dim(covs)[2]+1)*2+4))])-lr_out$RSS)

  pModelHet_each=apply(RSS_each_RSS,2,function(rer){pchisq(rer,1,lower.tail = FALSE)})
  colnames(pModelHet_each)=paste(rep("pvalue_",dim(covs)[2]),colnames(RSS_each_RSS),sep="")

  #logBF=rep(NA,length(cohort_count_filt))
  #logBF[(cohort_count_filt-2>dim(covs)[2])]=0.5*(TSS0_RSS[(cohort_count_filt-2>dim(covs)[2])]-(dim(covs)[2]+1)*log(cohort_count_filt[cohort_count_filt-2>dim(covs)[2]]))

  logout=data.frame(MARKERNAME=beta_pop$MARKERNAME,
                    lr_out[,1:((dim(covs)[2]+1)*2)],
                    chisq_association=TSS0_RSS,
                    ndf_association=dim(covs)[2]+1,
                    pvalue_association=pModelTest,
                    chisq_heter=TSS_RSS,
                    ndf_heter=dim(covs)[2],
                    pvalue_heter=pModelHet,
                    chisq_residual=RSS,
                    ndf_residual=lr_out$df_RSS,
                    pvalue_residual=pResidHet,
                    chisq_env=RSS_pc_RSS,
		    ndf_env=dim(envs)[2],
		    pvalue_env=pResPC,
		    chisq_PC=RSS_env_RSS,
		    ndf_PC=pcCount,
		    pvalue_PC=pResEnv)

  logout<-left_join(dta$marker_inf_pop,logout,by="MARKERNAME")
  logout<-left_join(logout,dta$cohort_count_filt,by="MARKERNAME")

  return(logout)
}


#'@title The weighted linear regression.
#'@description Fit the linear regression with weights using GWAS data across cohorts.
#'The contribution of each GWAS is weighted by the estimated inverse variance of the reference allele effect at the corresponding variant.
#'@param Y BETA values of the genetic variants across all cohorts.
#'@param X PCs.
#'@param W The vector of weights, each component referring to the estimated inverse variance of the allele effect at the corresponding variants. The length of the vector is the number of cohorts.
#'@return Output the data frame containing information of fitting the regression model: the estimated coefficients, standard errors, the deviance of the different constrained  models and the corresponding degrees of freedom.
#'@author Siru Wang
#'@export
lr_w<-function(Y,X,W){
  X=as.matrix(X)
  if((length(Y)-dim(X)[2])<1){
    df_lr_w=NULL
  }else{
    m0<-lm(Y~1,weights=W)
    m1<-lm(Y~X,weights=W)

    RSS=sum(W*(m1$residuals^2))
    TSS0=sum(W*Y^2)
    TSS=sum(W*m0$residuals^2)


    if(dim(X)[2]!=1){
      RES_each<-apply(as.matrix(1:dim(X)[2]),1,function(idx_col,X){
        m_each<-lm(Y~X[,-idx_col],weights=W)
        return(m_each$residuals)
      },X)
    }else{
      RES_each=NULL
    }

    if(!is.null(RES_each)){
      RSS_each=colSums(W*RES_each^2)
      names(RSS_each)=sprintf(rep("non_cov%s",length(RSS_each)),1:length(RSS_each))
    }else{
      RSS_each=TSS
      names(RSS_each)=sprintf(rep("non_cov%s",length(RSS_each)),1:length(RSS_each))
    }

    sum_m1<-summary(m1)
    coef_se<-c(t(sum_m1$coefficients[,1:2]))
    names(coef_se)=sprintf(rep(c("beta%s","se%s"),NROW(sum_m1$coefficients)),rep(0:(NROW(sum_m1$coefficients)-1),each=2))
    coef.se=as.data.frame(as.list(coef_se),col.names=names(coef_se))
    RSS.each=as.data.frame(as.list(RSS_each),col.names=names(RSS_each))
    df_lr_w=data.frame(coef.se,RSS=RSS,df_RSS=m1$df.residual,TSS0=TSS0,TSS=TSS,RSS.each)
  }
  return(df_lr_w)
}

#'@title Linear regression model with weights for MR-MEGA approach
#'@references
#'\insertRef{magi2017trans}{env.MRmega}
#'@details The meta-regression model was proposed by \insertCite{magi2017trans;textual}{env.MRmega}.
#'@description Fit the meta-regression model using the normalized statistic data across all cohorts.
#'@param dta A list containing BETA, SE^{-2},PCs, cohort_count_filt which may be used for Bayesian factor.
#'@param qt Whether the trait is quantitative or not. This method is only used for quantitative traits.
#'@param ncores The the number of cores which would be used for running in parallel.
#'@param pcCount The number of axes of genetic variation.
#'@return Output a file containing names of genetic variants, estimated coefficients, standard errors, chisq value of the association, the number of degrees of freedom of the association, p-value of the association, chisq value of the heterogeneity due to different ancestry, ndf of the heterogeneity due to different ancestry, p-value of the heterogeneity due to different ancestry, chisq value of the residual heterogeneity, ndf of the residual heterogeneity,  p-value of the residual heterogeneity.
#'@author Siru Wang
#'@export
#'@import data.table
#'@import parallel
#'@import doParallel
#'@import foreach
#'@import dplyr
#'@importFrom Rdpack reprompt
MR_mega_lrw<-function(dta,qt=TRUE,ncores,pcCount){
  #print(names(dta))
  beta_pop<-dta$beta_pop
  invse2_pop<-dta$invse2_pop
  pcs<-data.matrix(dta$pcs)
  colnames(pcs)=paste0("pc",1:pcCount)

  #Run linear regression in parallel and return results of hypotheses testings
  ncores <- min(c(ncores, detectCores(logical = TRUE)))
  print("Linear regression model for each genetic varaiant")
  if(ncores>1){
    cl<-makeCluster(ncores,type="FORK")#shared memory
    registerDoParallel(cl)
    lr_out<-foreach(i_snp=1:dim(beta_pop)[1],.combine=rbind)%dopar%{
      #cat("i_snp=",i_snp,"\n")
      check_w=!is.na(unlist(invse2_pop[i_snp,-1]))
      if(pcCount>1){
        lr_out<-lr_w(unlist(beta_pop[i_snp,-1])[check_w],as.matrix(pcs)[check_w,],unlist(invse2_pop[i_snp,-1])[check_w])
      }else{
        lr_out<-lr_w(unlist(beta_pop[i_snp,-1])[check_w],as.matrix(pcs)[check_w],unlist(invse2_pop[i_snp,-1])[check_w])
      }
      return(lr_out)
    }
  }else{
    lr_out<-foreach(i_snp=1:dim(beta_pop)[1],.combine=rbind)%do%{
      cat("i_snp=",i_snp,"\n")
      check_w=!is.na(unlist(invse2_pop[i_snp,-1]))
      if(pcCount>1){
        lr_out<-lr_w(unlist(beta_pop[i_snp,-1])[check_w],as.matrix(pcs)[check_w,],unlist(invse2_pop[i_snp,-1])[check_w])
      }else{
        lr_out<-lr_w(unlist(beta_pop[i_snp,-1])[check_w],as.matrix(pcs)[check_w],unlist(invse2_pop[i_snp,-1])[check_w])
      }
      return(lr_out)
    }
  }

  TSS_RSS=abs(lr_out$TSS-lr_out$RSS)
  RSS=abs(lr_out$RSS)
  TSS0_RSS=abs(lr_out$TSS0-lr_out$RSS)

  pModelHet=pchisq(TSS_RSS,pcCount,lower.tail = FALSE)
  #df_pModelHet=pcCount
  pResidHet=pchisq(RSS,lr_out$df_RSS,lower.tail = FALSE)
  #df_pResidHet=lr_out$df_RSS
  pModelTest=pchisq(TSS0_RSS,pcCount+1,lower.tail = FALSE)
  #df_pModelTest=pcCount+1
  if(pcCount==1){
    RSS_each_RSS=TSS_RSS
    pModelHet_each=pModelHet
  }else{
    RSS_each_RSS=abs(as.matrix(lr_out[,-(1:((pcCount+1)*2+4))])-RSS)
    pModelHet_each=apply(RSS_each_RSS,2,function(rer){pchisq(rer,1,lower.tail = FALSE)})
  }
  #pModelHet_each=apply(RSS_each_RSS)

  #logBF=0.5*(TSS0_RSS-(pcCount+1)*log(cohort_count_filt[cohort_count_filt-2>pcCount]))
  logout<-data.frame(MARKERNAME=beta_pop$MARKERNAME,
                    lr_out[,1:((pcCount+1)*2)],
                    chisq_association=TSS0_RSS,
                    ndf_association=pcCount+1,
                    pvalue_association=pModelTest,
                    chisq_heter=TSS_RSS,
                    ndf_heter=pcCount,
                    pvalue_heter=pModelHet,
                    chisq_residual=RSS,
                    ndf_residual=lr_out$df_RSS,
                   pvalue_residual=pResidHet)

  logout<-left_join(dta$marker_inf_pop,logout,by="MARKERNAME")
  logout<-left_join(logout,dta$cohort_count_filt,by="MARKERNAME")
  return(logout)
}


