# R source code 
# Main function for integrative models with no CNV-methylation association

# Arguments to input in the function "Same_Level_Original" :

# -Gene = matrix of genes, number of columns = number of genes
# -methy = matrix of methylation sites, number of columns = number of methylation sites
# -CNV = matrix of CNVs, number of columns = number of CNVs
# -y = matrix containing survival time (y), censoring indicator (Event), and other clinical variables
# -Gene_CNV = matrix of three columns. The first column contains the CNV identifiers,  
#the second column the gene to which the CNV maps, and the third column the associated
#functional network
# -Gene_Methy has the same format as Gene_CNV expect that the first column contains
#the identifiers for methylations sites
# -Gene_Pathway = matrix of two columns. First column contains the gene identifiers of genes and
#the second column gives the associated functional network
# -multi_methy = TRUE means that a multivariate model is used for the association between 
#methylation sites and CNVs, if FALSE univariate models are fitted
# -intra = TRUE means that an Integrative-gene scenario is considered
# -pathway = TRUE means that an Integrative-network scenario is considered
# -nfolds = number of folds used in the cross-validation procedure of the glmnet function called in  
#the function "regression" defined in UtilFunctions.R
# -alpha = penalty to use in the glmnet function, alpha = 1 corresponds to the lasso penalty

 Same_Level_Original <- function(Gene,methy,CNV,y,Gene_CNV,Gene_Methy,
                                multi_methy= FALSE, intra =FALSE, pathway=FALSE,
                                nfolds=10,alpha =1)
{
  n = nrow(y) 
   
  geneWithBoth = colnames(Gene)[which(colnames(Gene)%in%intersect(colnames(Gene),intersect(Gene_Methy[,2],Gene_CNV[,2])))]
  
  geneWithCNV = colnames(Gene)[which(colnames(Gene)%in%intersect(colnames(Gene),Gene_CNV[,2]))]
  geneWithCNV = geneWithCNV[-which(geneWithCNV%in%intersect(geneWithBoth,geneWithCNV))]
  
  geneWithMethy = colnames(Gene)[which(colnames(Gene)%in%intersect(colnames(Gene),Gene_Methy[,2]))]
  geneWithMethy = geneWithMethy[-which(geneWithMethy%in%intersect(geneWithBoth,geneWithMethy))]
  
  geneWithout = (1:ncol(Gene))[-which(colnames(Gene)%in%c(geneWithBoth,geneWithCNV,geneWithMethy))]
  
  RNA_modulated_methy = list() 
  RNA_modulated_CNV= list()
  RNA_modulated_other= list()
  
######################################################################
##########   1. Analysis for genes with methylation and CNV data
######################################################################  
  res3 =list()
  for( g in geneWithBoth){
    if (intra) {sub_CNV = as.character(Gene_CNV[Gene_CNV[,2] == g,1] )
    }else{ 
      if (pathway){
        path = Gene_Pathway[which(Gene_Pathway[,1]%in%g),2]
        gene_in_path = as.character(Gene_Pathway[which(Gene_Pathway[,2]%in%path),1])
        sub_CNV = as.character(Gene_CNV[which(Gene_CNV[,2]%in%gene_in_path),1])
      }else{      
        sub_CNV = colnames(CNV)}
    }    
    ind_CNV =  which(colnames(CNV)%in%sub_CNV)
    if (length(ind_CNV)!=0){ 
      part_CNV = as.matrix(CNV[,ind_CNV])
      colnames(part_CNV) = colnames(CNV)[ind_CNV]
      X = part_CNV
    }else{
      X = NULL
    }
    
    sub_methy = as.character(Gene_Methy[Gene_Methy[,2] == g,1] )
    ind_methy =  which(colnames(methy)%in%sub_methy)
    part_methy = as.matrix(methy[,sub_methy])
    colnames(part_methy) = colnames(methy)[ind_methy]
    X= cbind(part_methy,X)
    
    y_tmp = as.matrix(Gene[,g])
    colnames(y_tmp) = g
    res3[[g]] = regression(x=X,y=y_tmp,alpha=1, nfolds =nfolds)
    
    
    select = res3[[g]]$active_index 
    
    ind = which(colnames(part_methy)%in%names(select))
    ind_tmp =  which(names(select)%in%colnames(part_methy))
    if ((length(ind)!=1)&(length(select)!=0))
    	{lol= part_methy[,ind]%*%(res3[[g]]$active[ind_tmp])}
    if (length(ind)==1){lol= part_methy[,ind]*res3[[g]]$active[ind_tmp]}
    if (length(ind)==0){lol = as.matrix(methy[,1]*0)}
    RNA_modulated_methy[[g]] = lol 
    

    ind = which(colnames(part_CNV)%in%names(select))
    ind_tmp =  which(names(select)%in%colnames(part_CNV))
    if ((length(ind)!=1)&(length(ind)!=0))
    	{lol= part_CNV[,ind]%*%(res3[[g]]$active[ind_tmp])}
    if (length(ind)==1){lol= part_CNV[,ind]*res3[[g]]$active[ind_tmp]}
    if (length(ind)==0){lol =  as.matrix(CNV[,1]*c(0))}
    RNA_modulated_CNV[[g]] = lol

    
    RNA_modulated_other[[g]]= Gene[,g] - RNA_modulated_methy[[g]] - 
    	RNA_modulated_CNV[[g]] - res3[[g]]$inter
    
  }
  
  A1 = matrix(unlist(RNA_modulated_methy),
              ncol=length(RNA_modulated_methy),
              byrow=FALSE)
  colnames(A1) = paste(geneWithBoth,"methy",sep="_")
  
  A2 = matrix(unlist(RNA_modulated_CNV),
              ncol=length(RNA_modulated_CNV),
              byrow=FALSE)
  colnames(A2) = paste(geneWithBoth,"CNV",sep="_")

  A3 = matrix(unlist(RNA_modulated_other),
              ncol=length(RNA_modulated_other),
              byrow=FALSE)
  colnames(A3) = paste(geneWithBoth,"other",sep="_")
  
######################################################################
##########   2. Analysis for genes with methylation data   #######
######################################################################  
  RNA_modulated_CNV = list()
  RNA_modulated_methy = list()
  RNA_modulated_other =list() 
  resB=list()
  if (length(geneWithMethy)!=0){
    for( g in geneWithMethy){
      if (intra) {sub_CNV = as.character(Gene_CNV[Gene_CNV[,2] == g,1] )
      }else{ 
        if (pathway){
          path = Gene_Pathway[which(Gene_Pathway[,1]%in%g),2]
          gene_in_path = as.character(Gene_Pathway[which(Gene_Pathway[,2]%in%path),1])
          sub_CNV = as.character(Gene_CNV[which(Gene_CNV[,2]%in%gene_in_path),1])
        }else{      
          sub_CNV = colnames(CNV)}
      }
      ind_CNV =  which(colnames(CNV)%in%sub_CNV)
      if (length(ind_CNV)!=0){ 
        part_CNV = as.matrix(CNV[,ind_CNV])
        colnames(part_CNV) = colnames(CNV)[ind_CNV]
        X = part_CNV
      }else{
        X = NULL
      }
      
      sub_methy = as.character(Gene_Methy[Gene_Methy[,2] == g,1] )
      ind_methy =  which(colnames(methy)%in%sub_methy)
      part_methy = as.matrix(methy[,sub_methy])
      colnames(part_methy) = colnames(methy)[ind_methy]
      X= cbind(part_methy,X)
      
      y_tmp = as.matrix(Gene[,g])
      colnames(y_tmp) = g
      resB[[g]] = regression(x=X,y=y_tmp)
      
      select = resB[[g]]$active_index 
      
      ind = which(colnames(part_methy)%in%names(select))
      ind_tmp =  which(names(select)%in%colnames(part_methy))
      if ((length(ind)!=1)&(length(select)!=0))
      	{lol= part_methy[,ind]%*%(resB[[g]]$active[ind_tmp])}
      if (length(ind)==1){lol= part_methy[,ind]*resB[[g]]$active[ind_tmp]}
      if (length(ind)==0){lol = as.matrix(methy[,1]*0)}
      RNA_modulated_methy[[g]] = lol   
      
      if (length(ind_CNV)!=0){
        ind = which(colnames(part_CNV)%in%names(select))
        ind_tmp =  which(names(select)%in%colnames(part_CNV))
        if ((length(ind)!=1)&(length(select)!=0))
        		{lol= part_CNV[,ind]%*%(resB[[g]]$active[ind_tmp])}
        if (length(ind)==1){lol= part_CNV[,ind]*resB[[g]]$active[ind_tmp]}
        if (length(ind)==0){lol = as.matrix(CNV[,1]*0)}
        RNA_modulated_CNV[[g]] = lol  
      }else{
        RNA_modulated_CNV[[g]] = as.matrix(CNV[,1]*c(0))
      }
      
      RNA_modulated_other[[g]]= Gene[,g]-RNA_modulated_methy[[g]]-
        RNA_modulated_CNV[[g]]-resB[[g]]$inter
    }
    
    B1 = matrix(unlist(RNA_modulated_methy),ncol=length(RNA_modulated_methy),byrow=FALSE)
    colnames(B1) = paste(names(resB),"Methy",sep="_")
    
    B2 = matrix(unlist(RNA_modulated_CNV),ncol=length(RNA_modulated_CNV),byrow=FALSE)
    colnames(B2) = paste(names(resB),"CNV",sep="_")
    
    B3 =  matrix(unlist(RNA_modulated_other),ncol=length(RNA_modulated_other),byrow=FALSE)
    colnames(B3) = paste(names(resB),"other",sep="_")
  }else{
    B1 = rep(0,n)
    B2 = B1
    B3 = B1 
  }
  
############################################################
##########   3. Analysis for genes with CNV   ##########
############################################################  
  
  RNA_modulated_CNV = list()
  RNA_modulated_other = list()
  
  resC=list()
  if(length( c(geneWithCNV,colnames(Gene)[geneWithout])) == 0 ){
    C1 = rep(0,nrow(y))
    C2 = rep(0,nrow(y))
  }else{
    if (intra){
      for (g in geneWithCNV){
        sub_CNV = as.character(Gene_CNV[Gene_CNV[,2] == g,1] )
        ind_CNV =  which(colnames(CNV)%in%sub_CNV)
        if (length(ind_CNV)!=0){ 
          part_CNV = as.matrix(CNV[,ind_CNV])
          colnames(part_CNV) = colnames(CNV)[ind_CNV]
          X = part_CNV
          y_tmp = as.matrix(Gene[,g])
          colnames(y_tmp) = g
          resC[[g]] = regression(x=X,y=y_tmp,alpha =1 ,nfolds= nfolds)
          select = resC[[g]]$active_index 
          ##gives us all variables selected to explain the gene expression
          ind = which(colnames(part_CNV)%in%names(select))
          ind_tmp =  which(names(select)%in%colnames(part_CNV))
          if ((length(ind)!=1)&(length(select)!=0))
          	{lol= part_CNV[,ind]%*%(resC[[g]]$active[ind_tmp])}
          if (length(ind)==1){lol= part_CNV[,ind]*resC[[g]]$active[ind_tmp]}
          if (length(ind)==0){lol = as.matrix(methy[,1]*0)}
          RNA_modulated_CNV[[g]] = lol             
        } else {
          RNA_modulated_CNV[[g]] = as.matrix(methy[,1]*0)
          resC[[g]]$inter = lm(Gene[,g]~1)$coef[1]
        }
        RNA_modulated_other[[g]]= Gene[,g]-RNA_modulated_CNV[[g]]-resC[[g]]$inter
      }
    }else{
      for(g in c(geneWithCNV,colnames(Gene)[geneWithout])){        
        if (pathway){
          path = Gene_Pathway[which(Gene_Pathway[,1]%in%g),2]
          gene_in_path = as.character(Gene_Pathway[which(Gene_Pathway[,2]%in%path),1])
          sub_CNV = as.character(Gene_CNV[which(Gene_CNV[,2]%in%gene_in_path),1])
        }else{      
          sub_CNV = colnames(CNV)
        }
        ind_CNV = which(colnames(CNV)%in%sub_CNV)
        if (length(ind_CNV)!=0){ 
          part_CNV = as.matrix(CNV[,ind_CNV])
          colnames(part_CNV) = colnames(CNV)[ind_CNV]
          X = part_CNV
          y_tmp = as.matrix(Gene[,g])
          colnames(y_tmp) = g
          resC[[g]] = regression(x=X,y=y_tmp,alpha =1 ,nfolds= nfolds)
          select = resC[[g]]$active_index 
          ##gives us all variables selected to explain the gene expression
          ind = which(colnames(part_CNV)%in%names(select))
          ind_tmp =  which(names(select)%in%colnames(part_CNV))
          if ((length(ind)!=1)&(length(select)!=0))
          	{lol= part_CNV[,ind]%*%(resC[[g]]$active[ind_tmp])}
          if (length(ind)==1){lol= part_CNV[,ind]*resC[[g]]$active[ind_tmp]}
          if (length(ind)==0){lol = as.matrix(methy[,1]*0)}
          RNA_modulated_CNV[[g]] = lol             
        }else{
          RNA_modulated_CNV[[g]] = as.matrix(methy[,1]*0)
          resC[[g]]$inter = lm(Gene[,g]~1)$coef[1]
        }  
        RNA_modulated_other[[g]]= Gene[,g]-RNA_modulated_CNV[[g]]-resC[[g]]$inter
      }
    }
    C1 = matrix(unlist(RNA_modulated_CNV),ncol=length(RNA_modulated_CNV),byrow=FALSE)
    colnames(C1) = paste(names(resC),"CNV",sep="_")
 
    C2 =  matrix(unlist(RNA_modulated_other),ncol=length(RNA_modulated_other),byrow=FALSE)
    colnames(C2) = paste(names(resC),"other",sep="_")
  }
 
#########################################################
##########   4. Final analysis     ##############################
#########################################################

  indCNV_RNA = unique(unlist(lapply(res3,function(y){names(y$active_index)})))
  indCNV_RNAbis = unique(unlist(lapply(resC,function(y){names(y$active_index)}),
  	lapply(resB,function(y){names(y$active_index)})))
  indCNV_y = sort(unique(c(indCNV_RNA,indCNV_RNAbis)))
  if(length(indCNV_y)!=0){ CNVy = CNV[,-which(colnames(CNV)%in%indCNV_y)]
  }else{
    CNVy =rep(0,n)
  }
  

  R2 = rep(0,4)
  R2.adj = rep(0,4)
  cIndex = rep(0,4)
  res = coxph(Surv(y$y,y$Event)~ 1+y$Age)
  cIndex[1] = summary(res)$concordance[1]
  logNull = res$log[1]
  logAge = res$log[2]
  R2[1] = 1-exp((-2/n)*(logAge-logNull))
  R2.adj[1] = 1-(1-R2[1])*(n-1)/(n-2)
  if (intra){X = cbind(y$Age,A1,A2,A3,B1,B2,B3,C1,C2,Gene[,geneWithout], CNVy)
               Xmeth1 = X[,which(apply(X,2,sum)!=0)]
               colnames(Xmeth1)[1] ="Age"
               
               X.2 = cbind(y$Age,A1,A2,A3,B1,B2,B3,C1,C2,Gene[,geneWithout])
               X.meth.2 = X.2[,which(apply(X.2,2,sum)!=0)]
               colnames(X.meth.2)[1] ="Age"
    }else{
      X = cbind(y$Age,A1,A2,A3,B1,B2,B3,C1,C2, CNVy)
      Xmeth1 = X[,which(apply(X,2,sum)!=0)]
      colnames(Xmeth1)[1] ="Age"
      
      X.2 = cbind(y$Age,A1,A2,A3,B1,B2,B3,C1,C2)
      X.meth.2 = X.2[,which(apply(X.2,2,sum)!=0)]
      colnames(X.meth.2)[1] ="Age"
    }
    
    cv.fit <- cv.glmnet(x=Xmeth1,
                        Surv(time=y$y, event=y$Event,type="right"), 
                        family = "cox",alpha=alpha,maxit=2000,
                        penalty.factor=c(1,rep(1,ncol(Xmeth1))),
                        nfolds=nfolds)
    
    fit <- glmnet(x=Xmeth1,
                  Surv(time=y$y, event=y$Event,type="right"),
                  family = "cox",alpha=alpha,
                  penalty.factor=c(1,rep(1,ncol(Xmeth1)-1)))
    
    coef <- coef(fit, s = cv.fit$lambda.min)
    ind_coef = which(colnames(Xmeth1)%in%rownames(coef)[which(coef!=0)])
    cox = coxph(Surv(time=y$y, event=y$Event,type="right")~ Xmeth1[,ind_coef])
    cIndex[3] = summary(cox)$concordance[1]
    logMod1 = cox$loglik[2]
    R2[3] = 1-exp((-2/n)*(logMod1-logNull))
    p = as.numeric(length(coef[coef!=0]))
    R2.adj[3] = 1-(1-R2[3])*(n-1)/(n-p-1)


    cv.fit.2 <- cv.glmnet(x=X.meth.2,
                    Surv(time=y$y, event=y$Event,type="right"), 
                    family = "cox",alpha=alpha,maxit=2000,
                    penalty.factor=c(1,rep(1,ncol(Xmeth1))),
                    nfolds=nfolds)

    fit.2 <- glmnet(x=X.meth.2,
              Surv(time=y$y, event=y$Event,type="right"),
              family = "cox",alpha=alpha,
              penalty.factor=c(1,rep(1,ncol(Xmeth1)-1)))

    coef.2 <- coef(fit.2, s = cv.fit.2$lambda.min)
    ind_coef.2 = which(colnames(X.meth.2)%in%rownames(coef.2)[which(coef.2!=0)])
    cox = coxph(Surv(time=y$y, event=y$Event,type="right")~ X.meth.2[,ind_coef.2])
    cIndex[4] = summary(cox)$concordance[1]
    logMod1 = cox$loglik[2]
    R2[4] = 1-exp((-2/n)*(logMod1-logNull))
    p = as.numeric(length(coef.2[coef.2!=0]))
    R2.adj[4] = 1-(1-R2[4])*(n-1)/(n-p-1)


    X =cbind(y$Age,Gene)
    cv.fit4 <- cv.glmnet(x=X,
                         Surv(time=y$y, event=y$Event,type="right"), 
                         family = "cox",alpha=alpha,maxit=1000,
                         penalty.factor=c(1,rep(1,ncol(X)-1)),
                         nfolds=nfolds)
    
    fit4 <- glmnet(x=X, Surv(time=y$y, event=y$Event,type="right"), 
                   family = "cox",alpha=alpha,
                   penalty.factor=c(1,rep(1,ncol(X)-1)))
    
    coef4 <- coef(fit4, s = cv.fit4$lambda.min)
    ind_coef4 = which(colnames(X)%in%rownames(coef4)[which(coef4!=0)])
    cox4 = coxph(Surv(time=y$y, event=y$Event,type="right")~ X[,ind_coef4])
    cIndex[2] = summary(cox4)$concordance[1]
    logMod4 = cox4$loglik[2]
    R2[2]= 1-exp((-2/n)*(logMod4-logNull))
    p = as.numeric(length(coef4[coef4!=0]))
    R2.adj[2] = 1-(1-R2[2])*(n-1)/(n-p-1)
  
  
  return(list(R2=R2, R2.adj=R2.adj, cIndex=cIndex, Xmeth1=Xmeth1, X=X,
  	coef=coef, coef4=coef4, coefNo=coef.2))
}
