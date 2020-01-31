
library(narray)
library(biclust)
library(reshape)
library(RColorBrewer)
library(ggplot2)

#############################
######### TEMPLATE ##########
#############################

#' medoid_evaluation evaluates the medoid template function
#' @noRd
medoid_evaluation <- function(fun_mat,a,b,const_a,const_b){
  n=dim(fun_mat)[1]# numero di righe
  m=dim(fun_mat)[2] # numero di colonne
  p=dim(fun_mat)[3] # lunghezza griglia temporale

  fun_per_medoid=NULL
  for(j in 1:dim(fun_mat)[1]){
    fun_per_medoid=rbind(fun_per_medoid,fun_mat[j,,])
  }
  # calcolo le distanze euclidee tra tutte le funzioni
  #dist <- distNumeric(fun_per_medoid, fun_per_medoid, method = "se")
  distance <- as.matrix(stats::dist(fun_per_medoid))
  distance <- distance^2
  diag(distance) <- NA
  sum_dist <- colMeans(distance,na.rm=T)
  # trovo la funzione medoide con somma delle distanze minima
  rep_curve <- which(sum_dist==min(sum_dist,na.rm=T))[1]
  medoid_fun <- fun_per_medoid[rep_curve,]

  # costruisco l'array da fun_medoid n x m x p
  new_fun = array(medoid_fun, dim=c(1,1,p))
  new_fun = narray::rep(new_fun, n, along = 1)
  new_fun = narray::rep(new_fun, m, along = 2)

  new_fun
}

#' medoid_evaluation_add evaluates the medoid template function
#' @noRd
#' @keywords Internal
medoid_evaluation_add <- function(fun_mat,logr,logc,a,b,const_a,const_b){
  n=dim(fun_mat)[1]# numero di righe
  m=dim(fun_mat)[2] # numero di colonne
  p=dim(fun_mat)[3] # lunghezza griglia temporale

  fun_mat <- fun_mat[logr,logc,]

  fun_per_medoid=NULL
  for(j in 1:dim(fun_mat)[1]){
    fun_per_medoid=rbind(fun_per_medoid,fun_mat[j,,])
  }

  # calcolo le distanze euclidee tra tutte le funzioni
  distance <- as.matrix(stats::dist(fun_per_medoid))
  distance <- distance^2
  diag(distance) <- NA


  sum_dist <- colMeans(distance,na.rm=T)
  # trovo la funzione medoide con somma delle distanze minima
  rep_curve <- which(sum_dist==min(sum_dist,na.rm=T))[1]
  medoid_fun <- fun_per_medoid[rep_curve,]

  # costruisco l'array da fun_medoid n x m x p
  new_fun = array(medoid_fun, dim=c(1,1,p))
  new_fun = narray::rep(new_fun, n, along = 1)
  new_fun = narray::rep(new_fun, m, along = 2)

  new_fun
}

#' template_evaluation evaluates the template function
#' @noRd
#' @keywords Internal
template_evaluation <- function(fun_mat,a,b,const_a,const_b){
  n=dim(fun_mat)[1]# numero di righe
  m=dim(fun_mat)[2] # numero di colonne
  p=dim(fun_mat)[3] # lunghezza griglia temporale

  count_null <- apply(fun_mat, c(1,2), function(x) sum(is.na(x)))
  not_null <- count_null < dim(fun_mat)[3]

  #print(dim(fun_mat))
  fun_mean=colMeans(fun_mat, dims = 2,na.rm=T)  # 1 x p # calcolo la funzione media della fun_matrice

  #calcolo alpha - componente riga
  alpha_fun=NULL

  for(j in 1:n){
    alpha_fun=rbind(alpha_fun,colSums(fun_mat[j,,],na.rm=T)/sum(not_null[j,]))
  }
  #dim(alpha_fun) #  n x p
  alpha_fun=alpha_fun-matrix(fun_mean,nrow=n,ncol=length(fun_mean),byrow=TRUE)

  if(const_a){
    alpha_fun = rowMeans(alpha_fun,na.rm=T) # 1 x n
    alpha_fun = array(alpha_fun,dim=c(n,1,1)) # n x 1 x 1
    alpha_fun = narray::rep(alpha_fun, p, along = 3) # n x p
  }

  # calcolo beta - componente colonna
  beta_fun=NULL
  for(j in 1:m){
    beta_fun=rbind(beta_fun,(colSums(fun_mat[,j,],dims = 1,na.rm=T)/sum(not_null[,j]))) # m x p
  }
  beta_fun=beta_fun-matrix(fun_mean,nrow=m,ncol=length(fun_mean),byrow=TRUE) # m x p

  if(const_b){
    beta_fun=rowMeans(beta_fun,na.rm=T) # 1 x m
    beta_fun = array(beta_fun,dim=c(1,m,1)) # 1 x m x 1
    beta_fun = narray::rep(beta_fun, p, along = 3) # m x p
  }

  # costruisco l'array da fun_mean n x m x p
  fun_mean_mat = array(fun_mean, dim=c(1,1,p))
  fun_mean_mat = narray::rep(fun_mean_mat, n, along = 1)
  fun_mean_mat = narray::rep(fun_mean_mat, m, along = 2)

  # costruisco l'array da beta_fun n x m x p
  beta_fun_mat = array(beta_fun,dim=c(1,m,p))
  beta_fun_mat = narray::rep(beta_fun_mat, n, along = 1)
  beta_fun_mat[is.na(beta_fun_mat)] <- 0

  # costruisco l'array da alpha_fun n x m x p
  alpha_fun_mat = array(alpha_fun,dim=c(n,1,p))
  alpha_fun_mat = narray::rep(alpha_fun_mat, m, along = 2)
  alpha_fun_mat[is.na(alpha_fun_mat)] <- 0

  new_fun = fun_mean_mat+b*beta_fun_mat+a*alpha_fun_mat
  new_fun
}

#' template_evaluation_add evaluates the template function
#' @noRd
#' @keywords Internal
template_evaluation_add <- function(fun_mat,logr,logc,a,b,const_a,const_b){
  n=dim(fun_mat)[1]# numero di righe
  m=dim(fun_mat)[2] # numero di colonne
  p=dim(fun_mat)[3] # lunghezza griglia temporale

  count_null <- apply(fun_mat, c(1,2), function(x) sum(is.na(x)))
  not_null <- count_null < dim(fun_mat)[3]

  fun_mean=colMeans(fun_mat[logr,logc,], dims = 2,na.rm=T)  # 1 x p # calcolo la funzione media della fun_matrice

  #calcolo alpha - componente riga
  alpha_fun=NULL
  for(j in 1:n){
    alpha_fun=rbind(alpha_fun,colSums(fun_mat[j,logc,],na.rm=T)/sum(not_null[j,logc],na.rm=T))
  }
  #dim(alpha_fun) #  n x p
  alpha_fun=alpha_fun-matrix(fun_mean,nrow=n,ncol=length(fun_mean),byrow=TRUE)

  if(const_a){
    alpha_fun = rowMeans(alpha_fun,na.rm=T) # 1 x n
    alpha_fun = array(alpha_fun,dim=c(n,1,1)) # n x 1 x 1
    alpha_fun = narray::rep(alpha_fun, p, along = 3) # n x p

  }


  # calcolo beta - componente colonna
  beta_fun=NULL
  for(j in 1:m){
    beta_fun=rbind(beta_fun,(colSums(fun_mat[logr,j,],na.rm=T)/sum(not_null[logr,j],na.rm=T))) # m x p
  }
  beta_fun=beta_fun-matrix(fun_mean,nrow=m,ncol=length(fun_mean),byrow=TRUE) # m x p


  if(const_b){
    beta_fun=rowMeans(beta_fun,na.rm=T) # 1 x m
    beta_fun = array(beta_fun,dim=c(1,m,1)) # 1 x m x 1
    beta_fun = narray::rep(beta_fun, p, along = 3) # m x p
  }


  # costruisco l'array da fun_mean n x m x p
  fun_mean_mat = array(fun_mean, dim=c(1,1,p))
  fun_mean_mat = narray::rep(fun_mean_mat, n, along = 1)
  fun_mean_mat = narray::rep(fun_mean_mat, m, along = 2)

  # costruisco l'array da beta_fun n x m x p
  beta_fun_mat = array(beta_fun,dim=c(1,m,p))
  beta_fun_mat = narray::rep(beta_fun_mat, n, along = 1)
  beta_fun_mat[is.na(beta_fun_mat)] <- 0

  # costruisco l'array da alpha_fun n x m x p
  alpha_fun_mat = array(alpha_fun,dim=c(n,1,p))
  alpha_fun_mat = narray::rep(alpha_fun_mat, m, along = 2)
  alpha_fun_mat[is.na(alpha_fun_mat)] <- 0

  new_fun = fun_mean_mat+b*beta_fun_mat+a*alpha_fun_mat
  new_fun
}

#############################
######### SCORE ##########
#############################

#' evaluate_mat_dist evaluates distance matrix
#' @noRd
#' @keywords Internal
evaluate_mat_dist <- function(fun_mat,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter){
  n=dim(fun_mat)[1]# numero di righe
  m=dim(fun_mat)[2] # numero di colonne
  p=dim(fun_mat)[3] # lunghezza griglia temporale

  if(shift.alignement!=F){
    mat_dist=warping_function(fun_mat,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)
  }

  else{
    if(template.type=='mean'){    new_fun <- template_evaluation(fun_mat,a,b,const_a,const_b)}
    if(template.type=='medoid'){    new_fun <- medoid_evaluation(fun_mat,a,b,const_a,const_b)}

    mat_dist <- rowSums((fun_mat-new_fun)^2,dims = 2)/p
  }
  mat_dist
}

#' evaluate_mat_dist_add evaluates distance matrix
#' @noRd
#' @keywords Internal
evaluate_mat_dist_add <- function(fun_mat,logr,logc,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter){
  n=dim(fun_mat)[1]# numero di righe
  m=dim(fun_mat)[2] # numero di colonne
  p=dim(fun_mat)[3] # lunghezza griglia temporale

  if(shift.alignement!=F){
    mat_dist=warping_function_add(fun_mat,logr,logc,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)
  }

  else{
    if(template.type=='mean'){new_fun <- template_evaluation_add(fun_mat,logr,logc,a,b,const_a,const_b)}
    if(template.type=='medoid'){new_fun <- medoid_evaluation_add(fun_mat,logr,logc,a,b,const_a,const_b)}

    mat_dist <- rowSums((fun_mat-new_fun)^2,dims = 2)/p
  }
  mat_dist
}

#' ccscore_fun find submatrices score
#' @noRd
#' @keywords Internal
ccscore_fun<-function(mat_dist){
  #n=dim(mat_dist)[1]# numero di righe
  #m=dim(mat_dist)[2] # numero di colonne

  score_fun = mean(mat_dist,na.rm=T)

  score_fun
}

#' rowscore_fun find row score
#' @noRd
#' @keywords Internal
rowscore_fun<-function(mat_dist){
  #m=dim(mat_dist)[2] # numero di colonne

  score_fun = rowMeans(mat_dist,na.rm=T) #rowSums(mat_dist,na.rm=T)/(m)

  score_fun

}

#' colscore_fun find col score
#' @noRd
#' @keywords Internal
colscore_fun<-function(mat_dist){
  #n=dim(mat_dist)[1]# numero di righe

  score_fun = colMeans(mat_dist,na.rm=T)#colSums(mat_dist,na.rm=T)/(n)

  score_fun

}

#' addrowscore_fun find row score when adding a row
#' @noRd
#' @keywords Internal
addrowscore_fun<-function(mat_dist){

  #m=dim(mat_dist)[2] # numero di colonne

  score_fun = rowMeans(mat_dist,na.rm=T)#rowSums(mat_dist,na.rm=T)/(m)

  score_fun

}

#' addcolscore_fun find col score when adding a col
#' @noRd
#' @keywords Internal
addcolscore_fun<-function(mat_dist){

  score_fun = colMeans(mat_dist,na.rm=T)#colSums(mat_dist,na.rm=T)/(n)

  score_fun
  # vettore lungo tanto quanto le colonne (giorni)
  # mi serve per capire che colonna togliere
}


#############################
######### WARPING ###########
#############################
#' warping_function alignment function
#' @noRd
#' @keywords Internal
warping_function<- function(fun_mat,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter){

  # fun_mat_align <- fun_mat
  # new_fun <- template_evaluation(fun_mat_align,a,b,const_a,const_b)
  #
  warping.shift <- function(coeff) {
    st <- x.reg + coeff

    template.t <- new_fun[i,j,]
    b <- !is.na(template.t)

    data.t <- stats::approx(st, fun_mat[i,j,], xout = x.out)$y
    a <- !is.na(data.t)

    sel <- a & b
    data.t <- data.t[sel]
    template.t <- new_fun[i,j,sel]

    distance = sum((data.t-template.t)^2)/sum(sel)
    ##print(distance)

    distance

  }


  ## definisco lower and upper warp
  n = dim(fun_mat)[1]
  m = dim(fun_mat)[2]
  p = dim(fun_mat)[3]

  # pensare se p iniziale o del dominio nuovo
  x <- seq(1,p)
  x.reg <- x
  x.out <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length = p)

  min.temp <- diff(range(x,na.rm=T))
  lower.warp <- -shift.max*min.temp
  upper.warp <- shift.max*min.temp

  dist_mat = matrix(10,nrow=n,ncol=m)
  align_mat = matrix(0,nrow=n,ncol=m)

  count_null <- apply(fun_mat, c(1,2), function(x) sum(is.na(x)))
  not_null <- count_null < dim(fun_mat)[3]

  fun_mat_align <- fun_mat
  iter <- 0
  conv=F
  # iterazioni di allineamento fino a iter.max o convergenza
  while(iter<max.iter & conv==FALSE){
    iter <- iter+1
    ##print(iter)
    dist_mat_old <- dist_mat

    # calcolo i nuovi template
    if(template.type=='mean'){new_fun <- template_evaluation(fun_mat_align,a,b,const_a,const_b)}
    if(template.type=='medoid'){new_fun <- medoid_evaluation(fun_mat_align,a,b,const_a,const_b)}


    if(shift.alignement==T){

        for(i in 1:n){
          for(j in 1:m){
            if(not_null[i,j]){
              result = stats::optim(c(0),warping.shift,method='Brent',lower=lower.warp,upper=upper.warp)
              #result = optim(c(0),warping.shift,method=optim.method)#,lower=lower.warp[2],upper=upper.warp[2])
              dist_mat[i,j] = result$value
              align_mat[i,j] =  result$par
              new_x = seq(1,p) + align_mat[i,j]
              fun_mat_align[i,j,] <- stats::approx(new_x, fun_mat[i,j,], xout = x.out)$y
            }

            if(!not_null[i,j]){
              dist_mat[i,j] = NA
              align_mat[i,j] =  NA
              fun_mat_align[i,j,] <- rep(NA,dim(fun_mat)[3])
            }

        }
      }
    }

    ##print(paste0('dist_mat_old:',sum(dist_mat_old),' - ','dist_mat:',sum(dist_mat)))
    if(sum(dist_mat_old,na.rm=T)<=sum(dist_mat,na.rm=T)){conv=TRUE; iter=iter-1}


  }

  dist_mat_old
}

#' warping_function_add alignment function
#' @noRd
#' @keywords Internal
warping_function_add <- function(fun_mat,logr,logc,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter){

  # fun_mat_align <- fun_mat
  # new_fun <- template_evaluation(fun_mat_align,a,b,const_a,const_b)


  #
  warping.shift <- function(coeff) {
    st <- x.reg + coeff

    template.t <- new_fun[i,j,]
    b <- !is.na(template.t)

    data.t <- stats::approx(st, fun_mat[i,j,], xout = x.out)$y
    a <- !is.na(data.t)

    sel <- a & b
    data.t <- data.t[sel]
    template.t <- new_fun[i,j,sel]

    distance = sum((data.t-template.t)^2)/sum(sel)
    #print(distance)

    distance

  }


  ## definisco lower and upper warp
  n = dim(fun_mat)[1]
  m = dim(fun_mat)[2]
  p = dim(fun_mat)[3]

  # pensare se p iniziale o del dominio nuovo
  x <- seq(1,p)
  x.reg <- x
  x.out <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length = p)

  min.temp <- diff(range(x,na.rm=T))
  lower.warp <- -shift.max*min.temp
  upper.warp <- shift.max*min.temp

  dist_mat = matrix(10,nrow=n,ncol=m)
  align_mat = matrix(0,nrow=n,ncol=m)

  count_null <- apply(fun_mat, c(1,2), function(x) sum(is.na(x)))
  not_null <- count_null < dim(fun_mat)[3]

  fun_mat_align <- fun_mat
  iter <- 0
  conv=F
  # iterazioni di allineamento fino a iter.max o convergenza
  while(iter<max.iter & conv==FALSE){
    iter <- iter+1
    #print(iter)
    dist_mat_old <- dist_mat

    # calcolo i nuovi template
    if(template.type=='mean'){new_fun <- template_evaluation_add(fun_mat_align,logr,logc,a,b,const_a,const_b)}
    if(template.type=='medoid'){new_fun <- medoid_evaluation_add(fun_mat_align,logr,logc,a,b,const_a,const_b)}

    if(shift.alignement==T){

        for(i in 1:n){
          for(j in 1:m){
            if(not_null[i,j]){
              result = stats::optim(c(0),warping.shift,method='Brent',lower=lower.warp,upper=upper.warp)
              #result = optim(c(0),warping.shift,method=optim.method)#,lower=lower.warp[2],upper=upper.warp[2])
              #if(is.na(result$value)){print('!!!!!!!!!!')}
              dist_mat[i,j] = result$value
              align_mat[i,j] =  result$par
              new_x = seq(1,p) + align_mat[i,j]
              fun_mat_align[i,j,] <- stats::approx(new_x, fun_mat[i,j,], xout = x.out)$y}

            if(!not_null[i,j]){
              dist_mat[i,j] = NA
              align_mat[i,j] =  NA
              fun_mat_align[i,j,] <- rep(NA,dim(fun_mat)[3])
            }
          }

      }
    }

    #print(paste0('dist_mat_old:',sum(dist_mat_old),' - ','dist_mat:',sum(dist_mat)))
    if(sum(dist_mat_old,na.rm=T)<=sum(dist_mat,na.rm=T)){conv=TRUE; iter=iter-1}


  }

  dist_mat_old
}

########################################
######### addition and deletion ###########
########################################

#' cc1_fun algorithm 1 from CC: Single Node Deletion
#' @noRd
#' @keywords Internal
cc1_fun<-function(fun_mat,logr,logc,delta,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter){

  dist_mat <- evaluate_mat_dist(fun_mat[logr,logc,],template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)
  score_while <- ccscore_fun(dist_mat)

  #i <- 1
  while(score_while>delta)
  {

    di<-rowscore_fun(dist_mat)
    dj<-colscore_fun(dist_mat)
    mdi<-which.max(di)
    mdj<-which.max(dj)

    ifelse(di[mdi]>dj[mdj] ,logr[logr][mdi]<-FALSE ,logc[logc][mdj]<-FALSE)

    ##print(logr)
    if (!(sum(logr)>1 & sum(logc)>1))
      break

    #print(sum(logr))
    #print(sum(logc))
    dist_mat <- evaluate_mat_dist(fun_mat[logr,logc,],template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)
    #print(dist_mat)
    score_while <- ccscore_fun(dist_mat)
    #print(score_while)


  }
  ifelse(sum(logr)>1 & sum(logc)>1,ret<-list(logr,logc),ret<-list(0,warning(paste('No fun_matirx with score smaller', delta,'found'))))
  ret
}

#' cc2_fun algorithm 2 from CC: Multiple Node Deletion
#' @noRd
#' @keywords Internal
cc2_fun<-function(fun_mat,logr,logc,delta,theta,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter){
  mdi<-1
  mdj<-1
  dist_mat <- evaluate_mat_dist(fun_mat[logr,logc,],template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)
  h<-ccscore_fun(dist_mat)

  while(h>delta & (sum(mdi,na.rm=T)+sum(mdj,na.rm=T))>0)
  {

    if(sum(logr)>100)
    {
      di<-rowscore_fun(dist_mat)
      mdi<-di>(theta*h)
      if(sum(mdi,na.rm=T) < (sum(logr[!is.na(mdi)],na.rm=T)-1))
      {
        logr[logr][mdi]<-FALSE
        dist_mat <- evaluate_mat_dist(fun_mat[logr,logc,],template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)
        h<-ccscore_fun(dist_mat)
      }
      else
      {
        print(warning(paste('Theta', theta,'to small!')))
        mdi <- 0
      }
    }
    else{mdi<-0}


    if(sum(logc)>100)
    {
      dj<-colscore_fun(dist_mat)
      mdj<-dj>(theta*h)
      if(sum(mdj,na.rm=T) < (sum(logc[!is.na(mdj)])-1))
      {
        logc[logc][mdj]<-FALSE
        dist_mat <- evaluate_mat_dist(fun_mat[logr,logc,],template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)
      }
      else
      {

        print(warning(paste('theta', theta,'to small!')))
        mdi <- 0
      }
    }
    else{mdj<-0}

    h<-ccscore_fun(dist_mat)
  }

  ret<-list(logr,logc)
  ret
}

#' cc3_fun algorithm 3 from CC:  Node Addition
#' @noRd
#' @keywords Internal
cc3_fun<-function(fun_mat,logr,logc,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter){
  br<-1
  ilogr<-rep(FALSE,length(logr))
  while(br>0)
  {

    # dist mat
    dist_mat <- evaluate_mat_dist(fun_mat[logr,logc,],template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)

    br1<-sum(logc)
    br2<-sum(logr)
    h<-ccscore_fun(dist_mat)
    #print(paste0('h: ',h))

    dist_mat_add <- evaluate_mat_dist_add(fun_mat,logr,logc,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)

    dj<-addcolscore_fun(dist_mat_add)
    #print(dj)

    mdj<-dj<=h
    logc[mdj]<-TRUE

    dist_mat <- evaluate_mat_dist(fun_mat[logr,logc,],template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)

    h<-ccscore_fun(dist_mat)

    dist_mat_add <- evaluate_mat_dist_add(fun_mat,logr,logc,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)

    di<-addrowscore_fun(dist_mat_add)

    mdi<-di<=h
    logr[mdi]<-TRUE

    br<-sum(logc)+sum(logr)-br1-br2
  }
  ret<-list(logr,logc)
  ret

}

########################################
######### Find biggest Bicluster: ###########
########################################
#' bigcc_fun Finds biggest Bicluster
#' @noRd
#' @keywords Internal
bigcc_fun<-function(fun_mat,delta,theta,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter){

  n=dim(fun_mat)[1]
  m=dim(fun_mat)[2]
  p=dim(fun_mat)[3]

  logr<-rep(TRUE,n)
  logr[rowSums(is.na(fun_mat[,,1]))==dim(fun_mat)[2]] <- FALSE
  logc<-rep(TRUE,m)
  logc[colSums(is.na(fun_mat[,,1]))==dim(fun_mat)[1]] <- FALSE

  # if(sum(logr)<=1 | sum(logc)<=1){
  #   ret<-list(0,warning(paste('Mo fun_matrix with score smaller than', delta,'found')))
  #   ret
  #   break
  # }

  #print('cc2')
  step1<-cc2_fun(fun_mat,logr,logc,delta,theta,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)
  #print('cc1')
  step2<-cc1_fun(fun_mat,step1[[1]],step1[[2]],delta,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)
  if(sum(step2[[1]])==0)
  {ret<-list(0,warning(paste('Mo fun_matrix with score smaller than', delta,'found')))
  }
  else{
    #print('cc3')

    ret<-cc3_fun(fun_mat,step2[[1]],step2[[2]],template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter)
  }
  ret
}

#' rep.row
#' @noRd
#' @keywords Internal
rep.row<-function(x,n){
  matrix(base::rep(x,each=n),nrow=n)
}
#' rep.col
#' @noRd
#' @keywords Internal
rep.col<-function(x,n){
  matrix(base::rep(x,each=n), ncol=n, byrow=TRUE)
}

#' warping_function_plot alignment function for figures
#' @noRd
#' @keywords Internal
warping_function_plot<- function(res,fun_mat,template.type,a,b,const_a,const_b,shift.alignement,shift.max, max.iter){

  # estrarre parametri da res

  # fun_mat_align <- fun_mat
  # new_fun <- template_evaluation(fun_mat_align,a,b,const_a,const_b)
  #
  warping.shift <- function(coeff) {
    st <- x.reg + coeff

    template.t <- new_fun[i,j,]
    b <- !is.na(template.t)

    data.t <- stats::approx(st, fun_mat[i,j,], xout = x.out)$y
    a <- !is.na(data.t)

    sel <- a & b
    data.t <- data.t[sel]
    template.t <- new_fun[i,j,sel]

    distance = sum((data.t-template.t)^2)/sum(sel)
    #print(distance)

    distance

  }


  ## definisco lower and upper warp
  n = dim(fun_mat)[1]
  m = dim(fun_mat)[2]
  p = dim(fun_mat)[3]

  # pensare se p iniziale o del dominio nuovo
  x <- seq(1,p)
  x.reg <- x
  x.out <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length = p)

  min.temp <- diff(range(x,na.rm=T))
  lower.warp <- -shift.max*min.temp
  upper.warp <- shift.max*min.temp

  dist_mat = matrix(10,nrow=n,ncol=m)
  align_mat = matrix(0,nrow=n,ncol=m)

  count_null <- apply(fun_mat, c(1,2), function(x) sum(is.na(x)))
  not_null <- count_null < dim(fun_mat)[3]


  fun_mat_align <- fun_mat
  iter <- 0
  conv=F
  # iterazioni di allineamento fino a iter.max o convergenza
  while(iter<max.iter & conv==FALSE){
    iter <- iter+1
    #print(iter)
    ##print(iter)
    dist_mat_old <- dist_mat
    fun_mat_align_old <- fun_mat_align
    # calcolo i nuovi template
    if(template.type=='mean'){new_fun <- template_evaluation(fun_mat_align,a,b,const_a,const_b)}
    if(template.type=='medoid'){new_fun <- medoid_evaluation(fun_mat_align,a,b,const_a,const_b)}

    if(shift.alignement==T){

        for(i in 1:n){
          for(j in 1:m){
            if(not_null[i,j]){
              result = stats::optim(c(0),warping.shift,method='Brent',lower=lower.warp,upper=upper.warp)
              #result = optim(c(0),warping.shift,method=optim.method)#,lower=lower.warp[2],upper=upper.warp[2])
              dist_mat[i,j] = result$value
              align_mat[i,j] =  result$par
              new_x = seq(1,p) + align_mat[i,j]
              fun_mat_align[i,j,] <- stats::approx(new_x, fun_mat[i,j,], xout = x.out)$y
            }

            if(!not_null[i,j]){
              dist_mat[i,j] = NA
              align_mat[i,j] =  NA
              fun_mat_align[i,j,] <- rep(NA,dim(fun_mat)[3])
            }


        }
      }
    }

    #print(paste0('dist_mat_old:',sum(dist_mat_old),' - ','dist_mat:',sum(dist_mat)))
    if(sum(dist_mat_old,na.rm=T)<=sum(dist_mat,na.rm=T)){conv=TRUE; iter=iter-1; fun_mat_align=fun_mat_align_old}


  }

  if(template.type=='mean'){new_fun <- template_evaluation(fun_mat_align,a,b,const_a,const_b)}
  if(template.type=='medoid'){new_fun <- medoid_evaluation(fun_mat_align,a,b,const_a,const_b)}

  coeff_mat=align_mat

  x.out=matrix(x,nrow=n*m,ncol=p,byrow = T)
  x.align=c(matrix(coeff_mat,nrow=n*m,1,byrow =T))
  x.out=x.out+x.align

  res=list(fun_mat_align=fun_mat_align,template=new_fun,x.out=x.out)
  return(res)
}
