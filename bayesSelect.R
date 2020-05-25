rm(list = ls())
library("rjags")
library("TruncatedNormal")



#######################################################
#######################################################

inistring <- function(varstr, mixed = F, indMea = NA,  par.a=1, par.b=1){
    
    betastr <- paste0( "b.", varstr )
    Xstr <- paste0("X.", varstr, "[j]")
    
    parstr <- paste0( betastr, "*", Xstr, collapse = " + ")
    
    if( is.na(varstr[1]) ){
      
      ##############################
      ######Non-mixed Model######
      ##############################
      if( mixed == F ){
        priorstr <- paste0( "beta0 ~ dnorm(0, 0001) ", " phi ~ dgamma(",par.a, ", ", par.b,")" )
        
        modelstring <- paste0(
          "Y[j] ~ dbeta(alpha[j], beta[j])
          alpha[j] <- mu[j] * phi
          beta[j] <- (1 - mu[j]) * phi
          logit( mu[j] ) <- beta0  
          r[j]=Y[j] - mu[j]
      }
          ")
        
        finalstring <- paste0(
          "model{
          for(j in 1:N ){
          ", modelstring, priorstr, "}")
        
    }
      #############################
      ######Mixed Model#############
      ################################
      else{
        print("run mixed models")
        
        if( is.na(indMea) == T ){
          priorstr <- paste0( "mu.beta0 ~ dnorm(0, 0.001) tau.beta0 ~ dgamma(1,1) ", " phi ~ dgamma(1, 1)" )
          mixpriorstr <- "beta0[i] ~ dnorm(mu.beta0, tau.beta0) "
          
          modelstring <- paste0(
            "Y[j] ~ dbeta(alpha[j], beta[j])
            alpha[j] <- mu[j] * phi
            beta[j] <- (1 - mu[j]) * phi
            logit( mu[j] ) <- beta0[i] 
            r[j]=Y[j] - mu[j]
        }
            " )
          
          finalstring <- paste0(
            "model{
            for(i in 1:(Nset-1) ){
            for(j in offsets[i]:(offsets[i+1]-1)){
            ",modelstring, mixpriorstr,"
            }
            ", priorstr, "}")
      }
        else{
          priorstr <- paste0( "mu.beta0 ~ dnorm(0, 0.001) ", " phi ~ dgamma(1, 1)" )
          mixpriorstr <- "beta0[i] ~ dnorm(mu.beta0, tau.beta0[i]) tau.beta0[i] ~ dgamma(1,1)"
          
          modelstring <- paste0(
            "Y[j] ~ dbeta(alpha[j], beta[j])
            alpha[j] <- mu[j] * phi
            beta[j] <- (1 - mu[j]) * phi
            logit( mu[j] ) <- beta0[i]*indMea[i] + mu.beta0*(1-indMea[i]) 
            r[j]=Y[j] - mu[j]
        }
            " )
          
          finalstring <- paste0(
            "model{
            for(i in 1:(Nset-1) ){
            for(j in offsets[i]:(offsets[i+1]-1)){
            ",modelstring, mixpriorstr,"
            }
            ", priorstr, "}")
          }
        }
    }
    else{
    
    
    
    #####################
    ######Non-mixed Model
    #####################
    if( mixed == F ){
        priorstr <- paste0( "beta0 ~ dnorm(0, 0001) ", paste0(  betastr, " ~ dnorm(0, 0.001)", collapse = " "), " phi ~ dgamma(",par.a, ", ", par.b,")" )
        
        modelstring <- paste0(
        "Y[j] ~ dbeta(alpha[j], beta[j])
        alpha[j] <- mu[j] * phi
        beta[j] <- (1 - mu[j]) * phi
        logit( mu[j] ) <- beta0 + ", parstr,"
        r[j]=Y[j] - mu[j]
    }
    ")
    
    finalstring <- paste0(
    "model{
    for(j in 1:N ){
    ", modelstring, priorstr, "}")
    
}
#############################
######Mixed Model#############
################################
    else{
    print("run mixed models")
    
    if( is.na(indMea) == T ){
        priorstr <- paste0( "mu.beta0 ~ dnorm(0, 0.001) tau.beta0 ~ dgamma(1,1) ",
        paste0(  betastr, " ~ dnorm(0, 0.001)", collapse = " "), " phi ~ dgamma(1, 1)" )
        mixpriorstr <- "beta0[i] ~ dnorm(mu.beta0, tau.beta0) "
        
        modelstring <- paste0(
        "Y[j] ~ dbeta(alpha[j], beta[j])
        alpha[j] <- mu[j] * phi
        beta[j] <- (1 - mu[j]) * phi
        logit( mu[j] ) <- beta0[i] + ", parstr,"
        r[j]=Y[j] - mu[j]
    }
    " )
    
    finalstring <- paste0(
    "model{
    for(i in 1:(Nset-1) ){
        for(j in offsets[i]:(offsets[i+1]-1)){
            ",modelstring, mixpriorstr,"
        }
    ", priorstr, "}")
}
else{
    priorstr <- paste0( "mu.beta0 ~ dnorm(0, 0.001) ",
    paste0(  betastr, " ~ dnorm(0, 0.001)", collapse = " "), " phi ~ dgamma(1, 1)" )
    mixpriorstr <- "beta0[i] ~ dnorm(mu.beta0, tau.beta0[i]) tau.beta0[i] ~ dgamma(1,1)"
    
    modelstring <- paste0(
    "Y[j] ~ dbeta(alpha[j], beta[j])
    alpha[j] <- mu[j] * phi
    beta[j] <- (1 - mu[j]) * phi
    logit( mu[j] ) <- beta0[i]*indMea[i] + mu.beta0*(1-indMea[i]) + ", parstr,"
    r[j]=Y[j] - mu[j]
}
" )

finalstring <- paste0(
"model{
for(i in 1:(Nset-1) ){
    for(j in offsets[i]:(offsets[i+1]-1)){
        ",modelstring, mixpriorstr,"
    }
", priorstr, "}")
}
}
}
    
return(finalstring)
}





autoBayes <- function( formula, data, offsets = NA, indMea = NA, DIC = F,
n.adapt = 200, n.iter = 200, n.chain = 1,
phi.a = 1, phi.b = 1){
    
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match( c("formula", "data"), names(mf), 0L)
    mf <- mf[ c(1L, m)]
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "numeric")
    X <- model.matrix(mt, mf)
    
    X <- X[,-1]
    colnames(X) <- gsub(":", ".", colnames(X))
    
    ##Setting
    n.var <- ncol(X)
    n.obs <- length(Y)
    
    factorCols <- which( sapply( X,  is.factor) )
    
    if(DIC == T){
        if( n.chain < 2)
        n.chain = 2
    }
    
    
    if( length(factorCols) != 0 ){
        print("please assign proper coding to factor variables")
    }
    else{
        varname <- colnames(X)
        betastr <- paste0( "b.", varname )
        Xstr <- paste0("X.", varname)
        
        varlist <- lapply( 1:length(varname),  function(i) return(X[,i])   )
        names(varlist) <- Xstr
    }
    
    
    
    
    if( length(offsets) == 1 ) {
        mixed = F
        jags.data <- c( list( N = n.obs, Y = Y), varlist   )
    }else{
        mixed = T
        jags.data <- c( list( Nset=length(offsets), offsets=offsets,  indMea=indMea, Y = Y), varlist )
    }
    
    ########################
    #######fitting#######
    ########################
    jags.script <- inistring(varname, mixed , indMea = indMea, par.a = phi.a, par.b = phi.b)
    cat( jags.script )
    jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
    
    
    if( length(offsets) == 1 ){
        jags.res <- coda.samples(jags.mod, c(betastr, 'mu', 'r'), n.iter = n.iter)
    }else{
        jags.res <- coda.samples(jags.mod, c('mu.beta0', betastr, 'mu', 'r'), n.iter = n.iter)
    }
    
    ###R square processing###
    mustring <- rep('mu', n.obs )
    rstring <- rep('r', n.obs)
    
    mustring <- sapply( 1:n.obs, function(i) paste0( mustring[i], '[', i, ']' )   )
    rstring <- sapply( 1:n.obs, function(i) paste0( rstring[i], '[', i, ']' )   )
    
    
    mu.r.val <- summary( jags.res[[1]][, c(mustring,rstring) ] )
    mean.mu <-  mu.r.val[[1]][mustring, 1]
    mean.r <-  mu.r.val[[1]][rstring,1]
    
    SSfit <- var(mean.mu)
    SSres <- var(mean.r)
    
    
    r2 <- SSfit/(SSfit + SSres)
    
    
    penalizedDIC <- NA
    if(DIC == T){
        dic.beta <- dic.samples(jags.mod, n.iter = n.iter)
        penalizedDIC <- sum(dic.beta$deviance) + sum(dic.beta$penalty)
    }
    
    
    
    return( list( 
    jags.res = jags.res,
    n.obs = n.obs,
    r2 = r2,
    jags.script = jags.script,
    jags.data = jags.data,
    varlist = varlist,
    betastr = betastr,
    varname = varname,
    X=X, Y=Y, offsets = offsets, chainsetting = c(n.adapt, n.iter, n.chain),
    fitted = mean.mu, residual = mean.r,
    DIC = penalizedDIC,
    phi.a = phi.a, phi.b = phi.b, indMea = indMea))
}




step0 <- function (obj, direction = "backward",  maxstep = 1000, criterion = "R2"){
    
    Y <- obj$Y
    X <- obj$X
    n.obs <- length(Y)
    indMea <- obj$indMea
    
    n.adapt <- obj$chainsetting[1]; n.iter <- obj$chainsetting[2]; n.chain <- obj$chainsetting[3]
    varlist <- obj$varlist
    
    offsets <- obj$offsets
    mixed <- ifelse( length(offsets) == 1, F, T)
    
    #####process string##########
    #########################
    varname <- obj$varname
    betastr <- paste0( "b.", varname )
    
    mustring <- rep('mu', n.obs )
    rstring <- rep('r', n.obs)
    
    mustring <- sapply( 1:n.obs, function(i) paste0( mustring[i], '[', i, ']' )   )
    rstring <- sapply( 1:n.obs, function(i) paste0( rstring[i], '[', i, ']' )   )
    
    
    ###inital no-confident invertal sign
    jag.sum <- summary(a$jags.res[[1]])
    
    
    allCI <- jag.sum[[2]]
    signdiff <- apply(sign( allCI[betastr, c(1,5)] ), 1, diff )
    finaldiff<- length( which(signdiff != 0) )

    
    #################################
    ##############################
    
    jags.res <- list()
    if( mixed == F ) {
        jags.data <- c( list( N = n.obs, Y = Y), varlist   )
    }else{
        jags.data <- c( list( Nset=length(offsets), offsets=offsets,  indMea=indMea, Y = Y), varlist )
    }
    
    interID <- which( grepl(".", colnames(a$X), fixed = T)  )
    noninterID <- which( !grepl(".", colnames(a$X), fixed = T)  )

    
    if( direction == "backward"){
        if( criterion == "R2"){
          max.r2.prev <- obj$r2
          steptrack <- 1
          while( steptrack <= maxstep  ){
            r2 <- c()
            difflen <- c()
            
            if( length(varname) > 1 ){
              hiera_stop = 0
              
              while( length(interID) > 0 ){
                
                r2 <- c()
                difflen <- c()
                interID <- which( grepl(".", varname, fixed = T)  )
      
                
                for( i in 1:length(interID)  ){
                  currentID <- interID[i]
                  
                  cat("interID:", interID, "\n")
                  cat("currentID:", currentID, "\n")
                  cat("varname:", varname, "\n")
                  
                  jags.script <- inistring( varname[ -currentID ], mixed, indMea, obj$phi.a, obj$phi.b )
                  jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
                  temp.name <- paste0('- ', varname[ currentID ] )
                  
                  if ( mixed == F )
                    jags.res <-  coda.samples(jags.mod, c('beta0', betastr[-currentID], 'mu', 'r'), n.iter = n.iter)[[1]]
                  else
                    jags.res <-  coda.samples(jags.mod, c('mu.beta0', betastr[-currentID], 'mu', 'r'), n.iter = n.iter)[[1]]
                  
                  
                  jag.sum <- summary( jags.res )
                  ######get r2#########
                  
                  mean.mu <-  jag.sum[[1]][mustring, 1]
                  mean.r <-  jag.sum[[1]][rstring, 1]
                  
                  SSfit <- var(mean.mu); SSres <- var(mean.r)
                  
                  temp.r2 <- SSfit/(SSfit + SSres)
                  r2 <- rbind(r2, temp.r2)
                  rownames(r2)[i] <- temp.name
                  
                  
                  ######get sign#######
                  allCI <- jag.sum[[2]]
                  print(betastr[-currentID])
                  signdiff <- apply(sign( allCI[betastr[-currentID], c(1,5)] ), 1, diff )
                  temp.diff<- length( which(signdiff != 0) )
                  
                  difflen <- rbind(difflen, temp.diff)
                  
                  
                }
                  
                  colnames(r2) <- "R2"
                  
                  max.id <- which.max(r2)
                  print(varname[interID])
                  print(r2)
                  
                  print(difflen)
                  
                  if( r2[max.id] >= max.r2.prev ){
                    
                    varname <-  varname[ -interID[max.id] ]
                    betastr <- betastr[ -interID[max.id] ]
                    max.r2.prev <- r2[max.id]
                    interID <- interID[-max.id]
                    finaldiff <- difflen[max.id]
                  }
                  else{
                    hiera_stop = 1
                    print("Delete any variable won't lead a better model!")
                    break
                  }
                  
              }
              
              if( hiera_stop == 1 )
                break
              
              
              for( i in 1:length(varname) ){
                jags.script <- inistring(varname[-i], mixed, indMea, obj$phi.a, obj$phi.b)
                jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
                temp.name <- paste0('- ', varname[i])
                
                if ( mixed == F )
                  jags.res <-  coda.samples(jags.mod, c('beta0', betastr[-i], 'mu', 'r'), n.iter = n.iter)[[1]]
                else
                  jags.res <-  coda.samples(jags.mod, c('mu.beta0', betastr[-i], 'mu', 'r'), n.iter = n.iter)[[1]]
                
                
                jag.sum <- summary( jags.res )
                ###  R square processing###
                mean.mu <-  jag.sum[[1]][mustring, 1]
                mean.r <-  jag.sum[[1]][rstring, 1]
                
                SSfit <- var(mean.mu)
                SSres <- var(mean.r)
                
                temp.r2 <- SSfit/(SSfit + SSres)
                r2 <- rbind(r2, temp.r2)
                rownames(r2)[i] <- temp.name
                
                
                ######get sign#######
                allCI <- jag.sum[[2]]
                
                if( length(betastr[-i]) > 1 )
                  signdiff <- apply( sign( allCI[betastr[-i], c(1,5)] ), 1, diff )
                else
                  signdiff <- diff(sign( allCI[betastr[-i], c(1,5)] ))
                
                temp.diff<- length( which(signdiff != 0) )
                
                difflen <- rbind(difflen, temp.diff)
              }
            }
            else{
              print("Only one variable left in the model")
              jags.script <- inistring(varname, mixed, indMea, obj$phi.a, obj$phi.b)
              jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
              if ( mixed == F )
                jags.res <-  coda.samples(jags.mod, c('beta0', betastr, 'mu', 'r'), n.iter = n.iter)[[1]]
              else
                jags.res <-  coda.samples(jags.mod, c('mu.beta0', betastr, 'mu', 'r'), n.iter = n.iter)[[1]]
              
              
              
              jag.sum <- summary( jags.res )
              ###  R square processing###
              mean.mu <-  jag.sum[[1]][mustring, 1]
              mean.r <-  jag.sum[[1]][rstring,1]
              
              SSfit <- var(mean.mu)
              SSres <- var(mean.r)
              
              temp.r2 <- SSfit/(SSfit + SSres)
              r2 <- rbind(r2, temp.r2)
              
              ######get sign#######
              allCI <- jag.sum[[2]]
              signdiff <- sign( allCI[betastr, c(1,5)] )
              temp.diff<- length( which(signdiff != 0) )
              
              final.diff <- temp.diff
              
              break
            }
            
            
            colnames(r2) <- "R2"
            
            max.id <- which.max(r2)
            print(varname)
            print(r2)
            
            if( r2[max.id] >= max.r2.prev ){
              varname <-  varname[-max.id]
              betastr <- betastr[-max.id]
              max.r2.prev <- r2[max.id]
              finaldiff <- difflen[max.id]
            }
            else{
              print("Delete any variable won't lead a better model!")
              break
            }
          }
        }
      
        if( criterion == "DIC"){
          min.dic.prev <- obj$DIC
          cat("initial DIC", min.dic.prev, "\n")
          
          steptrack <- 1
          
          
          while( steptrack <= maxstep  ){
            dic <- c()
            difflen <- c()
            
            if( length(varname) > 1 ){
              hiera_stop = 0
              
              while( length(interID) > 0 ){
                
                dic <- c()
                difflen <- c()
                interID <- which( grepl(".", varname, fixed = T)  )
                
                
                for( i in 1:length(interID)  ){
                  currentID <- interID[i]
                  
                  cat("interID:", interID, "\n")
                  cat("currentID:", currentID, "\n")
                  cat("varname:", varname, "\n")
                  
                  jags.script <- inistring( varname[ -currentID ], mixed, indMea, obj$phi.a, obj$phi.b )
                  jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
                  temp.name <- paste0('- ', varname[ currentID ] )
                  
                  if ( mixed == F )
                    jags.res <-  coda.samples(jags.mod, c('beta0', betastr[-currentID ], 'mu', 'r'), n.iter = n.iter)[[1]]
                  else
                    jags.res <-  coda.samples(jags.mod, c('mu.beta0', betastr[-currentID ], 'mu', 'r'), n.iter = n.iter)[[1]]
                  
                  ###  R square processing###
                  dic.beta <- dic.samples(jags.mod, n.iter = n.iter)
                  temp.dic <- sum(dic.beta$deviance) + sum(dic.beta$penalty)
                  
                  dic <- rbind(dic, temp.dic)
                  rownames(dic)[i] <- temp.name
                  
                  
                  ###diffsign
                  jag.sum <- summary( jags.res )
                  allCI <- jag.sum[[2]]
                  print(betastr[-currentID])
                  signdiff <- apply(sign( allCI[betastr[-currentID], c(1,5)] ), 1, diff )
                  temp.diff<- length( which(signdiff != 0) )
                  
                  difflen <- rbind(difflen, temp.diff)
                  
                }
                
                colnames(dic) <- "DIC"
                
                min.id <- which.min(dic)
                print(varname[interID])
                print(dic)
                
                if( dic[min.id] <= min.dic.prev ){
                  
                  varname <-  varname[ -interID[min.id] ]
                  betastr <- betastr[ -interID[min.id] ]
                  min.dic.prev <- dic[min.id]
                  interID <- interID[-min.id]
                  finaldiff <- difflen[min.id]
                }
                else{
                  hiera_stop = 1
                  print("Delete any variable won't lead a better model!")
                  break
                }
                
              }
              
              if( hiera_stop == 1 )
                break
              
              
              
              for( i in 1:length(varname) ){
                jags.script <- inistring(varname[-i], mixed, indMea, obj$phi.a, obj$phi.b)
                jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
                temp.name <- paste0('- ', varname[i])
                
                if ( mixed == F )
                  jags.res <-  coda.samples(jags.mod, c('beta0', betastr[-i], 'mu', 'r'), n.iter = n.iter)[[1]]
                else
                  jags.res <-  coda.samples(jags.mod, c('mu.beta0', betastr[-i], 'mu', 'r'), n.iter = n.iter)[[1]]
                
                ###  R square processing###
                dic.beta <- dic.samples(jags.mod, n.iter = n.iter)
                temp.dic <- sum(dic.beta$deviance) + sum(dic.beta$penalty)
                
                dic <- rbind(dic, temp.dic)
                rownames(dic)[i] <- temp.name
                
                jag.sum <- summary( jags.res )
                if( length(betastr[-i]) > 1 )
                  signdiff <- apply( sign( allCI[betastr[-i], c(1,5)] ), 1, diff )
                else
                  signdiff <- diff(sign( allCI[betastr[-i], c(1,5)] ))
                
                temp.diff<- length( which(signdiff != 0) )
                
                difflen <- rbind(difflen, temp.diff)
              }
            }
            else{
              print("Only one variable left in the model")
              jags.script <- inistring(varname, mixed, indMea, obj$phi.a, obj$phi.b)
              cat( jags.script )
              jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
              if ( mixed == F )
                jags.res <-  coda.samples(jags.mod, c('beta0', betastr, 'mu', 'r'), n.iter = n.iter)[[1]]
              else
                jags.res <-  coda.samples(jags.mod, c('mu.beta0', betastr, 'mu', 'r'), n.iter = n.iter)[[1]]
              
              ###  R square processing###
              
              dic.beta <- dic.samples(jags.mod, n.iter = n.iter)
              temp.dic <- sum(dic.beta$deviance) + sum(dic.beta$penalty)
              
              dic <- rbind(dic, temp.dic)
              
              jag.sum <- summary( jags.res )
              allCI <- jag.sum[[2]]
              signdiff <- sign( allCI[betastr, c(1,5)] )
              temp.diff<- length( which(signdiff != 0) )
              
              finaldiff <- temp.diff
              
              break
            }
            
            
            colnames(dic) <- "DIC"
            
            min.id <- which.min(dic)
            print(varname)
            print(dic)
            
            if( dic[min.id] <= min.dic.prev ){
              varname <-  varname[-min.id]
              betastr <- betastr[-min.id]
              max.dic.prev <- dic[min.id]
              finaldiff <- difflen[min.id]
            }
            else{
              print("Delete any variable won't lead a better model!")
              break
            }
          }
        }
      
      return( list( varname = varname, finaldiff = finaldiff ) )
    }
    
    
    if( direction == "forward"){
      if( criterion == "R2"){
        
        max.r2.prev <- 0
        fittedvar <- c()
        fittedbeta <- c()
        
        steptrack <- 1
        hiera_stop = 0
        while( steptrack <= maxstep  ){
          r2 <- c()
          difflen <- c()
            
            while( length(noninterID) > 0 ){
              r2 <- c()
              difflen <- c()
              noninterID <- which( !grepl(".", varname, fixed = T)  )
              
            for( i in 1:length(noninterID)  ){
              currentID <- noninterID[i]
              
              cat("noninterID:", noninterID, "\n")
              cat("currentID:", currentID, "\n")
              cat("varname:", varname, "\n")
              
              jags.script <- inistring( c(fittedvar, varname[currentID]), mixed, indMea, obj$phi.a, obj$phi.b)
              jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
              temp.name <- paste0('+', varname[currentID])
              
              if ( mixed == F )
                jags.res <-  coda.samples(jags.mod, c('beta0', fittedbeta, betastr[currentID], 'mu', 'r'), n.iter = n.iter)[[1]]
              else
                jags.res <-  coda.samples(jags.mod, c('mu.beta0', fittedbeta, betastr[currentID], 'mu', 'r'), n.iter = n.iter)[[1]]
              
              
              jag.sum <- summary( jags.res )
              ###  R square processing###
              mean.mu <-  jag.sum[[1]][mustring, 1]
              mean.r <-  jag.sum[[1]][rstring,1]
              
              SSfit <- var(mean.mu)
              SSres <- var(mean.r)
              
              temp.r2 <- SSfit/(SSfit + SSres)
              r2 <- rbind(r2, temp.r2)
              rownames(r2)[i] <- temp.name
              
              
              ### sign
              allCI <- jag.sum[[2]]
              
              if( length(c(fittedbeta, betastr[currentID])) > 1 )
                signdiff <- apply( sign( allCI[ c(fittedbeta, betastr[currentID]), c(1,5)] ), 1, diff )
              else
                signdiff <- diff(sign( allCI[ c(fittedbeta, betastr[currentID]), c(1,5)] ))
              
              temp.diff<- length( which(signdiff != 0) )
              difflen <- rbind(difflen, temp.diff)
              
              
            }
              
              colnames(r2) <- "R2"
              
              max.id <- which.max(r2)
              cat("fittedvar", fittedvar, "\n")
              print(r2)
              
              
              if( r2[max.id] >= max.r2.prev ){
                fittedvar <-c(fittedvar, varname[noninterID[max.id]])
                fittedbeta <- paste0("b.", fittedvar)
                varname <-  varname[ -noninterID[max.id] ]
                betastr <- betastr[ -noninterID[max.id] ]
                max.r2.prev <- r2[max.id]
                noninterID <- noninterID[-max.id]
                finaldiff <- difflen[max.id]
              }
              else{
                hiera_stop = 1
                
                
                print("Add any variable won't lead a better model!")
                break
              }
          }
          
            if( hiera_stop == 1 )
              break
              
            
          r2 <- c()
          if( length(varname) > 1 ){
            for( i in 1:length(varname) ){
              jags.script <- inistring( c(fittedvar, varname[i]), mixed, indMea, obj$phi.a, obj$phi.b)
              jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
              temp.name <- paste0('+', varname[i])
              
              if ( mixed == F )
                jags.res <-  coda.samples(jags.mod, c('beta0', c(fittedbeta, betastr[i]), 'mu', 'r'), n.iter = n.iter)[[1]]
              else
                jags.res <-  coda.samples(jags.mod, c('mu.beta0', c(fittedbeta, betastr[i] ), 'mu', 'r'), n.iter = n.iter)[[1]]
              
              
              jag.sum <- summary( jags.res )
              ###  R square processing###
              mean.mu <-  jag.sum[[1]][mustring, 1]
              mean.r <-  jag.sum[[1]][rstring,1]
              
              SSfit <- var(mean.mu)
              SSres <- var(mean.r)
              
              temp.r2 <- SSfit/(SSfit + SSres)
              r2 <- rbind(r2, temp.r2)
              rownames(r2)[i] <- temp.name
              
              ### sign
              allCI <- jag.sum[[2]]
              signdiff <- apply(sign( allCI[c(fittedbeta, betastr[i]), c(1,5)] ), 1, diff )
              temp.diff<- length( which(signdiff != 0) )
              
              difflen <- rbind(difflen, temp.diff)
            }
          }else if( length(varname) == 1 ){
            print("Only one more variable could be fitted")
            jags.script <- inistring(c(fittedvar, varname), mixed, indMea, obj$phi.a, obj$phi.b)
            jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
            if ( mixed == F )
              jags.res <-  coda.samples(jags.mod, c('beta0', c(fittedbeta, betastr), 'mu', 'r'), n.iter = n.iter)[[1]]
            else
              jags.res <-  coda.samples(jags.mod, c('mu.beta0', c(fittedbeta, betastr), 'mu', 'r'), n.iter = n.iter)[[1]]
            
            
            jag.sum <- summary( jags.res )
            ###  R square processing###
            mean.mu <-  jag.sum[[1]][mustring, 1]
            mean.r <-  jag.sum[[1]][rstring,1]
            
            SSfit <- var(mean.mu)
            SSres <- var(mean.r)
            
            temp.r2 <- SSfit/(SSfit + SSres)
            r2 <- rbind(r2, temp.r2)
            rownames(r2) <- varname
            
            ### sign
            allCI <- jag.sum[[2]]
            signdiff <- apply(sign( allCI[c(fittedbeta, betastr), c(1,5)] ), 1, diff )
            temp.diff<- length( which(signdiff != 0) )
            
            difflen <- temp.diff
            
          }else{
            break
          }
          
          colnames(r2) <- "R2"
          
          max.id <- which.max(r2)
          cat("fittedvar", fittedvar, "\n")
          print(r2)
          
          if( r2[max.id] >= max.r2.prev ){
            
            fittedvar <-c(fittedvar, varname[max.id])
            fittedbeta <- paste0("b.", fittedvar)
            varname <-  varname[-max.id]
            betastr <- betastr[-max.id]
            max.r2.prev <- r2[max.id]
            finaldiff <- difflen[max.id]
          }
          else{
            print("Add any variable won't lead a better model!")
            break
          }
        }
        
      }
      if( criterion == "DIC"){
        
        min.dic.prev <- Inf
        fittedvar <- c()
        fittedbeta <- c()
        cat("initial DIC", min.dic.prev, "\n")
        
        steptrack <- 1
        hiera_stop = 0
        while( steptrack <= maxstep  ){
          dic <- c()
          difflen <- c()
          
          while( length(noninterID) > 0 ){
            dic <- c()
            difflen <- c()
            noninterID <- which( !grepl(".", varname, fixed = T)  )
            
            for( i in 1:length(noninterID)  ){
              currentID <- noninterID[i]
              
              cat("noninterID:", noninterID, "\n")
              cat("currentID:", currentID, "\n")
              cat("varname:", varname, "\n")
              
              jags.script <- inistring( c(fittedvar, varname[currentID]), mixed, indMea, obj$phi.a, obj$phi.b)
              jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
              temp.name <- paste0('+', varname[currentID])
              
              if ( mixed == F )
                jags.res <-  coda.samples(jags.mod, c('beta0', betastr[i], 'mu', 'r'), n.iter = n.iter)[[1]]
              else
                jags.res <-  coda.samples(jags.mod, c('mu.beta0', betastr[i], 'mu', 'r'), n.iter = n.iter)[[1]]
              
              ### DIC processing###
              dic.beta <- dic.samples(jags.mod, n.iter = n.iter)
              temp.dic <- sum(dic.beta$deviance) + sum(dic.beta$penalty)
              
              dic <- rbind(dic, temp.dic)
              rownames(dic)[i] <- temp.name
              
              ##
              allCI <- jag.sum[[2]]
              print(c(fittedbeta, betastr[currentID]))
              
              if( length(c(fittedbeta, betastr[currentID])) > 1 )
                signdiff <- apply( sign( allCI[ c(fittedbeta, betastr[currentID]), c(1,5)] ), 1, diff )
              else
                signdiff <- diff(sign( allCI[ c(fittedbeta, betastr[currentID]), c(1,5)] ))
              
              temp.diff<- length( which(signdiff != 0) )
              difflen <- rbind(difflen, temp.diff)
              
            }
            
            colnames(dic) <- "DIC"
            
            min.id <- which.min(dic)
            cat("fittedvar", fittedvar, "\n")
            print(dic)

            
            
            if( dic[min.id] < min.dic.prev ){
              fittedvar <-c(fittedvar, varname[noninterID[min.id]])
              fittedbeta <- paste0("b.", fittedvar)
              varname <-  varname[ -noninterID[min.id] ]
              betastr <- betastr[ -noninterID[min.id] ]
              min.dic.prev <- dic[min.id]
              noninterID <- noninterID[-min.id]
              finaldiff <- difflen[min.id]
            }
            else{
              hiera_stop = 1
              print("Add any variable won't lead a better model!")
              break
            }
          }
          
          if( hiera_stop == 1 )
            break
          
          
          dic <- c()
          if( length(varname) > 1 ){
            for( i in 1:length(varname) ){
              jags.script <- inistring( c(fittedvar, varname[i]), mixed, indMea, obj$phi.a, obj$phi.b)
              jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
              temp.name <- paste0('+', varname[i])
              
              if ( mixed == F )
                jags.res <-  coda.samples(jags.mod, c('beta0', betastr[i], 'mu', 'r'), n.iter = n.iter)[[1]]
              else
                jags.res <-  coda.samples(jags.mod, c('mu.beta0', betastr[i], 'mu', 'r'), n.iter = n.iter)[[1]]
              
              ###  R square processing###
              dic.beta <- dic.samples(jags.mod, n.iter = n.iter)
              temp.dic <- sum(dic.beta$deviance) + sum(dic.beta$penalty)
              
              dic <- rbind(dic, temp.dic)
              rownames(dic)[i] <- temp.name
              
              ### sign
              allCI <- jag.sum[[2]]
              signdiff <- apply(sign( allCI[c(fittedbeta, betastr[i]), c(1,5)] ), 1, diff )
              temp.diff<- length( which(signdiff != 0) )
              difflen <- rbind(difflen, temp.diff)
              
            }
          }else if( length(varname) == 1 ){
            print("Only one more variable could be fitted")
            jags.script <- inistring(c(fittedvar, varname), mixed, indMea, obj$phi.a, obj$phi.b)
            jags.mod <- jags.model(textConnection(jags.script), data = jags.data, n.adapt = n.adapt, n.chain = n.chain)
            if ( mixed == F )
              jags.res <-  coda.samples(jags.mod, c('beta0', betastr, 'mu', 'r'), n.iter = n.iter)[[1]]
            else
              jags.res <-  coda.samples(jags.mod, c('mu.beta0', betastr, 'mu', 'r'), n.iter = n.iter)[[1]]
            
            ###  R square processing###
            dic.beta <- dic.samples(jags.mod, n.iter = n.iter)
            temp.dic <- sum(dic.beta$deviance) + sum(dic.beta$penalty)
            
            dic <- rbind(dic, temp.dic)
            rownames(dic) <- varname
            
            ### sign
            allCI <- jag.sum[[2]]
            signdiff <- apply(sign( allCI[c(fittedbeta, betastr), c(1,5)] ), 1, diff )
            temp.diff<- length( which(signdiff != 0) )
            difflen <- temp.diff
            
          }else{
            break
          }
          
          colnames(dic) <- "DIC"
          
          min.id <- which.min(dic)
          cat("fittedvar", fittedvar, "\n")
          print(dic)
          
          if( dic[min.id] < min.dic.prev ){
            
            fittedvar <-c(fittedvar, varname[min.id])
            fittedbeta <- paste0("b.", fittedvar)
            varname <-  varname[-min.id]
            betastr <- betastr[-min.id]
            min.dic.prev <- dic[min.id]
            finaldiff <- difflen[min.id]
          }
          else{
            print("Add any variable won't lead a better model!")
            break
          }
        }
      }
      
      return( list( varname = fittedvar, finaldiff = finaldiff ) )
    }
}


