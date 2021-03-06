rm(list = ls())

library("psych")
library("rjags")


setwd("/Users/Qingyan/Documents/BUaca/research/APOE")
apoe.data <- read.csv("blessed.data.csv")
edu.data <- read.csv("edu.csv")


###Data preprocessing
apoe.data <- apoe.data[ apoe.data$Blessed.Missing %in% NA, ]
apoe.mem <- apoe.data[complete.cases( apoe.data$Mem.Total.Score), ]
apoe.mem$edu <- edu.data[ match( apoe.mem$ID, edu.data$ID ), 3 ]
apoe.mem$edu <- as.integer( levels( apoe.mem$edu )[apoe.mem$edu] )
apoe.mem.edu <- apoe.mem[ -which( is.na(apoe.mem$edu) == TRUE ), ] 
#remove 1 obs
apoe.mem.edu <- apoe.mem.edu[ !is.na( apoe.mem.edu$Age.at.Enrollment), ]
#remove 1 obs
apoe.mem.edu <- apoe.mem.edu[ apoe.mem.edu$Memory.Test.Age >= apoe.mem.edu$Age.at.Enrollment, ]

#####summary
cat("total observation:", nrow(apoe.mem.edu))
cat("total patient:", length( unique( apoe.mem.edu$ID ) ) )

##Avoid 0 and 1 in Beta regression
max(apoe.mem.edu$Mem.Total.Score); min(apoe.mem.edu$Mem.Total.Score)
apoe.mem.edu$Mem.Total.Score[ which( apoe.mem.edu$Mem.Total.Score == 37 ) ] <- 36
apoe.mem.edu$Mem.Total.Score[ which( apoe.mem.edu$Mem.Total.Score == 0 ) ] <- 1



##
apoe.e2vse3 <- apoe.mem.edu[ -which( apoe.mem.edu$APOE == "e2.e4" | apoe.mem.edu$APOE == "e3.e4" ), ]
offset.e2vse3 <- match( unique( apoe.e2vse3$ID), apoe.e2vse3$ID )
offset.e2vse3 <- append( offset.e2vse3, length(apoe.e2vse3$ID)+1)
apoe.e2vse3$APOE <- factor(apoe.e2vse3$APOE)

par(mfrow=c(1,1))
barMem <- barplot(table(apoe.e2vse3$APOE),ylim = c(-10, 1000), main = "Mem.score Apoe Barplot")
text(x = barMem, y = table(apoe.e2vse3$APOE), label = table(apoe.e2vse3$APOE) )



###
indMea <- numeric( length(offset.e2vse3) - 1 )
for( i in 1:( length(offset.e2vse3)-1 ) ){
  tempID <- offset.e2vse3[i]:(offset.e2vse3[i+1]-1)  
  if( length(tempID) == 1){
  }
  else
    indMea[i] = 1
}




data.e2vse3 <- list( Nset=length(offset.e2vse3), offset=offset.e2vse3, indMea=indMea,
                              Y=apoe.e2vse3$Mem.Total.Score/37, X.age=apoe.e2vse3$Memory.Test.Age, 
                              X.age.bar=mean(apoe.e2vse3$Memory.Test.Age ),  X.edu=apoe.e2vse3$edu, X.edu.bar=mean(apoe.e2vse3$edu),
                              X.sex=as.numeric(apoe.e2vse3$Sex)-1, X.apoe=dummy.code( apoe.e2vse3$APOE ) )


model.beta.e2vse3 <-
  "model{
for(i in 1:(Nset-1) ){
for(j in offset[i]:(offset[i+1]-1)){

Y[j] ~ dbeta(alpha[j], beta[j])
alpha[j] <- mu[j] * phi 
beta[j] <- (1 - mu[j]) * phi
logit( mu[j] ) <- beta0[i]*indMea[i] + mu.beta0*(1-indMea[i]) + b.age*(X.age[j]-X.age.bar) + b.sex*X.sex[j] + b.edu*(X.edu[j]-X.edu.bar) + b.age2*(X.age[j]-X.age.bar)^2 + 
b.age.edu*(X.age[j]-X.age.bar)*(X.edu[j]-X.edu.bar) + b.e2*( X.apoe[j,1]||X.apoe[j,2] ) + 
b.e2in*(X.apoe[j,1]||X.apoe[j,2])*(X.age[j]-X.age.bar)
r[j] = Y[j] - mu[j]
}


beta0[i] ~ dnorm(mu.beta0, tau.beta0)

}


mu.beta0 ~ dnorm(0, 0.001)
b.age ~ dnorm(0, 0.001)
b.sex ~ dnorm(0, 0.001)
b.edu ~ dnorm(0, 0.001)
b.age2 ~ dnorm(0, 0.001)
b.age.edu ~ dnorm(0, 0.001)
b.e2 ~ dnorm(0, 0.001)
b.e2in ~ dnorm(0, 0.001)
tau.beta0 ~ dgamma(1, 1)
phi ~ dgamma(1, 1)
}"


jags.beta.e2vse3 <- jags.model(textConnection(model.beta.e2vse3), data = data.e2vse3, n.adapt = 4000, n.chains = 3)
update(jags.beta.e2vse3, 4000)
test.beta.e2vse3 <- coda.samples(jags.beta.e2vse3, c('mu.beta0', 'b.age', 'b.sex', 'b.edu', 'b.age2', 'b.age.edu', 'b.e2', 'b.e2in'), n.iter=4000)


par(mfrow=c(2,2))
summary( test.beta.e2vse3[[1]] )
traceplot( test.beta.e2vse3[[1]] )
autocorr.plot(test.beta.e2vse3[[1]])
gelman.plot(test.beta.e2vse3)



###########################e4 vs e3################
apoe.e4vse3 = apoe.mem.edu[ -which( apoe.mem.edu$APOE == "e2.e2" | apoe.mem.edu$APOE == "e2.e3" | apoe.mem.edu$APOE == "e2.e4" ), ]
apoe.e4vse3$APOE <- factor( apoe.e4vse3$APOE )

offset.e4vse3 <- match( unique(apoe.e4vse3$ID), apoe.e4vse3$ID )
offset.e4vse3 <- append( offset.e4vse3, length(apoe.e4vse3$ID)+1)


par(mfrow=c(1,1))
barMem <- barplot(table(apoe.e4vse3$APOE),ylim = c(-10, 1000), main = "Mem.score Apoe Barplot")
text(x = barMem, y = table(apoe.e4vse3$APOE), label = table(apoe.e4vse3$APOE) )


##
###
indMea <- numeric( length(offset.e4vse3) - 1 )
for( i in 1:( length(offset.e4vse3)-1 ) ){
  tempID <- offset.e4vse3[i]:(offset.e4vse3[i+1]-1)  
  if( length(tempID) == 1){
  }
  else
    indMea[i] = 1
}


data.e4vse3 <- list( Nset=length(offset.e4vse3), offset=offset.e4vse3, indMea=indMea,
                              Y=apoe.e4vse3$Mem.Total.Score/37, X.age=apoe.e4vse3$Memory.Test.Age, 
                              X.age.bar=mean(apoe.e4vse3$Memory.Test.Age ),  X.edu=apoe.e4vse3$edu, X.edu.bar=mean(apoe.e4vse3$edu),
                              X.sex=as.numeric(apoe.e4vse3$Sex)-1, X.apoe=dummy.code( apoe.e4vse3$APOE ) )

model.beta.e4vse3 <-
  "model{
for(i in 1:(Nset-1) ){
for(j in offset[i]:(offset[i+1]-1)){

Y[j] ~ dbeta(alpha[j], beta[j])
alpha[j] <- mu[j] * phi 
beta[j] <- (1 - mu[j]) * phi
logit( mu[j] ) <- beta0[i]*indMea[i] + mu.beta0*(1-indMea[i]) + b.age*(X.age[j]-X.age.bar) + b.sex*X.sex[j] + b.edu*(X.edu[j]-X.edu.bar) + b.age2*(X.age[j]-X.age.bar)^2 + 
b.age.edu*(X.age[j]-X.age.bar)*(X.edu[j]-X.edu.bar) + b.e4*X.apoe[j,2]  + b.e4in*X.apoe[j,2]*(X.age[j]-X.age.bar)
r[j] = Y[j] - mu[j]
}

beta0[i] ~ dnorm(mu.beta0, tau.beta0)

}


mu.beta0 ~ dnorm(0, 0.001)
b.age ~ dnorm(0, 0.001)
b.sex ~ dnorm(0, 0.001)
b.edu ~ dnorm(0, 0.001)
b.age2 ~ dnorm(0, 0.001)
b.age.edu ~ dnorm(0, 0.001)
b.e4 ~ dnorm(0, 0.001)
b.e4in ~ dnorm(0, 0.001)
tau.beta0 ~ dgamma(1, 1)
phi ~ dgamma(1, 1)
}"

jags.beta.e4vse3 <- jags.model(textConnection(model.beta.e4vse3), data = data.e4vse3, n.adapt = 4000, n.chains = 3)
update(jags.beta.e4vse3, 4000)
test.beta.e4vse3 <- coda.samples(jags.beta.e4vse3, c('mu.beta0', 'b.age', 'b.sex', 'b.edu', 'b.age2', 'b.age.edu', 'b.e4', 'b.e4in'), n.iter=4000)



par(mfrow=c(2,2))
summary( test.beta.e4vse3[[1]] )
traceplot( test.beta.e4vse3[[1]] )
autocorr.plot(test.beta.e4vse3[[1]])
gelman.plot(test.beta.e4vse3)
