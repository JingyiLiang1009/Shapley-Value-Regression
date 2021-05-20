# Shapley Value Regression.
#
# This is a package for relative importance calculation based on Shapley Value Regression.
#
#' @param y is the data set of dependent variable
#' @param x is the data set of independent variables
#
# In this package, users do not need to point out the number of independent variables, it will be caught automatically while inputting x.
#


error2 <- function(message) {
  stop(simpleError(message))
}

shapleyvalue <- function(y,x){
  k <- ncol(x)
  times_wt <- NULL
  for (i in 1:k-1){
    times_wt[i] <- choose(k-1,i)
  }
  data <- data.frame(matrix(ncol=k,nrow=2^(k-1)))
  name <- colnames(x)

  if (k==1){
    error2("No need for using shapley regression!")
  }
  else if (k==2){
    for (i in 1:k){
      n <- name[i]
      lm <- lm(y~x[,n])
      lm2 <- lm(y~x[,-i])
      lm_all <- lm(y~.,data=x)
      r1 <- summary(lm)$r.squared
      r2 <- summary(lm2)$r.squared
      r3 <- summary(lm_all)$r.squared
      data[1,i] <- r1
      data[2,i] <- r3-r2
    }
    weight <- c(1/2,1/2)
    data <- as.matrix(data)
    dat <- weight%*%data
    dat <- dat %>% as.data.frame()
    sum <- rowSums(dat)
    uni_value <- dat/sum
    dat <- rbind(round(dat,4),round(uni_value,4))
    colnames(dat) <- name
    rownames(dat) <- c("Shapley Value","Standardized Shapley Value")
    return(dat)
  }
  else{
    for (i in 1:k){
      num <- 1
      for (j in 1:k){
        if(j==1){
          n <- name[i]
          lm <- lm(y~x[,n])
          data[1,i] <- summary(lm)$r.squared
        }
        else if(j==k){
          xx <- combn(k, k)
          xlab <- NULL
          for(l in 1:j){
            x_lab[l] <- xx[l]
            if(x_lab[l]==i){
              x1 <- x[,x_lab[l]]
              x_reslab <- xx[-l]
              x_res <- x[x_reslab]
              lm <- lm(y~.,data=x_res)
              r <- summary(lm)$r.squared
              lm2 <- lm(y~.+x1,data=x_res)
              r2 <- summary(lm2)$r.squared
            }
          }
          data[2^(k-1),i] <- r2-r
        }
        else{
          times <- choose(k-1,j-1)
          times2 <- choose(k,j)
          combine <- combn(k,j)
          n <- NULL
          nn <- combine
          for(m in 1:j){
            nn[m,] <- combine[m,]==i
          }
          for(m in 1:times2){
            n[m] <- 1 %in% nn[,m]
          }
          xx <- combine[,which(n)]
          for (o in 1:times){
            num <- num+1
            x_lab <- NULL
            for(l in 1:j){
              x_lab[l] <- xx[,o][l]
              if(x_lab[l]==i){
                x1 <- x[,x_lab[l]]
                x_reslab <- xx[,o][-l]
                x_res <- x[x_reslab]
                lm <- lm(y~.,data=x_res)
                r <- summary(lm)$r.squared
                lm2 <- lm(y~.+x1,data=x_res)
                r2 <- summary(lm2)$r.squared
              }
            }
            data[num,i] <- r2-r
          }
        }
      }
    }
    weight <- 1/k
    for (q in 1:k-1){
      weight <- c(weight,rep(1/k/times_wt[q],times_wt[q]))
    }
    data <- as.matrix(data)
    weight <- weight %>% as.matrix() %>% t()
    dat2 <- weight%*%data
    dat2 <- dat2 %>% as.data.frame()
    sum <- rowSums(dat2)
    uni_value <- dat2/sum
    dat2 <- rbind(round(dat2,4),round(uni_value,4))
    colnames(dat2) <- name
    rownames(dat2) <- c("Shapley Value","Standardized Shapley Value")
    return(dat2)
  }
}
