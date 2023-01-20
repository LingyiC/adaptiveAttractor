# 011223

library("cafr")

## ===
## Compared to the original version, the algorithm stops when the attractor diverged.
## ===
findAttractor.2 <- function(data, seed, a = 5, seedTopN = 100, maxIter = 100, epsilon = 1e-7, bin = 6, 
    so = 3, rankBased = FALSE, negateMI = TRUE, verbose = TRUE) 
{   
    vec <- data[seed, ]
    dataIn <- data
    if (rankBased) {
        vec <- rank(vec)
        for (i in 1:nrow(data)) {
            dataIn[i, ] <- rank(data[i, ])
        }
    }
    mi <- getAllMIWz(dataIn, vec, bin = bin, so = so, negateMI = negateMI)
    premi <- mi
    w <- abs(mi)^a/sum(abs(mi)^a)
    w[mi < 0] <- 0
    metagene <- w %*% data
    if (rankBased) 
        metagene <- rank(metagene)
    c <- 0
    while (c < maxIter) {
        mi <- getAllMIWz(dataIn, metagene, bin = bin, so = so, 
            negateMI = negateMI)
        delta <- sum((mi - premi)^2)
        if (verbose) {
            cat("Iteration ", (c + 1), "\tDelta = ", delta, "\n", 
                sep = "")
            print(mi[order(mi, decreasing = T)[1:20]])
            flush.console()
        }
        if (delta < epsilon) {
            break
        }
        if(which(names(sort(mi, decreasing = T)) == seed) > seedTopN){
          if (verbose) {
            cat("Iteration ", (c + 1), "\tDelta = ", delta, "\n", "Diverged!!!",
                            sep = "")
          }
          break
        }
        premi <- mi
        mi[mi < 0] <- 0
        w <- abs(mi)^a/sum(abs(mi)^a)
        w[mi < 0] <- 0
        metagene <- w %*% data
        if (rankBased) 
            metagene <- rank(metagene)
        c <- c + 1
    }
    if (c >= maxIter) 
        return(NULL)
    return(sort(mi, decreasing = T))
}

## 011223 version In this version, we make sure the adaptive can work using any exponent.max/exponent.min settings

findAttractor.adaptive <- function(data, seed, exponent.max = 5, exponent.min = 2, step.large = 1, step.small = 0.1, dominantThreshold = 0.2, seedTopN = 100, compareNth = 10, verbose = FALSE, epsilon = 1e-7, minimize = FALSE, first.idx = 1, second.idx = 2){
  genesOrder <- rownames(data)
  history_matrix <- matrix(data = NA, nrow = nrow(data), ncol = 1)
  rownames(history_matrix) <- genesOrder
  history_exponent <- c()
  
## initial
  i <- 1
  a = exponent.max
  tmp <- findAttractor.2(data, seed, a, epsilon = epsilon, seedTopN = seedTopN, verbose = verbose)

  if(is.null(tmp) == TRUE){
    history_matrix[, 1] <- rep(0, length(genesOrder))
  }else{
    history_matrix[, 1] <- tmp[genesOrder] 
  }
    
  
  if(is.null(tmp)){
      Failed = FALSE
      notConverge = TRUE
      diverged = FALSE
      increaseStop = FALSE 
      history_exponent <- a
      if((exponent.max-exponent.min) >= step.large){ 
        a <- a - step.large
        skipRound1 = FALSE 
      }else{ # skip round 2 
        a <- a - step.small
        skipRound1 = TRUE
      }
      
  }else if(which(names(tmp) == seed) > seedTopN){
      Failed = TRUE # if attractor diverged using a = exponent.max, then we skip round 1 and round 2. 
      notConverge = FALSE
      diverged = TRUE 
      increaseStop = FALSE    
      skipRound1 = FALSE
      history_exponent <- a
  }else {
      Failed = FALSE
      notConverge = FALSE
      diverged = FALSE
      increaseStop = FALSE   
      history_exponent <- a
      if((exponent.max-exponent.min) >= step.large){
        a <- a - step.large
        skipRound1 = FALSE
      }else if((exponent.max-exponent.min) < step.large && (exponent.max-exponent.min) > 0 ){ # skip round 1
        exponent.first <- exponent.min - step.small #it skipped round 1, so we make the exponent.first = exponent.min - step.small
        a <- a - step.small
        skipRound1 = TRUE
        skipRound2 = FALSE
      }else if(exponent.max == exponent.min){
        final.exponent <- exponent.max
        skipRound1 = TRUE
        skipRound2 = TRUE        
      }
  }



## Round 1, step.large 
cat("1st round start... \n")
  while(Failed == FALSE && skipRound1 == FALSE && a >= exponent.min ){
    tmp <- findAttractor.2(data, seed, a, epsilon = epsilon, seedTopN = seedTopN, verbose = verbose)
    if(is.null(tmp)){
      notConverge = TRUE
      diverged = NULL
      increaseStop = NULL
      exponent.first <- a
      a <- a + step.large - step.small
      skipRound2 = FALSE
      break
    }
    if(which(names(tmp) == seed) > seedTopN){
      notConverge = FALSE
      diverged = TRUE
      increaseStop = NULL
      exponent.first <- a
      a <- a + step.large - step.small
      skipRound2 = FALSE
      break
    }

    if(minimize == FALSE){
      if(tmp[compareNth] < sort(history_matrix[, i], decreasing = T)[compareNth]){
        notConverge = FALSE
        diverged = FALSE
        increaseStop = TRUE
        exponent.first <- a
        a <- a + step.large - step.small
        skipRound2 = FALSE
        break
      }}else if(minimize == TRUE){
      if((tmp[second.idx] - tmp[first.idx])< (sort(history_matrix[, i], decreasing = T)[second.idx] - sort(history_matrix[, i], decreasing = T)[first.idx])){
        notConverge = FALSE
        diverged = FALSE
        increaseStop = TRUE
        exponent.first <- a
        a <- a + step.large - step.small
        skipRound2 = FALSE
        break
      }      
    }
    i <- i + 1
    history_matrix <- cbind(history_matrix, tmp[genesOrder])
    history_exponent <- c(history_exponent, a)
    a <- a - step.large

    if(a < exponent.min && (a + step.large) != exponent.min ){
      increaseStop = NULL
      diverged = NULL
      notConverge = NULL
      exponent.first <- exponent.min - step.small #it skipped round 1, so we make the exponent.first = exponent.min - step.small
      a <- a + step.large - step.small
      skipRound2 = FALSE
      break
    }else if((a + step.large) == exponent.min){
      increaseStop = NULL
      diverged = NULL
      notConverge = NULL
      exponent.first <- exponent.min
      skipRound2 = TRUE
      final.exponent <- exponent.min
    }
  }

cat("1st round end. \n")



## Round 2, step.small 
cat("2nd round start... \n")
  if(Failed == FALSE && skipRound2 == FALSE && (skipRound1 = TRUE || diverged == TRUE || increaseStop == TRUE || notConverge == TRUE) ){
    while(a > exponent.first  && a < exponent.max){
      
      tmp <- findAttractor.2(data, seed, a, epsilon = epsilon, seedTopN = seedTopN, verbose = verbose)
      if(is.null(tmp)){
        notConverge = TRUE
        diverged = NULL
        increaseStop = NULL
        final.exponent <- a + step.small
        break
      }
      if(which(names(tmp) == seed) > seedTopN){
        diverged = TRUE
        increaseStop = NULL
        notConverge = NULL
        final.exponent <- a + step.small
        break
      }
      if(minimize == FALSE){
        if(tmp[compareNth] < sort(history_matrix[, i], decreasing = T)[compareNth]){
          increaseStop = TRUE
          diverged = NULL
          notConverge = NULL
          final.exponent <- a + step.small
          break
        }
      }else if(minimize == TRUE){
        if((tmp[second.idx] - tmp[first.idx])< (sort(history_matrix[, i], decreasing = T)[second.idx] - sort(history_matrix[, i], decreasing = T)[first.idx])){
          increaseStop = TRUE
          diverged = NULL
          notConverge = NULL
          final.exponent <- a + step.small
          break
        }      
      }
      i <- i + 1
      history_matrix <- cbind(history_matrix, tmp[genesOrder])
      history_exponent <- c(history_exponent, a)
      a <- a - step.small
      if(a <= exponent.first){
        final.exponent <- a + step.small
      }
    }  
  }
    
cat("2nd round end. \n")

## OUTPUT
  if(Failed == FALSE){
      # if(a < exponent.min){
      #   final.exponent <- exponent.min
      # }
      attractor.final <- sort(history_matrix[, ncol(history_matrix)], decreasing = T)

      if(attractor.final[1] - attractor.final[2] < dominantThreshold){ ## Dominant threshold
          coexpressed <- names(attractor.final[1:compareNth])
      }else{
          coexpressed <- c()
      }
      scannedSeeds <- unique(c(names(which(attractor.final >= attractor.final[seed])), coexpressed))

      return(list(attractor.final = attractor.final, scannedSeeds = scannedSeeds, history_matrix = history_matrix, history_exponent = history_exponent, final.exponent = final.exponent))
  }else if(Failed == TRUE){
      attractor.final <- rep(0, length(genesOrder))
      names(attractor.final) <- genesOrder
      return(list(attractor.final = attractor.final, scannedSeeds = seed, history_matrix = history_matrix, history_exponent = history_exponent, final.exponent = NULL))
  }


  
}


