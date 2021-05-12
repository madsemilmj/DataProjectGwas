#' Assign_LTFH function
#'
#' This function assigns phenotypes based on LTFH article. This implementation uses Gibs-sampling to efficiently assign genetic liabilities to each individual.
#' @param Pheno_data Is the Phenotype data that comes out of the simulation when running SimulerData
#' @param valT is the dignificance level, usually 0.05
#' @param h2 Is the .. usually 0.5
#' @param with_sib Binary, defines wheter or not we should include siblings. Defaults to 1 = yes.
#' @param Burn_in Defines the number of samples to discard in the gibs-sampling. Defaults to 100 (found by simulating convergence tests)
#' @keywords Assign LTFH
#' @export
#' @importFrom dplyr %>%
#' @examples
#' Assign_LTFH(Pheno_data = tibble::tibble(ID = 1:100, Child = rbinom(100,1,0.3), Mom = rbinom(100,1,0.3), Dad = rbinom(100,1,0.3), Nr_sib = rep(0,100), Sib_status = rep(0,100)), valT=0.05, h2=0.5, with_sib = 1, Burn_in = 100)

Assign_LTFH <- function(Pheno_data, valT, h2, with_sib = 1, Burn_in = 100){
  #profvis::profvis({
  Pheno_data <- Pheno_data %>%
    tibble::as_tibble()%>%
    dplyr::select(ID, Child, Mom, Dad, Nr_sib, Sib_status)
  if (with_sib == 0){
    Pheno_data$Nr_sib <- 0
    Pheno_data$Sib_status <- 0
  }

  Combos <- Pheno_data %>%
    tibble::as_tibble() %>%
    dplyr::distinct(Child, Mom, Dad, Nr_sib, Sib_status)

  sib_status_comb <- Combos %>%
    dplyr::distinct(Nr_sib,Sib_status) %>%
    as.matrix()

  Combos <- as.matrix(Combos)

  for (i in 1:nrow(sib_status_comb)){
    row <- sib_status_comb[i,]
    ##If c(0,1,0) for that sib_status_comb is present delete c(0,0,1) for that comb
    if (nrow(as.matrix(Combos[which(Combos[,1] == 0 & Combos[,2] == 1 & Combos[,3] == 0 & Combos[,4] == row[1] & Combos[,5] == row[2]),])) > 0 ){
      if (nrow(as.matrix(Combos[which(Combos[,1] == 0 & Combos[,2] == 0 & Combos[,3] == 1 & Combos[,4] == row[1] & Combos[,5] == row[2]),])) > 0 ){
        Combos <- Combos[-which(Combos[,1] == 0 & Combos[,2] == 0 & Combos[,3] == 1 & Combos[,4] == row[1] & Combos[,5] == row[2]),]
      }
    }
    ##If c(1,1,0) for that sib_status_comb is present delete c(1,0,1) for that comb
    if (nrow(as.matrix(Combos[which(Combos[,1] == 1 & Combos[,2] == 1 & Combos[,3] == 0 & Combos[,4] == row[1] & Combos[,5] == row[2]),])) > 0){
      if (nrow(as.matrix(Combos[which(Combos[,1] == 1 & Combos[,2] == 0 & Combos[,3] == 1 & Combos[,4] == row[1] & Combos[,5] == row[2]),])) > 0){
        Combos <- Combos[-which(Combos[,1] == 1 & Combos[,2] == 0 & Combos[,3] == 1 & Combos[,4] == row[1] & Combos[,5] == row[2]),]
      }
    }
  }

  means <- numeric(nrow(Combos))
  SEs <- numeric(nrow(Combos))
  SDs <- numeric(nrow(Combos))
  for(j in 1:nrow(Combos)){
    row <- Combos[j,]
    nr_sib <- row[4]
    case_sib <- row[5]
    LowerBound1 <- ifelse(row[1:3]>0, qnorm(1-valT), -Inf)
    UpperBound1 <- ifelse(row[1:3]>0, Inf, qnorm(1-valT))
    LowerBound2 <- rep(-Inf, nr_sib)
    UpperBound2 <- rep(qnorm(1-valT), nr_sib)
    c <- case_sib
    while (c > 0){
      LowerBound2[c] <- qnorm(1-valT)
      UpperBound2[c] <- Inf
      c <- c-1
    }

    LowerBound <- c(LowerBound1,LowerBound2)
    UpperBound <- c(UpperBound1,UpperBound2)
    CovMatrix <- LinearTransformation(CovarianceMatrix(h2,nr_sib))
    CondInfo <- ConditionalDistribution(CovMatrix)
    sigma <- CondInfo$sigma
    mu <- CondInfo$mu
    #Drawing from truncated
    S <- 6000
    my_sample <- matrix(ncol = (length(LowerBound)+1), nrow = S)
    vals <- rep(1,(length(LowerBound)+1))
    for (i in 1:S){
      #Child gen
      vals[1] <- rnorm(1,mu[1,] %*% vals[2:length(vals)], sqrt(sigma[1]))
      #REST
      for (p in 2:(length(LowerBound)+1)){
        mu_prod <- mu[p,] %*% vals[-p]
        l1 <- pnorm(LowerBound[p-1], mu_prod, sqrt(sigma[p]))
        u1 <- pnorm(UpperBound[p-1], mu_prod, sqrt(sigma[p]))
        g1 <- runif(1,l1,u1)
        vals[p] <- qnorm(g1,mu_prod, sqrt(sigma[p]))
      }
      my_sample[i,] <- vals
    }
    real_sample <- my_sample[Burn_in:S,]
    eps_g_mean <- mean(real_sample[,1])
    SEs[j] <- sem(real_sample[,1])
    means[j] <- eps_g_mean
    SDs[j] <- sd(real_sample[,1])
  }

  Ny_combos <- cbind(Combos,means,SDs)%>%
    tibble::as_tibble()

  Ny_pheno <- dplyr::left_join(tibble::as_tibble(Pheno_data),Ny_combos,by = c("Child"="Child", "Mom"= "Mom", "Dad" = "Dad", "Nr_sib" = "Nr_sib", "Sib_status" = "Sib_status"))

  #ADDING THOSE THAT WE DELETED
  Missing <- Ny_pheno %>%
    dplyr::filter(is.na(means))%>%
    dplyr::distinct(Child, Mom, Dad, Nr_sib, Sib_status) %>%
    dplyr::rename_with(tolower)%>%
    dplyr::rowwise() %>%
    dplyr::mutate(means = dplyr::filter(Ny_combos, Child == child & Mom == mom+1 & Dad == dad-1 & Nr_sib == nr_sib & Sib_status == sib_status)$means)%>%
    dplyr::mutate(SDs = dplyr::filter(Ny_combos, Child == child & Mom == mom+1 & Dad == dad-1 & Nr_sib == nr_sib & Sib_status == sib_status)$SDs)
  result <- dplyr::left_join(Ny_pheno,Missing,by = c("Child"="child", "Mom"= "mom", "Dad" = "dad", "Nr_sib" = "nr_sib", "Sib_status" = "sib_status")) %>%
    dplyr::mutate(means = dplyr::coalesce(means.x, means.y)) %>%
    dplyr::mutate(SDs = dplyr::coalesce(SDs.x, SDs.y)) %>%
    dplyr::mutate(FID = ID) %>%
    dplyr::rename(Pheno = means) %>%
    dplyr::rename(IID = ID) %>%
    dplyr::select(FID,IID,Child, Mom, Dad, Nr_sib, Sib_status, Pheno, SDs)
  #})
  return(result)
}
