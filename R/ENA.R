
ENA <- function(mgmModel, data, factors, activity, symptoms){
  NM_new <- mgmModel

  p <- ncol(data)
  # Apply thresholding
  for(i in 1:p) {
    NM_i <- NM_new[[i]]
    if(data_mgm_type[i]=="c") {
      NM_i$model$`1`[abs(NM_i$model$`1`) < NM_i$tau] <- 0
    } else {
      NM_i$model[abs(NM_i$model) < NM_i$tau] <- 0
    }
    NM_new[[i]] <- NM_i
  }


  # Apply AND-rule
  for(i in 1:p) {
    NM_i <- NM_new[[i]]
    if(data_mgm_type[i]=="c") {
      NM_i$model$`1`[-1,1][FullModel$pairwise$wadj[-i, i] == 0] <- 0
    } else {
      NM_i$model[-1,1][FullModel$pairwise$wadj[-i, i] == 0] <- 0
    }
    NM_new[[i]] <- NM_i
  }

  new_means <- NULL
  data_means <- apply(data, 2, mean) #means from data

  for (i in 1:ncol(data)){

    scaled_betas <- as.numeric(NM_new[[i]]$model)
  }

  for (i in 1:ncol(data)){

    unscaled_betas <- scaled_betas[-1] * sd(data[,i]) / sd(as.numeric(unlist(data[,-i])))
    unscaled_intercept <- data_means[i] - sum(data_means[-i] * unscaled_betas)
    new_means[i] <- unscaled_intercept + sum(unscaled_betas * data_means[-i])
  }

  #plot(data_means, new_means)

  round(data_means, 10) == round(new_means, 10)

  data_means_conditioned_protective <- data_means
  data_means_conditioned_protective[c(factors)] <- activity

  conditioned_means_protective <- NULL

  for (i in 1:ncol(data)){
    if(data_mgm_type[i] == "c"){
      scaled_betas <- as.numeric(NM_new[[i]]$model$'2')
    }else {
      scaled_betas <- as.numeric(NM_new[[i]]$model)
    }
    unscaled_betas <- scaled_betas[-1] * sd(data[,i]) / sd(as.numeric(unlist(data[,-i])))
    unscaled_intercept <- data_means[i] - sum(data_means[-i] * unscaled_betas)
    conditioned_means_protective[i] <- unscaled_intercept + sum(unscaled_betas * data_means_conditioned_protective[-i])
  }

  bESA <- sum(new_means[c(symptoms)])
  fESA <- sum(conditioned_means_protective[c(symptoms)])

  dESA <- sum(new_means[c(symptoms)]) - sum(conditioned_means_protective[c(symptoms)])

  res <- list(bESA, fESA, dESA)

  ESA <- list("baseline ENA" = bESA, "simulation ENA" = fESA, "diference ENA" = dESA)

  return(ESA)

}
