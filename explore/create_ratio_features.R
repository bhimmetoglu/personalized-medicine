## Engineer features for Variants Data
# Author(s): B. Himmetoglu (burakhmmtgl@gmail.com)
# 9/20/2017

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(caret)

# Fill NA's with default ratio
fill_na <- function(x, default_ratio){
  x[is.na(x)] <- default_ratio
  x
}

## Create n_Var_Class/n_Var on folds
# For a given class, compute : 
# (number of the given variation in that class) / (number of occurances of the given variation)
# This ratio is computed on out-of-sample fashion to prevent overfitting
create_ratios_folds <- function(train, nfolds, seed=1234){
  # Create n folds
  set.seed(seed)
  folds <- createFolds(train$Class, k = nfolds)
  
  # Initiate a list where holdout data frames from each fold will be stored
  list_out <- list()
  
  # Loop over folds 
  for (ifold in 1:length(folds)){
    train_in <- train[-folds[[ifold]],]
    train_ho <- train[folds[[ifold]],]
    
    # Group by gene Variation and Class
    df1 <- train %>% 
      group_by(Variation, Class) %>%
      summarise(n_Var_Class = n())
    
    # Group by Variation
    df2 <- train %>% 
      group_by(Variation) %>%
      summarise(n_Var = n())
    
    # Join and calculate ratio
    df3 <- left_join(df1,df2, by = "Variation")
    df3 <- df3 %>% mutate(freq = n_Var_Class/n_Var)
    
    # Spread the ratios to columns
    df3$n_Var <- NULL # Not needed
    df3$n_Var_Class <- NULL # Not needed
    df4 <- spread(df3, key = Class, value = freq, fill = 0)
    
    # Join df4 with hold-out data
    train_ho <- left_join(train_ho, 
                          df4 %>% ungroup(Variation),
                          by = "Variation")
    
    # Impute NAs with 0
    train_ho <- train_ho %>% mutate_at(names(df4)[-1], fill_na, default_ratio = 0.0)
    
    # Add to list
    list_out[[ifold]] <- train_ho
  }
  
  # Return train_out
  bind_rows(list_out)
}