library(tidyverse)
library(lmtest)
library(doParallel)
registerDoParallel(cores=10)

# Load SNES and drop na's
# SNES can be downloaded from https://www.kaggle.com/datasets/parsabg/stocknewseventssentiment-snes-10
data <- readr::read_csv('./data-augmented.csv')
data <- data %>% drop_na()

# Granger test orders
# Note: Higher order tests (e.g. 30 and above) run significantly slower than lower order ones
test_orders <- c(1, 3, 7, 14, 30)

# Minimum standard deviation for a time series to be considered in the granger test
# Saves compute time by skipping low variance time series
min_SD <- 1

# Fetch tickers from SNES
tickers <- data %>%
  distinct(Symbol) %>%
  pull(Symbol)

# Fetch event types from SNES (used in computing the granger matrix)
event_types <- data %>% 
  select(c(matches("News")),-`News - Mergers and Acquisitions`, `Adj Close`, `Volume`) %>% 
  colnames()

# Helper function for differencing a time series to make it stationary
difference <- function(v) {
  v <- v - lag(v)
  v[!is.na(v)]
}

# Multithreaded granger computation (one ticker per CPU core)
output <- foreach(i=1:length(tickers), .combine = rbind, .verbose=TRUE) %dopar% {
  # Create a dataframe containing data about the current ticker only
  ticker <- tickers[i]
  ticker_df <- data %>% filter(Symbol == ticker)
  
  # Create an empty dataframe to append test results to
  res <- data.frame(Symbol=character(0),
                    order=numeric(0),
                    event_type_1=character(0),
                    event_type_2=character(0),
                    Res.Df=double(0),
                    Df=double(0),
                    F=double(0),
                    P=double(0))
  
  # Run pairwise granger test for each test order and all event types within this ticker
  for (order in test_orders) {
    for (event_type_i in event_types) {
      ts_i <- ticker_df[[event_type_i]]
      
      # Skip if ts_i is empty
      if(length(ts_i) == 0) next
      # Skip if ts_i is low variance
      if(sd(ts_i) < min_SD) next
      
      # Replace ts_i with its differenced version
      ts_i <- difference(ts_i)
      
      for (event_type_j in event_types) {
        if(event_type_i != event_type_j) { # Don't test with self
          ts_j <- ticker_df[[event_type_j]]
          
          # Skip if ts_j is empty
          if(length(ts_j) == 0) next
          # Skip if ts_j is low variance
          if(sd(ts_j) < min_SD) next
          
          # Replace ts_j with its differenced version
          ts_j <- difference(ts_j)
          
          tryCatch({
            print(paste("Computing granger causality of order [",order,"] between [",event_type_i,"] and [",event_type_j,"] for [",ticker,"]"))
            granger <- grangertest(ts_i ~ ts_j, order = order)
            res[nrow(res) + 1,] <- c(ticker, order, event_type_i, event_type_j, 
                     granger[[1]][2], granger[[2]][2], granger[[3]][2], granger[[4]][2])
          }
          , error = function(e) {print("An error occurred while running granger test")}
          )
        }
      }
    }
  }
  
  if(nrow(res) > 0) {
    # Fix column types (for some reason in the loop above numeric and dbl columns convert to chr)
    res$order <- as.double(res$order)
    res$Res.Df <- as.double(res$Res.Df)
    res$Df <- as.double(res$Df)
    res$F <- as.double(res$F)
    res$P <- as.double(res$P)
    
    # Write a CSV file with all pairwise granger results for this ticker
    write.csv(res, paste("./causal-data-differenced/",ticker,".csv",sep=""))
  }
}

# Read all individual CSVs and combine them into a single dataframe
causal_csvs <- list.files("causal-data-differenced", "*.csv")
causal_data <- lapply(paste("causal-data-differenced/", causal_csvs, sep=""), readr::read_csv, col_types = "ncnccdddd")
causal_data <- bind_rows(causal_data)

# Add GICS sector and sub industries to the final dataframe
ticker_industry <- data %>% 
  select(Symbol, `GICS Sector`, `GICS Sub-Industry`) %>% 
  distinct()
causal_data <- causal_data %>% 
  inner_join(ticker_industry)

# Write final data as a single CSV
write.csv(causal_data, "causal_data_differenced.csv")
