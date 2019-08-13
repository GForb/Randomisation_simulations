##################################################
## Project: OPTYC
##
## Script purpose: 
##  Simulations to evaluate different randomsiation options
##
## Author: Gordon Forbes
##################################################

library(tidyverse)
library(randomizeR)
library(Minirand)
library(furrr)
library(tictoc)
plan(multiprocess)

## Section: Funtions implemting randomisation methods   ##################################################

# Simple randomisation
rand_simple <- function(patients, trt_seq, ratio) {
  n <- nrow(patients)
  patients %>% mutate(allocation = sample(trt_seq, n, replace = TRUE, prob = ratio/sum(ratio)) %>% 
                                    factor(levels = trt_seq), 
                      rand_method = "simple")
}

#Random permuted block
rand_perm_block <- function(patients, rand_param, trt_seq) {
  patients %>% mutate(rand_method = "perm_block", allocation = genSeq(rand_param) %>%  
    getRandList %>% 
    as.character %>% 
   factor(levels = trt_seq))
}

# Funtions for stratified permuted block randomisation
make_strata <- function(data, ...) {
  vars <- enquos(...)
  strata_index <- data %>% 
    select(!!!vars) %>% 
    group_by(!!!vars) %>% 
    summarise_all(count) %>% 
    rowid_to_column(".strata")
  inner_join(data, strata_index) %>% arrange(id)
}  

rand_one_strata <- function(patients, strat_var, strat_no, rand_param, trt_seq) {
  strat_var <- enquo(strat_var)
  strata <- patients %>% filter(!!strat_var == strat_no)
  rand_seq <- genSeq(rand_param) %>%  getRandList %>% as.character %>% .[1:nrow(strata)]
  strata %>% mutate(rand_method = "strat_perm_block", allocation = factor(rand_seq, levels = trt_seq))
}

rand_strat_perm_block <- function(patients, rand_param, trt_seq, strat_vars) {
  patients <- make_strata(patients, !!!strat_vars)
  strat_nos <- patients %>% select(.strata) %>% distinct  %>% .[[1]] 
  map_dfr(strat_nos, ~rand_one_strata(patients = patients,
                                    strat_var = .strata,
                                    strat_no = .,
                                    trt_seq = trt_seq,
                                    rand_param = rand_param)) %>% 
    select(-.strata)

}


# Minimization
rand_minimise <- function(patients, min_vars, trt_seq, ratio , covwt = rep(1/ncol(covmat), ncol(covmat)), random_component = 0.8) {
  covmat <- patients %>% select(!!!min_vars)
  nsample <- nrow(covmat)
  result <- rep(100, nsample) 
  result[1] = sample(trt_seq, 1, replace = TRUE, prob = ratio/sum(ratio))
  for (j in 2:nsample) {
    # get treatment assignment sequentiall for all subjects
    result[j] <- Minirand(covmat=covmat, 
                          j, 
                          covwt=covwt, 
                          ratio=ratio,
                          ntrt=length(trt_seq), 
                          trtseq=trt_seq, 
                          method="Range", 
                          result=result, 
                          p = random_component)
    
  }
  mutate(patients, rand_method = "minimisation", allocation = result)
}



rand_minimise_uneq_ratio <- function(patients, min_vars, trt_seq, ratio , covwt = rep(1/ncol(covmat), ncol(covmat)), random_component = 0.8) {
  n_arms = sum(ratio)
  breaks = c(0, cumsum(ratio))
  fake_trt_seq = 1:n_arms
  fake_ratio = rep(1, n_arms)
  patients %>% 
    rand_minimise(min_vars, fake_trt_seq, fake_ratio, random_component) %>% 
    mutate(allocation = cut(allocation, breaks, labels = trt_seq))
}

## Section: Funtions: Implementing simulation ##################################################
create_patients <- function(n) {
  tibble(id = 1:!!n, 
         gender = rbinom(!!n, 1, 0.5),   # 50/50 gender
         severity = rbinom(!!n, 1, 0.25),      # 75/25 severity
         site = cut(runif(!!n), c(0, 0.4, 0.8, 1), labels = FALSE) -1,  # Only 20% recrutied from thrid site
         referal = cut(runif(!!n), c(0, 0.5, 0.8, 1), labels = FALSE) -1) # Imbalance in referal routes
}

run_ranomisations <- function(data, trt_seq, ratio, rand_param, strat_vars, min_vars, random_component = 0.8) {
  bind_rows(
    data %>% rand_simple(trt_seq, ratio),
    data %>% rand_perm_block(rand_param, trt_seq),
    data %>% rand_strat_perm_block(rand_param, trt_seq, strat_vars = vars(severity, referal)),
    data %>% rand_minimise_uneq_ratio(min_vars =  vars(severity, referal, site, gender),
                                     trt_seq = trt_seq, 
                                     ratio = ratio, 
                                    random_component = random_component)
  ) 
  
}

analyse_sim <- function(allocations, true) {
  counts <- allocations %>% group_by(rand_method, allocation) %>% summarise(n = n())
  
  allocation_summaries <- allocations %>% select(-id) %>% group_by(rand_method, allocation) %>% summarise_all(mean) %>%
    inner_join(counts) %>% select(rand_method, allocation, n, everything())
  
}



sim_rep <- function(trt_seq, ratio, rand_param, strat_vars, min_vars, n, random_component = 0.8) {
  allocations <- run_ranomisations(create_patients(n),
                                   trt_seq,
                                   ratio,
                                   rand_param, 
                                   strat_vars = strat_vars,
                                   min_vars =  min_vars,
                                   random_component)
  analyse_sim(allocations, true)
}

sim_analysis <- function(allocation_summaries, true, sum_vars) {
  allocation_differences <- allocation_summaries %>% 
    mutate(n = case_when(allocation == "C" ~ n - true["n"]/2,
                         allocation != "C" ~ n - true["n"]), 
           gender = gender - true["gender"],
           severity = severity - true["severity"],
           referal = referal - true["referal"],
           site = site - true["site"])
  
  max_differences <- allocation_differences %>% group_by(sim_no, rand_method) %>% select(-allocation) %>% summarise_all(~max(abs(.)))
  medians <- max_differences  %>%  group_by(rand_method) %>% summarize_at(sum_vars, median) %>% mutate(percentile = 50)
  p75 <- max_differences %>% group_by(rand_method) %>% summarize_at(sum_vars, ~quantile(.,probs = 0.75)) %>% mutate(percentile = 75)
  p95 <- max_differences %>% group_by(rand_method) %>% summarize_at(sum_vars, ~quantile(.,probs = 0.95)) %>% mutate(percentile = 95)
  p99 <- max_differences %>% group_by(rand_method) %>% summarize_at(sum_vars, ~quantile(.,probs = 0.99)) %>% mutate(percentile = 99)
  bind_rows(medians, p75, p95, p99) %>% select(percentile, everything())
  
}

## Section: Simulation   ##################################################

# Setting up random permuted block structure
n <- 60
trt_seq <- c("A", "B", "C")
ratio <- c(2,2,1)
block_sizes <- c(5, 10)
rand_param <-  rpbrPar(n, block_sizes, K = length(ratio), ratio = ratio , filledBlock = FALSE) 

true <- c(n = n*0.4, gender = 0.5, severity = 0.25, site = 0.8, referal = 0.7)


#Sims
tic()
set.seed(46416)
sim_results <- future_map_dfr(1:5000, ~sim_rep(trt_seq = trt_seq, 
                              ratio = ratio, 
                              rand_param = rand_param, 
                              strat_vars = vars(severity, referal), 
                              min_vars = vars(gender, severity, site, referal),
                              n = n,
                              random_component = 0.9), .progress = TRUE, .id = "sim_no")
toc()

#Saving
paste("rand_sim", gsub(":", "-", Sys.time()), ".rdata", sep="") %>% 
  save(sim_results, file = .)

## Section: Processing simulation results   ##################################################



sum_vars <-  vars(n, gender, severity, site, referal)
analysed_results <- sim_analysis(sim_results, true, sum_vars) %>% 
  mutate(rand_method = factor(rand_method, labels = c("Minimisation", 
                                                      "Rand. Permuted block", 
                                                      "Simple",
                                                      "Strat. permuted block")),
                              percentile = factor(percentile, labels =c("Median", 
                                                                       "75th",
                                                                       "95th",
                                                                       "99th")))

analysed_results_long <- analysed_results %>% gather(key = "variable", value = "imbalance", !!!sum_vars)

## Section: Plotting results   ##################################################

theme_set(theme_light(base_size = 15))


  analysed_results_long %>% filter(variable == "n") %>% ggplot(aes(x = imbalance, y = percentile, colour = rand_method)) +
    geom_point(size = 3) +
    labs(x = "Max. imbalance", y = "Centile of simulated data sets", colour = "Randomisation method", title = "Imbalance in group size") +
    scale_colour_brewer(type = "div", palette = 6)
  
  analysed_results_long %>% 
    filter(variable != "n" & (percentile == "Median" | percentile == "99th")) %>% 
    ggplot(aes(x = imbalance, y = variable, colour = rand_method)) +
    geom_jitter(size = 2, width = 0, height = 0.2) +
    labs(x = "Max. imbalance", colour = "Randomisation method", title = "Imbalance in covariates") +
    scale_colour_brewer(type = "div", palette = 6) +
    facet_grid(rows = vars(percentile))



# note stratified permuted blocks is on severity and referal
