rm(list=ls())

library(maxLik)
library(foreach)
library(doParallel)
library(dplyr)

# Function to constrain behavior within range ----

est_constrain <- function(x){
  
  if (x > 90){
    u <- 90
  } else if (x < 0){
    u <- 0
  } else{
    u <- x
  }
  
  return(u)
}


allocate_constrain <- function(x,radius){
  
  if (x > radius){
    u <- radius
  } else if (x < 0){
    u <- 0
  } else{
    u <- x
  }
  
  return(u)
}


# Function to transfer degree into radian ----

degree_to_radian <- function(d){
  
  r <- d * pi/180
  
  return(r)
  
}

# Set path ----
phi_gt_thres <- 38 # 10:40

orig_path <- 'D:/SZU/research/projects/egoism/experiment/data/R/update_2/learning_phase'
path <- file.path(orig_path, paste0('phi_gt_thres_',phi_gt_thres), 'diff_lamdas')

# Check if the directory exists; if not, create it
if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE)
}

# Create subdirectories
lik_path <- file.path(path, 'lik_all')
allocate_path <- file.path(path, 'allocate_offer_prescreen')
fb_path <- file.path(path, 'FB_prescreen')
phi_gt_path <- file.path(path, 'phi_gt_prescreen')

if (!dir.exists(lik_path)) {
  dir.create(lik_path, recursive = TRUE)
}

if (!dir.exists(allocate_path)) {
  dir.create(allocate_path, recursive = TRUE)
}

if (!dir.exists(fb_path)) {
  dir.create(fb_path, recursive = TRUE)
}

if (!dir.exists(phi_gt_path)) {
  dir.create(phi_gt_path, recursive = TRUE)
}

# Step 1: Pre-screen process ----
## 1. Generate all task parameters ----
sim_gn_n <- 100
dn <- 4
subj_num <- 30
bias <- 2
target_num <- 15
tn <- dn * target_num

# dc controls egoistic feedback
dc_set <- 5

sim_time <- 200
tp_num <- 300

# Range of task parameters
# lamda controls range of allocation offer phi+lamda:38 ~ 68
lamda_rng <- seq(0, 30, 5)

library(gtools)

# list all the combinations of four lamdas for four allocation plans
lamda_seq <- permutations(n = length(lamda_rng), r = dn, v = lamda_rng)
tp_all <- lamda_seq

colnames(tp_all) <- paste0("lamda_", 1:4)

# Select tp with different lamdas
keep <- apply(tp_all, 1, function(x){
  
  diffs <- abs(outer(x, x, FUN = "-"))
  diffs <- diffs[!diag(TRUE, length(x))]
  all(diffs >= 5)
  
})

tp <- as.data.frame(tp_all[keep, ])

fn_tp <- file.path(path, 'tp.txt')
write.table(tp, fn_tp, row.names = F, col.names = T)


# !!!!!!!!!!! for test
# tp <- tp[1:5,]

## 2. Simulate learning task and estimates ----
numCores <- detectCores()
registerDoParallel(numCores)


results <- foreach (st = 1:sim_time, .combine = rbind) %dopar% {
  
  # mn is indicator for all combination of task parameter
  # max (mn) = nrow(tp)[dc_range x lamda_range; tp_num] x number of simulated gn
  mn <- 1
  
  # Store likelihood (averaged across 30 subjects) for all combinations of task parameters
  # for each simulation
  Loglik_all <- matrix(0, nrow(tp) * sim_gn_n, 7)
  
  # Store stimulus data
  stimulus_data_list <- list()
  
  for (n in 1:dim(tp)[1]){
    
    phi_gt <- matrix(0, sim_gn_n, target_num)
    
    for (i in 1:sim_gn_n){
      
      ### Generate allocation offer ----
      allocate_offer <- matrix(0, target_num, dn)
      allocate_offer_new <- matrix(0, target_num, dn)
      
      for (target in 1:target_num){
        
        set.seed(1000 * i + target)
        phi_gt[i, target] <- sample(0:phi_gt_thres, size = 1, replace = T)
        
        for (tau in 1:dn){
          
          allocate_offer[target, tau] <- sapply(phi_gt[i, target] + rpois(1, lambda = tp[n,tau]), est_constrain)
          
        }
        
        # reorder allocation offer
        new_indx <- sample(seq(1,dn,1), dn, replace = F)
        allocate_offer_new[target, ] <- allocate_offer[target, new_indx]
        
      }
      
      ### Simulate feedback ----
      
      FB <- matrix(0, target_num, dn)
      
      # Generate dc for each set of simulated group norms for all targets
      set.seed(i)
      dc <- rpois(1, lambda = dc_set)
      
      for (target in 1:target_num){
        
        set.seed(10 * i + target)
        
        for (tau in 1:dn){
          
          if (allocate_offer_new[target, tau] <= phi_gt[i, target] + dc){
            FB[target, tau] <- 1
          } else {
            FB[target, tau] <- 0
          }
        }
    
      }
      
      # Store stimulus data
      stimulus_data_list[[mn]] <- list(
        lamda_1 = tp$lamda_1[n],
        lamda_2 = tp$lamda_2[n],
        lamda_3 = tp$lamda_3[n],
        lamda_4 = tp$lamda_4[n],
        gn = i,
        dc = dc,
        phi_gt = phi_gt[i, ],
        allocate_offer_new = allocate_offer_new,
        FB = FB
      )

      ### Simulate behaviors ----
      
      phi_gn <- matrix(0, 1, target_num)
      
      param_est <- matrix(0, subj_num, 4)
      param_gt <- matrix(0, subj_num, 3)
      Loglik_sum <- matrix(0, 1, subj_num)
      
      for (subj in 1:subj_num){
        
        seed <- 100 * st + subj
        set.seed(seed)
        alpha_gn <- runif(1, min = 0, max = 1)
        beta_gn <- runif(1, min = 0, max = 1)
        beta_0 <- sample(seq(-bias,bias,1), 1, replace = T)
        
        param_gt[subj, 1] <- alpha_gn
        param_gt[subj, 2] <- beta_gn
        param_gt[subj, 3] <- beta_0
        
        behav_est <- matrix(0, 1, target_num)
        
        for (target in 1:target_num) {
          
          phi_est <- matrix(0, 1, dn + 1)
          phi_est[1] <- 60
          
          for (tau in 1:dn){
            
            phi_est[tau+1] <- phi_est[tau] + alpha_gn * FB[target, tau] * (allocate_offer_new[target, tau] - phi_est[tau])
            
          }
          
          phi_gn[target] <- phi_est[dn+1]
          behav_est[target] <- round(est_constrain(beta_gn * phi_gn[target] + beta_0), 2) # phi_gn[target] <- phi_gn'
        }
        
        
        ### likelihood of ground truth ----
        
        phi_gn_mf <- matrix(0, 1, target_num)
        Loglik <- matrix(0, 1, target_num)

        
        for (target in 1:target_num){
          
          phi_est_mf <- matrix(0, 1, 4)
          phi_est_mf[1] <- 60
          
          for (tau in 1:dn){
            
            phi_est_mf[tau+1] <- phi_est_mf[tau] + alpha_gn * FB[target, tau] * (allocate_offer_new[target, tau] - phi_est_mf[tau])
            
          }
          
          phi_gn_mf[target] <- phi_est_mf[tau+1]
          
          u <- behav_est[target] - beta_gn * phi_gn_mf[target] - beta_0
          Loglik[target] <- log10(dnorm(u))
          
        }
        
        Loglik_sum[subj] <- sum(Loglik)
        
      }
      
      Loglik_all[mn,1] <- tp$lamda_1[n]
      Loglik_all[mn,2] <- tp$lamda_2[n]
      Loglik_all[mn,3] <- tp$lamda_3[n]
      Loglik_all[mn,4] <- tp$lamda_4[n]
      Loglik_all[mn,5] <- i
      Loglik_all[mn,6] <- dc
      Loglik_all[mn,7] <- sum(Loglik_sum)
      
      # mn iterates for all combinations of task parameters
      mn <- mn + 1
      
    }
  }
  
  # Loglik_all stores likelihood (summed over30 subjects) for all combinations of task parameters (nrow(tp)*sim_gn_n)
  # for each simulation; number of row: all combination of tp
  
  colnames(Loglik_all) <- c(paste0("lamda_",1:4), "gn", "dc", "ll")
  Loglik_all <- as.data.frame(Loglik_all)
  
  # Normalize log likelihood across all tp for each simulation
  z_ll_st <- (Loglik_all$ll - mean(Loglik_all$ll))/sd(Loglik_all$ll)

  # Dimension of results: sim_time (differs in subjects) x four stored results
  list(ll = Loglik_all,
       z_ll_st = z_ll_st, 
       st = st, 
       bias = bias,
       stimulus_data_list = stimulus_data_list)
  
}

# Write the results outside the foreach loop
cat("Saving results...\n")

# results: sim_time x four stored variables
for (i in 1:nrow(results)) {
  
  # Save z_ll result of each simulation
  fn_z_ll_st <- paste0(lik_path, '/z_ll_with_bias_',
                       results[i,]$bias, '_st_', results[i,]$st,'.txt')
  write.table(results[i,]$z_ll_st, fn_z_ll_st, row.names = F)
  
  cat(sprintf("Saved st %d z_ll\n", results[i,]$st))
}


# results: sim_time x four stored variables
# stimulus_data, length: nrow(tp) x sim_gn
cat("Saving stimulus data...\n")
stimulus_data <- results[1,]$stimulus_data_list
dc_all <- matrix(0, length(stimulus_data), 1)

# dc is the same for each set of lamdas
# dc_all, length: nrow(tp)(repeat) x sim_gn (seq)
for (j in 1:length(stimulus_data)){
  dc_all[j] <- stimulus_data[[j]]$dc
}

fn_dc <- file.path(path, 'dc_all.txt')
write.table(dc_all, fn_dc, col.names = F, row.names = F)


for (j in 1:length(stimulus_data)) {
  
  lamda_1 <- stimulus_data[[j]]$lamda_1
  lamda_2 <- stimulus_data[[j]]$lamda_2
  lamda_3 <- stimulus_data[[j]]$lamda_3
  lamda_4 <- stimulus_data[[j]]$lamda_4
  gn <- stimulus_data[[j]]$gn
  dc <- stimulus_data[[j]]$dc
  
  # Save phi_gt
  fn_phi_gt <- paste0(phi_gt_path, '/phi_gt_lamdas_', lamda_1, '_', lamda_2, '_', lamda_3, '_', lamda_4,
                      '_dc_', dc, '_gn_', gn, '.txt')
  write.table(stimulus_data[[j]]$phi_gt, fn_phi_gt, row.names = F, col.names = F)
  
  # Save allocate_offer_new
  fn_allocate <- paste0(allocate_path, '/allocate_offer_lamdas_', lamda_1, '_', lamda_2, '_', lamda_3, '_', lamda_4, 
                        '_dc_', dc, '_gn_', gn, '.txt')
  write.table(stimulus_data[[j]]$allocate_offer_new, fn_allocate, row.names = F, col.names = F)
  
  # Save FB
  fn_fb <- paste0(fb_path, '/FB_lamdas_', lamda_1, '_', lamda_2, '_', lamda_3, '_', lamda_4, 
                  '_dc_', dc, '_gn_', gn, '.txt')
  write.table(stimulus_data[[j]]$FB, fn_fb, row.names = F, col.names = F)
  
  if (j %% 1000 == 0) {
    cat(sprintf("Saved %d/%d stimulus files\n", j, length(stimulus_data)))
  }
}

cat("Simulation complete!\n")


# ==============================================================================
# Step 2: Read data and analyze ----
# ==============================================================================

rm(list=ls())

library(dplyr)

## Helper function to expand matrix ----

expand_matrix <- function(orig_matrix, n, mode){
  
  if (mode == 1){
    expand_matrix <- orig_matrix[rep(seq_len(nrow(orig_matrix)),each = n),]
  } else if (mode == 2){
    expand_matrix <- orig_matrix[rep(seq_len(nrow(orig_matrix)),n),]
  }
  
}

## 1. Set parameters ----

sim_gn_n <- 100
dn <- 4  # Changed back to 4 to match pre-screen
subj_num <- 30
bias <- 0
target_num <- 15
tn <- dn * target_num

sim_time <- 200

tp_num <- 300


## 2. Set paths ----

phi_gt_thres <- 38 # 10:40

orig_path <- 'D:/SZU/research/projects/egoism/experiment/data/R/update_2/learning_phase'
path <- file.path(orig_path, paste0('phi_gt_thres_',phi_gt_thres), 'diff_lamdas')

lik_path <- file.path(path, 'lik_all')
allocate_path <- file.path(path, 'allocate_offer_prescreen')
fb_path <- file.path(path, 'FB_prescreen')
phi_gt_path <- file.path(path, 'phi_gt_prescreen')

## 3. Read tp data ----
fn_tp <- file.path(path, 'tp.txt')
tp <- as.data.frame(read.table(fn_tp, header = T))

fn_dc <- file.path(path, 'dc_all.txt')
dc_all <- read.table(fn_dc, header = F)

# !!!!!!!!!! for test
# tp <- tp[1:5,]

## 4. Read allocation offer from saved allocation offer files ----

mn <- 1
allocate_offer_new_all <- matrix(0, nrow(tp) * sim_gn_n, target_num * dn+1)
prosocial_allocate_offer <- matrix(0, nrow(tp) * sim_gn_n, target_num+1)

cat("Reading allocation offer files...\n")

for (n in 1:dim(tp)[1]){
  
  for (i in 1:sim_gn_n){
    
    lamda_1 <- tp$lamda_1[n]
    lamda_2 <- tp$lamda_2[n]
    lamda_3 <- tp$lamda_3[n]
    lamda_4 <- tp$lamda_4[n]
    gn <- i
    
    # dc is the same for each set of lamdas
    # dc_all, length: nrow(tp)(repeat) x sim_gn (seq)
    dc <- dc_all[i+(n-1)*sim_gn_n,]
    
    # Read saved FB
    fn_allocate <- paste0(allocate_path, '/allocate_offer_lamdas_', lamda_1, '_', lamda_2, '_', lamda_3, '_', lamda_4, 
                          '_dc_', dc, '_gn_', gn, '.txt')
    allocate_offer_new <- as.matrix(read.table(fn_allocate))
    
    allocate_offer_new_all[mn,1] <- mn
    allocate_offer_new_all[mn,2:ncol(allocate_offer_new_all)] <- c(t(allocate_offer_new))
    
    # Count zeros in FB for each target
    for (target in 1:target_num){
      prosocial_allocate_offer[mn, target+1] <- sum(allocate_offer_new[target,]>=60)
    }
    
    # fb_num_zeros: nrow(tp)(rep) x sim_gn(seq)
    prosocial_allocate_offer[mn, 1] <- mn
    
    mn <- mn + 1
  }
  
  if (n %% 100 == 0) {
    cat(sprintf("Processed %d/%d parameter combinations\n", n, dim(tp)[1]))
  }
}

colnames(prosocial_allocate_offer) <- c("ID", paste0("prosocial_allocate_offer_", 1:target_num))
colnames(allocate_offer_new_all) <- c("ID", paste0("allocate_offer_", 1:(target_num * dn)))

## 5. Calculate frequency of rejection from saved FB files ----

mn <- 1
fb_num_zeros <- matrix(0, nrow(tp) * sim_gn_n, target_num+1)
phi_gt_all <- matrix(0, nrow(tp) * sim_gn_n, target_num+1)

cat("Reading FB and phi_gt files...\n")

for (n in 1:dim(tp)[1]){
  
  for (i in 1:sim_gn_n){
    
    lamda_1 <- tp$lamda_1[n]
    lamda_2 <- tp$lamda_2[n]
    lamda_3 <- tp$lamda_3[n]
    lamda_4 <- tp$lamda_4[n]
    gn <- i
    
    # dc is the same for each set of lamdas
    # dc_all, length: nrow(tp)(repeat) x sim_gn (seq)
    dc <- dc_all[i+(n-1)*sim_gn_n,]
    
    # Read saved FB
    fn_fb <- paste0(fb_path, '/FB_lamdas_', lamda_1, '_', lamda_2, '_', lamda_3, '_', lamda_4, 
                    '_dc_', dc, '_gn_', gn, '.txt')
    FB <- as.matrix(read.table(fn_fb))
    
    # Read saved phi_gt
    fn_phi_gt <- paste0(phi_gt_path, '/phi_gt_lamdas_', lamda_1, '_', lamda_2, '_', lamda_3, '_', lamda_4,
                        '_dc_', dc, '_gn_', gn, '.txt')
    phi_gt <- as.vector(as.matrix(read.table(fn_phi_gt)))
    
    # Count zeros in FB for each target
    for (target in 1:target_num){
      fb_num_zeros[mn, target+1] <- sum(FB[target, ] == 0)
    }
    
    # fb_num_zeros: nrow(tp)(rep) x sim_gn(seq)
    fb_num_zeros[mn, 1] <- mn
    phi_gt_all[mn, 1] <- gn
    phi_gt_all[mn, 2:(target_num+1)] <- phi_gt
    
    mn <- mn + 1
  }
  
  if (n %% 100 == 0) {
    cat(sprintf("Processed %d/%d parameter combinations\n", n, dim(tp)[1]))
  }
}

## 6. Create tp_all dataframe ----

tp_exp <- expand_matrix(tp, n = sim_gn_n, mode = 1)
tp_all <- cbind(1:(dim(tp)[1]*sim_gn_n), tp_exp, phi_gt_all, dc_all)

colnames(tp_all) <- c("ID", "lamda_1", "lamda_2", "lamda_3", "lamda_4",
                      "gn", paste0("phi_gt_", 1:target_num), "dc")
colnames(fb_num_zeros) <- c("ID", paste0("fb_zero_gn_", 1:target_num))

### Read z_ll files ----

cat("Reading z_ll files...\n")

z_ll <- matrix(0, nrow(tp) * sim_gn_n, sim_time)

for (i in 1:sim_time){
  
  fn <- paste0(lik_path, '/z_ll_with_bias_', bias, '_st_', i, '.txt')
  
  # z_ll: [nrow(tp)(rep) x sim_gn(seq)] x sim_time
  z_ll[,i] <- as.matrix(read.table(fn, header = T))
  
  if (i %% 50 == 0) {
    cat(sprintf("Read %d/%d z_ll files\n", i, sim_time))
  }
}

z_ll <- as.data.frame(z_ll)
z_ll$ID <- 1:(dim(tp)[1]*sim_gn_n)

# z_ll represents the ranking of each tp combination
# z_ll_m represents the average ranking of each tp combination of all simulation
z_ll$z_ll_m <- rowMeans(z_ll[, 1:sim_time])

### Merge all data ----

tp_all <- as.data.frame(tp_all)
z_ll <- as.data.frame(z_ll)
fb_num_zeros <- as.data.frame(fb_num_zeros)

data_all <- merge(merge(merge(merge(tp_all, z_ll, by = "ID"), 
                        fb_num_zeros, by = "ID"), 
                  prosocial_allocate_offer, by = "ID"), 
                  allocate_offer_new_all, by = "ID")

data <- data_all[order(data_all$z_ll_m, decreasing = T),]

data_new <- data %>%
  select("ID", starts_with("lamda_"), "gn", "dc",
         starts_with("phi_gt_"), 
         starts_with("fb_zero_gn_"),
         starts_with("prosocial_allocate_offer_"), 
         starts_with("allocate_offer_"), "z_ll_m") %>%
  mutate(sum_fb_zero_gn = rowSums(select(., starts_with("fb_zero_gn_")), na.rm = TRUE))

nrow(data_new)

fn_data <- paste0(path, '/lik_all_with_bias_',bias,'_data_new.csv')
write.csv(x = data_new, file = fn_data, row.names = F, fileEncoding = "UTF-8")

## 6. Search for task parameter ----

cat("Filtering results...\n")

data_new <- read.csv(file = fn_data, header = T, sep = ",", 
                     stringsAsFactors = F, na.strings = "", 
                     fileEncoding = "UTF-8")

thresh <- quantile(data_new$z_ll_m, probs = 0.85)

result <- data_new %>%
  filter(z_ll_m <= thresh) %>%
  # More than 1 reject for each demonstrator
  filter(if_all(starts_with("fb_zero_gn_"), ~ .x >= 1 & .x <= 2)) %>%
  # filter(if_all(starts_with("prosocial_allocate_offer_"), ~ .x >= 1 & .x <= 2)) %>%
  filter(sum_fb_zero_gn >= dn * target_num * 0.1 & sum_fb_zero_gn <= dn * target_num * 0.5)

cat(sprintf("Found %d qualifying parameter combinations\n", nrow(result)))

# Save filtered results ----

fn_lik_thres<- paste0(path, '/lik_all_with_bias_',bias,'_thres_0.85_zero_fb_1~0.5_z_ll_m.txt')
write.table(result, fn_lik_thres, row.names=F)

cat("Saved filtered results\n")


# ==============================================================================
# Step 3: Save stimulus for selected parameters (directly from saved files) ----
# ==============================================================================

# Function to constrain behavior within range ----

est_constrain <- function(x){
  
  if (x > 90){
    u <- 90
  } else if (x < 0){
    u <- 0
  } else{
    u <- x
  }
  
  return(u)
}

# Create output directories ----

output_path <- file.path(path, 'z_ll_m')
output_fb_path <- file.path(output_path, 'FB')
output_allocate_path <- file.path(output_path, 'allocate_offer')
output_behav_path <- file.path(output_path, 'behav_est')
output_param_path <- file.path(output_path, 'gt_param')

if (!dir.exists(output_fb_path)) {
  dir.create(output_fb_path, recursive = TRUE)
}

if (!dir.exists(output_allocate_path)) {
  dir.create(output_allocate_path, recursive = TRUE)
}

if (!dir.exists(output_behav_path)) {
  dir.create(output_behav_path, recursive = TRUE)
}

if (!dir.exists(output_param_path)) {
  dir.create(output_param_path, recursive = TRUE)
}

# Process selected parameters ----

cat("Processing selected parameters...\n")

dn <- 4  # Match pre-screen dn
param_gt <- matrix(0, subj_num, 3)

for (j in 1:dim(result)[1]){
  
  lamda_1 <- result$lamda_1[j]
  lamda_2 <- result$lamda_2[j]
  lamda_3 <- result$lamda_3[j]
  lamda_4 <- result$lamda_4[j]
  gn <- result$gn[j]
  dc <- result$dc[j]
  
  # Read saved allocate_offer
  fn_allocate_saved <- paste0(allocate_path, '/allocate_offer_lamdas_', lamda_1, '_', lamda_2, '_', lamda_3, '_', lamda_4,
                              '_dc_', dc, '_gn_', gn, '.txt')
  allocate_offer_new <- as.matrix(read.table(fn_allocate_saved))
  
  # Read saved FB
  fn_fb_saved <- paste0(fb_path, '/FB_lamdas_', lamda_1, '_', lamda_2, '_', lamda_3, '_', lamda_4, 
                        '_dc_', dc, '_gn_', gn, '.txt')
  FB <- as.matrix(read.table(fn_fb_saved))
  
  # Read saved phi_gt
  fn_phi_gt_saved <- paste0(phi_gt_path, '/phi_gt_lamdas_', lamda_1, '_', lamda_2, '_', lamda_3, '_', lamda_4,
                            '_dc_', dc, '_gn_', gn, '.txt')
  phi_gn <- as.vector(as.matrix(read.table(fn_phi_gt_saved)))
  
  cat(sprintf("Processing %d/%d: lamda_1=%d, lamda_2=%d,lamda_3=%d,lamda_4=%d,
              gn=%d, dc=%f\n",j, dim(result)[1], lamda_1, lamda_2, lamda_3, lamda_4, gn, dc))
  
  ## Simulate behaviors using saved data ----
  
  behav_est_all <- matrix(0, subj_num, target_num)
  
  for (subj in 1:subj_num){
    
    set.seed(subj)  # Fixed seed for reproducibility
    alpha_gn <- runif(1, min = 0, max = 1)
    beta_gn <- runif(1, min = 0, max = 1)
    beta_0 <- sample(seq(-bias,bias,1), 1, replace = T)
    
    param_gt[subj, 1] <- alpha_gn
    param_gt[subj, 2] <- beta_gn
    param_gt[subj, 3] <- beta_0
    
    behav_est <- matrix(0, 1, target_num)
    phi_gn_est <- matrix(0, 1, target_num)
    
    for (target in 1:target_num) {
      
      phi_est <- matrix(0, 1, dn + 1)
      phi_est[1] <- 60
      
      for (tau in 1:dn){
        
        phi_est[tau+1] <- phi_est[tau] + alpha_gn * FB[target, tau] * (allocate_offer_new[target, tau] - phi_est[tau])
        
      }
      
      phi_gn_est[target] <- phi_est[dn+1]
      behav_est[target] <- round(est_constrain(beta_gn * phi_gn_est[target] + beta_0), 2)
    }
    
    behav_est_all[subj, ] <- behav_est
  }
  
  # Save to output directory
  fn_fb_out <- paste0(output_fb_path, '/lamdas_', lamda_1, '_', lamda_2, '_', lamda_3, '_', lamda_4,
                      '_dc_', dc, '_gn_', gn, '.txt')
  write.table(FB, fn_fb_out, row.names = F, col.names = F)
  
  fn_offer_out <- paste0(output_allocate_path, '/lamdas_', lamda_1, '_', lamda_2, '_', lamda_3, '_', lamda_4, 
                         '_dc_', dc, '_gn_', gn, '.txt')
  write.table(allocate_offer_new, fn_offer_out, row.names = F, col.names = F)
  
  fn_behav_est <- paste0(output_behav_path, '/lamdas_', lamda_1, '_', lamda_2, '_', lamda_3, '_', lamda_4,
                         '_dc_', dc, '_gn_', gn, '.txt')
  write.table(behav_est_all, fn_behav_est, row.names = F, col.names = F)
  
}

# Save ground truth parameters
fn_param_gt <- paste0(output_param_path, '/param_gt.txt')
write.table(param_gt, fn_param_gt, row.names = F, col.names = F)

cat("All processing complete!\n")
cat("Output saved in:", output_path, "\n")

