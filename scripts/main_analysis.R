# =========================================================================
# NOTE: THIS SCRIPT IS FOR REFERENCE ONLY.
# It details the initial cleaning and merging of the raw dunnhumby
# 'Carbo-Loading' data.
# It cannot be run without the original data files, which are not shared
# due to licensing restrictions.
# =========================================================================


library(dplyr)
library(Matrix)
library(igraph)

FILE_PATH <- "data_intermediate/replication_edges_clean.csv"

ALPHA_M <- 0.01
ALPHA_L <- 0.20

# STEP 1: DATA PREPERATION CREATING BIPARTITE AND BIADJACNEY MATRIX 
df <- read.csv(FILE_PATH, na.strings = c("", "NA")) %>%
  filter(!is.na(order_id), !is.na(product_name), !is.na(product_name)) %>%
  mutate(
    order_cat   = as.integer(factor(order_id)),
    product_cat = as.integer(factor(product_name))
  )

n_transactions <- n_distinct(df$order_id)
n_products     <- n_distinct(df$product_name)

A_b <- sparseMatrix(
  i = df$order_cat,
  j = df$product_cat,
  x = 1,
  dims = c(n_transactions, n_products)
)

products <- sort(unique(df$product_name))
colnames(A_b) <- products



# STEP 2: NULL MODEL — Bipartite Configuration Model (BiCM)

cn_matrix <- crossprod(A_b)   # observed co-purchases: t(A_b) %*% A_b
diag(cn_matrix) <- 0

d_p <- colSums(A_b)           # product degrees
d_t <- rowSums(A_b)           # basket sizes
m   <- sum(d_p)               # total edges

mean_dt  <- mean(d_t)
mean_dt2 <- mean(d_t^2)

# BiCM expected common neighbours μ_ij
expected_cn_bicm <- outer(d_p, d_p, function(dp_i, dp_j) {
  (dp_i * dp_j / m) * ((mean_dt2 - mean_dt) / mean_dt)
})


# For "significantly more": P(X >= cn_ij) < alpha_m, which is 1 - P(X <= cn_ij - 1) < alpha_m
A_m <- (1 - ppois(as.matrix(cn_matrix) - 1, lambda = expected_cn_bicm)) < ALPHA_M

# For "significantly less": P(X <= cn_ij) < alpha_l
A_l <- ppois(as.matrix(cn_matrix), lambda = expected_cn_bicm) < ALPHA_L

# Substitutes must also share ≥1 complement
has_common_complement <- (A_m %*% t(A_m)) > 0

A_c <- A_m * 1
A_s <- (A_l & has_common_complement) * 1
diag(A_s) <- 0



# STEP 3: ORIGINAL MEASURES FOR WEIGHTED NETWORKS 

#EQUATION 1 - ORIGINAL COMPLEMENTARITY MEASURE 
d_t <- rowSums(A_b)
d_t_inv <- ifelse(d_t > 0, 1/d_t, 0)
N <- t(A_b) %*% Diagonal(x = d_t_inv) %*% A_b

denom_vec <- t(A_b) %*% d_t_inv
norm_matrix_denom <- sqrt(outer(denom_vec, denom_vec, "*"))
norm_matrix_denom[norm_matrix_denom == 0] <- 1 

sim_o_matrix <- N / norm_matrix_denom

#Weighted adjancency matrix of complement unipartite network a_c *sim
W_c <- A_c * sim_o_matrix 


# EQUATION 2 - SUBSTITUABILITY WEIGHTED MEASURES 

# Compute row norms of W_c
row_norms <- sqrt(rowSums(W_c^2))
# Normalise each row of W_c
Wc_norm <- W_c / row_norms
# Row-wise cosine similarities = dot product of normalised rows
sim_s_matrix <- Wc_norm %*% t(Wc_norm)
# Weighted adjancency matrix of substitiute unipartite network  a_s * sim 
W_s <- sim_s_matrix * A_s

# Apply complementarity score threshold
non_zero_W_c <- W_c[W_c > 0]
threshold_c <- quantile(non_zero_W_c, probs = 0.35)
W_c[W_c < threshold_c] <- 0

# Apply substitutability score threshold
threshold_s <- 0
W_s[W_s < threshold_s] <- 0

# STEP 4: ROLE EXTRACTION USING INFOMAP

# Infomap requires edge list format
g_comp <- graph_from_adjacency_matrix(W_c, mode = "undirected", weighted = TRUE, diag = FALSE)
g_subs <- graph_from_adjacency_matrix(W_s, mode = "undirected", weighted = TRUE, diag = FALSE)

# Get the community assignments
comp_roles <- cluster_infomap(g_comp, e.weights = E(g_comp)$weight)
subs_roles <- cluster_infomap(g_subs, e.weights = E(g_subs)$weight)

print(comp_roles)
print(subs_roles)


#=====================================================================================================================================================
#COMPUTING THE NETWORK-LEVEL VARIABLES
#=====================================================================================================================================================

library(igraph)

# Convert adjacency matrices to igraph objects
g_comp <- graph_from_adjacency_matrix(W_c, mode = "undirected", weighted = TRUE, diag = FALSE)
g_subs <- graph_from_adjacency_matrix(W_s, mode = "undirected", weighted = TRUE, diag = FALSE)

# Compute centralities for complementarity network
node_comp <- strength(g_comp)                                     #NODE STRENGTH

original_weights_comp <- E(g_comp)$weight

# Stronger ties (higher weight) become shorter distances (lower value).
distance_weights_comp <- 1 / original_weights_comp

betw_comp <- betweenness(g_comp, 
                         directed = FALSE, 
                         weights = distance_weights_comp, 
                         normalized = TRUE)                                           # BETWEENNESS CENTRALITY 
clust_comp <- transitivity(g_comp, type = "local", isolates = "zero")                 #CLUSTERING COEFFICIENT

# Compute centralities for substitutability network
node_subs <- strength(g_subs)                                       #NODE STRENGTH 


original_weights_subs <- E(g_subs)$weight

distance_weights_subs <- 1 / original_weights_subs

betw_subs <- betweenness(g_subs, 
                         directed = FALSE, 
                         weights = distance_weights_subs, 
                         normalized = TRUE)                                         # BETWEENNESS CENTRALITY

clust_subs <- transitivity(g_subs, type = "local", isolates = "zero")               # CLUSTERING COEFFICIENT 

#  product-level centrality table
network_info <- data.frame(
  product = names(node_comp),
  node_comp = node_comp,
  betweenness_comp = betw_comp,
  clustering_comp = clust_comp,
  node_subs = node_subs,
  betweenness_subs = betw_subs,
  clustering_subs = clust_subs
)

# Save to file
write.csv(network_info, "network_variables_new.csv", row.names = FALSE)

#=====================================================================================================================================================
# EXPORTING THE PRODUCT PAIRS 
#=====================================================================================================================================================

# Ensure product names are stored
products <- colnames(W_c)

Wc_mat <- as.matrix(W_c)
Ws_mat <- as.matrix(W_s)


df_comp <- as.data.frame(as.table(Wc_mat)) %>%
  rename(product_i = Var1, product_j = Var2, complementarity_score = Freq)

df_subs <- as.data.frame(as.table(Ws_mat)) %>%
  rename(product_i = Var1, product_j = Var2, substitutability_score = Freq)

# Merge 
df_pairs <- df_comp %>%
  left_join(df_subs, by = c("product_i", "product_j"))


df_pairs <- df_pairs %>%
  filter(product_i != product_j)

# keep only nonzero scores
df_pairs <- df_pairs %>%
  filter(complementarity_score > 0 | substitutability_score > 0)

# Replace numeric indices with product names (if available)
if (!is.null(products)) {
  df_pairs$product_i <- products[as.integer(df_pairs$product_i)]
  df_pairs$product_j <- products[as.integer(df_pairs$product_j)]
}

# Deduplicate symmetric pairs (keep only one of i–j vs j–i)
df_pairs <- df_pairs %>%
  rowwise() %>%
  mutate(pair_id = paste(sort(c(product_i, product_j)), collapse = "_")) %>%
  ungroup() %>%
  group_by(pair_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(-pair_id)

# Inspect
head(df_pairs)

products <- sort(unique(df$product_name))
rownames(W_c) <- colnames(W_c) <- products
rownames(W_s) <- colnames(W_s) <- products

write.csv(df_pairs, "product_pairs_for_stage3.csv", row.names = FALSE)

pairs <- read.csv("product_pairs_for_stage2")





#=====================================================================================================================================================
# CALCULATING THE PRODUCT LEVEL VARIABLES
#=====================================================================================================================================================


library(readxl)
library(writexl)
library(dplyr)

pairs <- read.csv("product_pairs_for_stage3.csv", stringsAsFactors = FALSE)

merged <- readxl::read_excel("merged_master.xlsx")

#  BUILD PRODUCT-LEVEL TABLE

prod_info <- merged %>%
  group_by(product_description) %>%
  summarise(
    brand      = dplyr::first(brand),
    category   = dplyr::first(commodity),
    avg_price  = sum(dollar_sales, na.rm = TRUE) / sum(units, na.rm = TRUE),
    promo_rate = mean((coupon == 1) | !is.na(feature_desc) | !is.na(display_desc), na.rm = TRUE),
    .groups    = "drop"
  )

# STEP 3: MERGE INTO PAIRS 
pairs_attr <- pairs %>%
  # attach attributes for i
  left_join(prod_info, by = c("product_i" = "product_description")) %>%
  rename(
    brand_i    = brand,
    category_i = category,
    price_i    = avg_price,
    promo_i    = promo_rate
  ) %>%
  # attach attributes for j
  left_join(prod_info, by = c("product_j" = "product_description")) %>%
  rename(
    brand_j    = brand,
    category_j = category,
    price_j    = avg_price,
    promo_j    = promo_rate
  )

#  CREATE PAIR-LEVEL FEATURES 
pairs_attr <- pairs_attr %>%
  mutate(
    same_brand    = as.integer(brand_i    == brand_j),
    same_category = as.integer(category_i == category_j),
    average_price = (price_i + price_j) / 2, # <-- ADD THIS LINE
    price_difference    = abs(price_i - price_j),
    promo_overlap = as.integer(promo_i > 0 & promo_j > 0)
  )


writexl::write_xlsx(pairs_attr, "product_level_variables2.xlsx")


#=====================================================================================================================================================
# CALCULATING THE PRODUCT LEVEL VARIABLES + FINAL REG SET
#=====================================================================================================================================================


library(readxl)
library(writexl)
library(dplyr)


# STEP 1: LOAD FILES

pairs <- read_excel("product_level_variables2.xlsx")

network_info <- read.csv("network_variables_new.csv")


# merge info product_i

pairs_net <- pairs %>%
  left_join(network_info, by = c("product_i" = "product")) %>%
  rename_with(~ paste0(.x, "_i"), 
              c("node_comp","betweenness_comp","clustering_comp",
                "node_subs","betweenness_subs","clustering_subs"))


# merge info product_j

pairs_net <- pairs_net %>%
  left_join(network_info, by = c("product_j" = "product")) %>%
  rename_with(~ paste0(.x, "_j"), 
              c("node_comp","betweenness_comp","clustering_comp",
                "node_subs","betweenness_subs","clustering_subs"))


# pair level absolute differences

pairs_net <- pairs_net %>%
  mutate(
    strength_comp_diff = abs(node_comp_i - node_comp_j),
    betweenness_comp_diff = abs(betweenness_comp_i - betweenness_comp_j),
    clustering_comp_diff = abs(clustering_comp_i - clustering_comp_j),
    
    strength_subs_diff = abs(node_subs_i - node_subs_j),
    betweenness_subs_diff = abs(betweenness_subs_i - betweenness_subs_j),
    clustering_subs_diff = abs(clustering_subs_i - clustering_subs_j)
  )


# STEP 5: SAVE MASTER DATASET
write_xlsx(pairs_net, "dv_iv_variables2.xlsx")


#=====================================================================================================================================================
# ORINDARY LEAST SQUARES REGRESSION
#=====================================================================================================================================================

library(car)
library(stargazer)
library(dplyr)
library(readxl)


master_df <- read_excel("dv_iv_variables.xlsx") # Use the correct filename

regression_df <- master_df %>%
  dplyr::select(
    complementarity_score, substitutability_score,
    same_brand, same_category, average_price, price_difference, promo_overlap,
    strength_comp_diff, betweenness_comp_diff, clustering_comp_diff,
    strength_subs_diff, betweenness_subs_diff, clustering_subs_diff
  )

# COMP MODEL (regression analysis for complementarity score as a DV)

# Model 1: Product-level variables ONLY
comp_model_1 <- lm(complementarity_score ~ same_brand + same_category + average_price + price_difference + promo_overlap,
                   data = regression_df)

# Model 2: Network variables (Complement network ONLY)
comp_model_2 <- lm(complementarity_score ~ strength_comp_diff + betweenness_comp_diff + clustering_comp_diff,
                   data = regression_df)

# Model 3: Full model (Product + Complement network)
comp_model_3 <- lm(complementarity_score ~ same_brand + same_category + average_price + price_difference + promo_overlap +
                     strength_comp_diff + betweenness_comp_diff + clustering_comp_diff,
                   data = regression_df)

# SUBS MODEL (regression analysis for subsituability score as DV)

# Model 1: Product-level variables ONLY
subs_model_1 <- lm(substitutability_score ~ same_brand + same_category + average_price + price_difference + promo_overlap,
                   data = regression_df)

# Model 2: Network variables (Substitute network ONLY)
subs_model_2 <- lm(substitutability_score ~ strength_subs_diff + betweenness_subs_diff + clustering_subs_diff,
                   data = regression_df)

# Model 3: Full model (Product + Substitute network)
subs_model_3 <- lm(substitutability_score ~ same_brand + same_category + average_price + price_difference + promo_overlap +
                     strength_subs_diff + betweenness_subs_diff + clustering_subs_diff,
                   data = regression_df)


# Multicollineary checks VIF

cat("\n--- VIF CHECKS ---\n")
print(vif(comp_model_3))
print(vif(subs_model_3))

# Mean VIF helper
mean_vif <- function(model) mean(vif(model))
cat("\nMean VIF (Comp Model 3):", mean_vif(comp_model_3), "\n")
cat("Mean VIF (Subs Model 3):", mean_vif(subs_model_3), "\n")



# Complementarity models
stargazer(comp_model_1, comp_model_2, comp_model_3,
          type = "html",
          out = "3complementarity_regression_OLS.html",
          title = "Regression Results for Complementarity Score OLS", 
          dep.var.labels = "Complementarity Score",
          column.labels = c("Model 1: Product Only", "Model 2: Network Only", "Model 3: Full Model"),
          covariate.labels = c("Same Brand", "Same Category", "Average Price", "Price Difference", "Promotion Overlap",  
                               "Strength Comp Diff", "Betweenness Comp Diff", "Clustering Comp Diff"),
          standardized = TRUE, 
          no.space = TRUE,
          align = TRUE) 

# Substitutability models
stargazer(subs_model_1, subs_model_2, subs_model_3,
          type = "html",
          out = "3substitutability_regression_OLS.html",
          title = "Regression Results for Substitutability Score OLS", 
          dep.var.labels = "Substitutability Score",
          column.labels = c("Model 1: Product Only", "Model 2: Network Only", "Model 3: Full Model"),
          covariate.labels = c("Same Brand", "Same Category", "Average Price", "Price Difference", "Promotion Overlap",
                               "Strength Subs Diff", "Betweenness Subs Diff", "Clustering Subs Diff"),
          standardized = TRUE, 
          no.space = TRUE,
          align = TRUE) 


#=====================================================================================================================================================
# TOBIT REGRESSION
#=====================================================================================================================================================

if (!require("AER")) install.packages("AER")   # contains 'tobit'
library(AER)


tobit_comp <- tobit(
  complementarity_score ~ same_brand + same_category + average_price + price_difference+
    promo_overlap + strength_comp_diff + betweenness_comp_diff + clustering_comp_diff,
  left = 0, right = 1,
  data = regression_df
)

tobit_subs <- tobit(
  substitutability_score ~ same_brand + same_category + average_price + price_difference +
    promo_overlap + strength_subs_diff + betweenness_subs_diff + clustering_subs_diff,
  left = 0, right = 1,
  data = regression_df
)


# Define clean, professional labels for your variables
covariate_labels <- c("Same Brand", "Same Category", "Average Price", 
                      "Price Difference", "Promotion Overlap", 
                      "Degree Difference", "Betweenness Difference", 
                      "Clustering Difference")

# Generate the combined table
stargazer(tobit_comp, tobit_subs,   # List of models
          type = "html",             # Change to "latex" or "html" for your document
          title = "Tobit Regression Models for Complementarity and Substitutability",
          align = TRUE,
          no.space = TRUE,           # Removes extra space for a compact table
          dep.var.labels = c("Complementarity Score", "Substitutability Score"),
          covariate.labels = covariate_labels,
          out = "2tobit_models.html"   # Saves the output to a file
)








