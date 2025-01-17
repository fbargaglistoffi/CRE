set.seed(2021)
library(causalTree)

# Set Experiment Parameter
n_rules <- 2
sample_size <- 5000
effect_size <- 2
confounding <- "nc"
ratio_dis <- 0.5

# Set Ground Truth
{
  if (n_rules==1) {
    dr <- c("x1>0.5 & x2<=0.5")
  } else if (n_rules==2) {
    dr <- c("x1>0.5 & x2<=0.5", "x5>0.5 & x6<=0.5")
  } else if (n_rules==3) {
    dr <- c("x1>0.5 & x2<=0.5", "x5>0.5 & x6<=0.5", "x4>0.5")
  } else if (n_rules==4) {
    dr <- c("x1>0.5 & x2<=0.5", "x5>0.5 & x6<=0.5",
            "x4>0.5", "x5<=0.5 & x7>0.5 & x8<=0.5")
  } else {
    stop(paste("Synthtic dataset with", n_rules,"rules has not been
                implemented yet. Available 'n_rules' options: {1,2,3,4}."))
  }
}

# Generate Dataset
dataset <- generate_cre_dataset(n = sample_size,
                                rho = 0,
                                p = 10,
                                effect_size = effect_size,
                                n_rules = n_rules,
                                binary_covariates = TRUE,
                                binary_outcome = FALSE,
                                confounding = confounding)
y <- dataset[["y"]]
z <- dataset[["z"]]
X <- dataset[["X"]]
ite <- dataset[["ite"]]
X_names <- colnames(X)

# Preprocessing
subgroups <- honest_splitting(y, z, X, ratio_dis)
discovery <- subgroups[["discovery"]]
inference <- subgroups[["inference"]]
y_dis <- discovery$y
z_dis <- discovery$z
X_dis <- discovery$X
y_inf <- inference$y
z_inf <- inference$z
X_inf <- inference$X
data_dis <- as.data.frame(cbind(y_dis, X_dis))
data_inf <- as.data.frame(cbind(y_inf, z_inf, X_inf))

# Honest Causal Tree Model
fit.tree <- causalTree(y_dis ~ ., data = data_dis, treatment = z_dis,
                       split.Rule = "CT", cv.option = "CT",
                       split.Honest = T, cv.Honest = T, maxdepth = 3)
opt.cp <- fit.tree$cptable[,1][which.min(fit.tree$cptable[,4])]
pruned <- prune(fit.tree, opt.cp)
# Visualize
rpart.plot(pruned)

# Extract Causal Decision Rules
rules <- as.numeric(row.names(pruned$frame[pruned$numresp]))
rules.ctree <- vector("list",length(rules))
for (k in rules[-1]){
  sub <- as.data.frame(matrix(NA, nrow = 1,
                              ncol = nrow(as.data.frame(path.rpart(pruned, node=k, print.it = FALSE)))-1))
  capture.output(for (h in 1:ncol(sub)){
    sub[,h] <- as.character(print(as.data.frame(path.rpart(pruned,node=k,print.it=FALSE))[h+1,1]))
    sub_pop <- noquote(paste(sub , collapse = " & "))
  })
  subset <- with(data_inf, data_inf[which(eval(parse(text=sub_pop))),])
  if (length(unique(subset$z_inf))!= 1){
    rules.ctree[[k]] <- sub_pop
  }
}
leaves_dfs <- rep(FALSE,length(rules))
leaves_dfs[unique(pruned$where)] <- TRUE
cdr_pred <- unlist(rules.ctree[rules[leaves_dfs]])
for (i in 1:length(cdr_pred)){
  cdr_pred[i] <- gsub("< ", "<=", cdr_pred[i])
  cdr_pred[i] <- gsub(">=", ">", cdr_pred[i])
}
# Predict
ite_pred <- predict(pruned, X)

metrics_cdr <- evaluate(cdr,cdr_pred)
print(paste("Causal Decision Rules: ",
            "IoU=",round(metrics_cdr$IoU,2),
            ", Recall=",round(metrics_cdr$recall,2),
            ", Precision=",round(metrics_cdr$precision,2),
            sep=""))

em_pred <- extract_effect_modifiers(cdr_pred, X_names)
metrics_em <- evaluate(em,em_pred)
print(paste("Effect Modifiers:      ",
            "IoU=",round(metrics_em$IoU,2),
            ", Recall=",round(metrics_em$recall,2),
            ", Precision=",round(metrics_em$precision,2),
            sep=""))

#Estimation
rmse <- sqrt(mean((ite - ite_pred)^2))
print(paste("RMSE: ", rmse))
bias <- mean((ite - ite_pred))
print(paste("Bias: ", bias))
