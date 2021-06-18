library(data.table)
library(readr)


# hmm is 3 hidden satates and 3 agent responses
hmm = initHMM(1:3, 1:9, nAgentResponse = 3,
              startProbs = c(1/3,1/3,1/3),
              transProbs = matrix(c(0.9,0.7,0.5,0.4,0.05,0.01,0.1,0.05,0.01,0.09,0.2,0.4,0.5,0.9,0.4,0.5,0.15,0.05,0.01,0.1,0.1,0.1,0.05,0.59,0.4,0.8,0.94),9),
              emissionProbs = matrix(c(0.3,1/60,1/300,0.3,1/60,1/300,0.3,1/60,1/300,0.03,0.3,1/60,0.03,0.3,1/60,0.03,0.3,1/60,1/300,1/60,47/150,1/300,1/60,47/150,1/300,1/60,47/150),3))
# hmm2 is 2 hidden state and 2 agent responses
hmm2 = initHMM(1:2, 1:6, nAgentResponse = 2,
               startProbs = c(1/2,1/2),
               transProbs = matrix(c(0.8,0.6,0.2,0.03,0.2,0.4,0.8,0.97),4),
               emissionProbs = matrix(c(0.25,0.05,0.25,0.05,3/16,0.15,3/16,0.15,1/16,0.3,1/16,0.3),2))
setting_param_list = list(seed=1:20,n_sessions=c(10,50,100,200,500,1000),start_point=1:1)
other_param_list = list(original_hmm=hmm)
res = run_sim(runTest,setting_param_list,other_param_list,nworkers = 55)

find_mean <- function(some_list) {
  reduce(some_list, `+`) / length(some_list)
}

graph_dots_trans_mat = c()
graph_dots_emiss_mat = c()
graph_dots_start_probs = c()
for (n_sessions in setting_param_list$n_sessions) {
  # divide result according to column val
  res_sessions <- res[res$n_sessions==n_sessions,]
  # get list of objects
  res_sessions <- res_sessions %>% pull(res)
  # ignore "differences" vector (every second variable in the original list)
  res_sessions[seq(0,length(res_sessions),2)] <- NULL
  
  # get list of transition matrices
  mat1_list <- res_sessions %>% map(~.$transProbs)
  # get list of emission matrices
  mat2_list <- res_sessions %>% map(~.$emissionProbs)
  # get list of startProbs
  vec_list <- res_sessions %>% map(~.$startProbs)
  
  mat1_diff_vec = c()
  for (mat1 in mat1_list) {
    mat1_diff_vec = c(mat1_diff_vec,sum((mat1-hmm$transProbs)^2))
  }
  mat2_diff_vec = c()
  for (mat2 in mat2_list) {
    mat2_diff_vec = c(mat2_diff_vec,sum((mat2-hmm$emissionProbs)^2))
  }
  vec_diff_vec = c()
  for (vec in vec_list) {
    vec_diff_vec = c(vec_diff_vec,sum((vec-hmm$startProbs)^2))
  }
  
  fwrite(list(mat1_diff_vec),file=paste(c("exp1_trans_mat_boxplot_", n_sessions, ".txt"), collapse = ""))
  fwrite(list(mat2_diff_vec),file=paste(c("exp1_emiss_mat_boxplot_", n_sessions, ".txt"), collapse = ""))
  fwrite(list(vec_diff_vec),file=paste(c("exp1_start_prob_boxplot_", n_sessions, ".txt"), collapse = ""))
  
  
  mat1_mean <- find_mean(mat1_list)
  mat2_mean <- find_mean(mat2_list)
  vec_mean <- find_mean(vec_list)
  
  dif_trans = sum((mat1_mean-hmm$transProbs)^2)
  dif_emis = sum((mat2_mean-hmm$emissionProbs)^2)
  dif_start = sum((vec_mean-hmm$startProbs)^2)
  
  graph_dots_trans_mat= c(graph_dots_trans_mat,dif_trans)
  graph_dots_emiss_mat= c(graph_dots_emiss_mat,dif_emis)
  graph_dots_start_probs= c(graph_dots_start_probs,dif_start)
  
}

fwrite(list(graph_dots_trans_mat),file="exp1_trans_mat_boxplot.txt")
fwrite(list(graph_dots_emiss_mat),file="exp1_emiss_mat_boxplot.txt")
fwrite(list(graph_dots_start_probs),file="exp1_start_probs_boxplot.txt")


png(file = "exp1_results_trans_mat.jpg")
plot(c(10,50,100,200,500,1000), graph_dots_trans_mat,type = "o",xlab = "number of sessions",ylab = "sum of squared error")
dev.off()
png(file = "exp1_results_emiss_mat.jpg")
plot(c(10,50,100,200,500,1000), graph_dots_emiss_mat,type = "o",xlab = "number of sessions",ylab = "sum of squared error")
dev.off()
png(file = "exp1_results_startprob_vec.jpg")
plot(c(10,50,100,200,500,1000), graph_dots_start_probs,type = "o",xlab = "number of sessions",ylab = "sum of squared error")
dev.off()