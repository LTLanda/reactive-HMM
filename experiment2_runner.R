library(data.table)
library(readr)
library(tidyverse)


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
setting_param_list = list(seed=1:1,n_sessions=c(1000),start_point=1:100)
other_param_list = list(original_hmm=hmm)
res = run_sim(runTest,setting_param_list,other_param_list,nworkers = 3)
# divide result according to column val
res_1000_sessions <- res[res$n_sessions==1000,]
# get list of objects
res_1000_sessions <- res_1000_sessions %>% pull(res)
# ignore "differences" vector (every second variable in the original list)
res_1000_sessions[seq(0,length(res_1000_sessions),2)] <- NULL
# get list of transition matrices
mat1_list_res_1000_sessions <- res_1000_sessions %>% map(~.$transProbs)
# get list of emission matrices
mat2_list_res_1000_sessions <- res_1000_sessions %>% map(~.$emissionProbs)
# get list of startProbs
vec_list_res_1000_sessions <- res_1000_sessions %>% map(~.$startProbs)

# fill boxplot vectors
find_mean <- function(some_list) {
  reduce(some_list, `+`) / length(some_list)
}

get_row_perm = function(f,s,t){
  return (c(((f-1)*nAgentResponse+1):(f*nAgentResponse),(((s-1)*nAgentResponse+1):(s*nAgentResponse)),(((t-1)*nAgentResponse+1):(t*nAgentResponse))))
}

mat1_diff_vec = c()
m=Inf
best_p = -1
for (p in permn(1:3)) {
  r = get_row_perm(p[1],p[2],p[3])
  mat1_diff_vec = c()
  for (mat1 in mat1_list_res_1000_sessions) {
    mat1_diff_vec = c(mat1_diff_vec,sqrt(sum((mat1[r,p,drop=FALSE]-hmm$transProbs)^2)))
  }
  mat2_diff_vec = c()
  for (mat2 in mat2_list_res_1000_sessions) {
    mat2_diff_vec = c(mat2_diff_vec,sqrt(sum((mat2[p,,drop=FALSE]-hmm$emissionProbs)^2)))
  }
  vec_diff_vec = c()
  for (vec in vec_list_res_1000_sessions) {
    vec_diff_vec = c(vec_diff_vec,sqrt(sum((vec[p,drop=FALSE]-hmm$startProbs)^2)))
  }
  temp = mean(mat1_diff_vec)+mean(mat2_diff_vec)+mean(vec_diff_vec)
  if (temp<m) {
    best_p=p
  }
}
print(best_p)
r = get_row_perm(best_p[1],best_p[2],best_p[3])
mat1_diff_vec = c()
for (mat1 in mat1_list_res_1000_sessions[2:length(mat1_list_res_1000_sessions)]) {
  mat1_diff_vec = c(mat1_diff_vec,sqrt(sum((mat1[r,best_p,drop=FALSE]-mat1_list_res_1000_sessions[1]$hmm[r,best_p,drop=FALSE])^2)))
}
mat2_diff_vec = c()
for (mat2 in mat2_list_res_1000_sessions[2:length(mat2_list_res_1000_sessions)]) {
  mat2_diff_vec = c(mat2_diff_vec,sqrt(sum((mat2[best_p,,drop=FALSE]-mat2_list_res_1000_sessions[1]$hmm[best_p,,drop=FALSE])^2)))
}
fwrite(list(mat1_diff_vec),file="exp2_sorted_new_trans_mat_boxplot.txt")
fwrite(list(mat2_diff_vec),file="exp2_sorted_new_emiss_mat_boxplot.txt")
fwrite(list(vec_diff_vec),file="exp2_sorted_new_start_prob_boxplot.txt")


boxplot(mat1_diff_vec)
boxplot(mat2_diff_vec)
boxplot(vec_diff_vec)