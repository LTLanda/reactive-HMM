library(tidyverse)
library(magrittr) 
library(dplyr)
library(foreach)
library(doParallel)
library(profvis)
library(useful)
source('run_sim.R')
library(data.table)
library(readr)
library(parallel)
library(combinat)
library(e1071)
library(ROCR)
library(tree)
library(ISLR)
library(rpart)
library(rpart.plot)
library(randomForest)
library(gbm)
library(caret)
library(MASS)
library(class)
library(caret)
library(ROCR)
library(pROC)
library(MASS)
library(e1071)
library(party)
library(pROC)
library(RColorBrewer)
library(ROCR)
library(class)
library(rpart)
library(rattle)
library(rpart.plot)
library(ISLR)
library(lattice)
library(ggplot2)
library(plyr)
library(ggplot2)


initHMM = function(States, Symbols, nAgentResponse, startProbs=NULL, transProbs=NULL,
                   emissionProbs=NULL)
{
  nStates    = length(States)
  nSymbols   = length(Symbols)
  nAR = nAgentResponse
  S          = rep(1/nStates,nStates)
  T          = 1*diag(0, nAgentResponse*nStates, nStates) + array(1/(nStates),c(nAgentResponse*nStates,nStates))
  E          = array(1/(nSymbols),c(nStates,nSymbols))
  names(S)   = States
  dimnames(T)= list(from=1:(nAgentResponse*nStates),to=States)
  dimnames(E)= list(states=States,symbols=Symbols)
  if(!is.null(startProbs)){S[]  = startProbs[]}
  if(!is.null(transProbs)){T[,] = transProbs[,]}
  if(!is.null(emissionProbs)){E[,] = emissionProbs[,]}
  if(!is.null(nAgentResponse)){nAR = nAgentResponse}
  return(list(States=States,Symbols=Symbols,startProbs=S,transProbs=T,
              emissionProbs=E, nAgentResponse = nAR))
}

baumWelchRecursion = function(hmm, observation, agent_obs, bigger_i ,iteration)
{
  TransitionMatrix    = hmm$transProbs
  TransitionMatrix[,] = 0
  EmissionMatrix      = hmm$emissionProbs
  EmissionMatrix[,]   = 0
  StartProbs = hmm$startProbs
  StartProbs[] = 0
  likelihood = 0
  f = forward(hmm,  observation, agent_obs)
  b = backward(hmm, observation, agent_obs)
  nAgentResponse = hmm$nAgentResponse
  probObservations = -Inf
  for(i in 1:length(hmm$States))
  {
    for (k in 1:nAgentResponse) {
      j = f[i,length(observation),k]
      probObservations = log(exp(j) + exp(probObservations))
    }
  }
  sum_of_gammas = array(NA,c(length(hmm$States),nAgentResponse))
  for (state in hmm$States) {
    for (k in 1:nAgentResponse) {
      s = sum(exp(head(f[state,,k],-1)+head(b[state,,k],-1)))
      if(s == 0){
        sum_of_gammas[state,k] = -Inf
      }
      else{
        sum_of_gammas[state,k] = log(s)
      }
    }
  }
  for (kt in 1:nAgentResponse) {
    for(x in hmm$States)
    {
      for(y in hmm$States)
      {
        xi_sum = 0
        temp = -Inf
        for(i in 1:(length(observation)-1))
        {
          sum_k_tag = 0
          for (k in 1:nAgentResponse) {
            j = f[x,i,kt] + log(hmm$transProbs[(((x-1)*nAgentResponse+1) : (x*nAgentResponse))[kt],y]) +
              log(hmm$emissionProbs[y,(((observation[i+1]-1)*nAgentResponse+1 ): (observation[i+1]*nAgentResponse))[k]]) + b[y,i+1,k]
            sum_k_tag = sum_k_tag + exp(j)
          }
          xi_sum = xi_sum + sum_k_tag
        }
        if (sum_of_gammas[x,kt] == -Inf) {
          temp = 0
        }
        else{
          temp = exp(log(xi_sum) - sum_of_gammas[x,kt])
        }
        TransitionMatrix[(((x-1)*nAgentResponse+1):(x*nAgentResponse))[kt],y] = temp
      }
    }
  }
  
  for(x in hmm$States)
  {
    for(s in hmm$Symbols)
    {
      temp = -Inf
      for(i in 1:length(observation))
      {
        if(s == get_state(observation[i],agent_obs[i],nAgentResponse))
        {
          j = f[x,i,agent_obs[i]] + b[x,i,agent_obs[i]]
          temp = log(exp(j) + exp(temp))
        }
      }
      if (sum(exp(sum_of_gammas[x,])) == 0){
        temp = 0
      }
      else{
        temp = exp(temp - log(sum(exp(sum_of_gammas[x,]))))
      }
      EmissionMatrix[x,s] = temp
    }
  }
  
  for (x in hmm$States) {
    temp = -Inf
    for (k in 1:nAgentResponse) {
      j = f[x,1,k] + b[x,1,k]
      if(j > - Inf)
      {
        temp = j + log(1+exp(temp-j))
      }
    }
    if (probObservations==-Inf) {
      temp = 0
    }
    else{
      temp = exp(temp - probObservations)
    }
    StartProbs[x] = temp
  }
  
  for (x in hmm$States) {
    temp = 0
    for (k in 1:nAgentResponse) {
      temp = temp + exp(b[x,1,k])
    }
    temp = temp * StartProbs[x] * EmissionMatrix[x,get_state(observation[1],agent_obs[1],nAgentResponse)]
    likelihood = likelihood+temp
  }
  return(list(TransitionMatrix=TransitionMatrix,EmissionMatrix=EmissionMatrix,StartProbs=StartProbs, likelihood=likelihood))
}

Noisify <- function(data) {
  
  if (is.vector(data)) {
    noise <- rnorm(length(data), 0, 0.0005)
    noisified <- data + noise
  } else {
    length <- dim(data)[1] * dim(data)[2]
    noise <- matrix(rnorm(length, 0, 0.0005), dim(data)[1])
    noisified <- data + noise
  }
  return(noisified)
}

main = function(nTests, originalHmm){
  l <- (1:nTests) %>% map(runTestMain)
  
  mat1_list <- l %>% map(~.$hmm$transProbs)
  mat2_list <- l %>% map(~.$hmm$emissionProbs)
  vec_list <- l %>% map(~.$hmm$startProbs)
  
  find_mean <- function(some_list) {
    reduce(some_list, `+`) / length(some_list)
  }
  
  mat1_diff_vec = c()
  for (mat1 in mat1_list) {
    mat1_diff_vec = c(mat1_diff_vec,sum((mat1-originalHmm$transProbs)^2))
  }
  mat2_diff_vec = c()
  for (mat2 in mat2_list) {
    mat2_diff_vec = c(mat2_diff_vec,sum((mat2-originalHmm$emissionProbs)^2))
  }
  vec_diff_vec = c()
  for (vec in vec_list) {
    vec_diff_vec = c(vec_diff_vec,sum((vec-originalHmm$startProbs)^2))
  }
  
  fwrite(list(mat1_diff_vec),file="stability_of_estimates_trans_mat_boxplot.txt")
  fwrite(list(mat2_diff_vec),file="stability_of_estimates_emiss_mat_boxplot.txt")
  fwrite(list(vec_diff_vec),file="stability_of_estimates_start_prob_boxplot.txt")
  
  boxplot(mat1_diff_vec)
  boxplot(mat2_diff_vec)
  boxplot(vec_diff_vec)
  
  mat1_mean <- find_mean(mat1_list)
  mat2_mean <- find_mean(mat2_list)
  vec_mean <- find_mean(vec_list)
  
  obj_mean = list(mat1 = mat1_mean, mat2 = mat2_mean, vec = vec_mean)
  
  return(obj_mean)
}
runTestMain = function(seed){
  setting_param_list = list(seed=seed, n_sessions=1000, start_point=-1)
  other_param_list = list(original_hmm=hmm)
  seed <-  setting_param_list$seed
  n_sessions <-  setting_param_list$n_sessions
  startPoint <-  setting_param_list$start_point
  
  original_hmm <- other_param_list$original_hmm
  n_hidden_states <- length(original_hmm$States)
  n_agent_response <- length(original_hmm$Symbols)/length(original_hmm$States)
  print(seed)
  obs = simHMMwithSeed(original_hmm, n_sessions,seed)
  likelihood = -Inf
  bestpoint = 0
  print('finished simming starting fitting')
  startPointhmm = original_hmm
  if (startPoint == 1) {
    startPointhmm = initHMM(1:3, 1:9, nAgentResponse = 3,
                            startProbs = c(1/3,1/3,1/3),
                            transProbs = matrix(c(1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3),9),
                            emissionProbs = matrix(c(0.3,1/60,1/300,0.3,1/60,1/300,0.3,1/60,1/300,0.03,0.3,1/60,0.03,0.3,1/60,0.03,0.3,1/60,1/300,1/60,47/150,1/300,1/60,47/150,1/300,1/60,47/150),3))
    
  }
  if (startPoint == 2) {
    startPointhmm = initHMM(1:3, 1:9, 
                            startProbs = c(1/3,1/3,1/3),
                            transProbs = matrix(c(1/2,1/2,1/2,1/4,1/4,1/4,1/4,1/4,1/4,1/4,1/4,1/4,1/2,1/2,1/2,1/4,1/4,1/4,1/4,1/4,1/4,1/4,1/4,1/4,1/2,1/2,1/2),9))
    
  }
  if (startPoint == 3) {
    startPointhmm = initHMM(1:3, 1:9, 
                            startProbs = c(1/3,1/3,1/3),
                            transProbs = matrix(c(1/2,1/4,1/4,1/2,1/4,1/4,1/2,1/4,1/4,1/4,1/2,1/4,1/4,1/2,1/4,1/4,1/2,1/4,1/4,1/4,1/2,1/4,1/4,1/2,1/4,1/4,1/2),9))
    
  }
  
  temp = baumWelch(hmm = startPointhmm, obs$observation, obs$agent_obs,maxIterations = 100,lengths = obs$lengths)
  print(temp$hmm$likelihood)
  return(temp)
}

runTest = function(setting_param_list,other_param_list){
  seed <-  setting_param_list$seed
  n_sessions <-  setting_param_list$n_sessions
  startPoint <-  setting_param_list$start_point
  
  original_hmm <- other_param_list$original_hmm
  n_hidden_states <- length(original_hmm$States)
  n_agent_response <- original_hmm$nAgentResponse
  
  obs = simHMMwithSeed(original_hmm, n_sessions,seed)
  likelihood = -Inf
  bestpoint = 0
  
  if (startPoint == 1) {
    startPointhmm = initHMM(1:3, 1:9, nAgentResponse = 3,
                            startProbs = c(1/3,1/3,1/3),
                            transProbs = matrix(c(1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3),9),
                            emissionProbs = matrix(c(13/54,1/18,1/27,13/54,1/18,1/27,13/54,1/18,1/27,1/18,13/54,1/18,1/18,13/54,1/18,1/18,13/54,1/18,1/27,1/27,13/54,1/27,1/27,13/54,1/27,1/27,13/54),3))
    
  }
  if (startPoint > 1) {
    set.seed(startPoint)
    startPointhmm = initHMM(1:3, 1:9, nAgentResponse = 3,
                            startProbs = c(1/3,1/3,1/3),
                            transProbs = matrix(rnorm(27),9),
                            emissionProbs = matrix(rnorm(27),3))
    #normalization
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    
    startPointhmm$transProbs = range01(startPointhmm$transProbs)
    startPointhmm$transProbs[startPointhmm$transProbs==0] = 0.0000001
    startPointhmm$transProbs = (startPointhmm$transProbs/apply(startPointhmm$transProbs,1,sum))
    startPointhmm$emissionProbs = range01(startPointhmm$emissionProbs)
    startPointhmm$emissionProbs[startPointhmm$emissionProbs==0] = 0.0000001
    startPointhmm$emissionProbs = (startPointhmm$emissionProbs/apply(startPointhmm$emissionProbs,1,sum))
    
  }
  temp = baumWelch(startPointhmm, obs$observation, obs$agent_obs,maxIterations = 200,lengths = obs$lengths)
  print(temp$hmm$likelihood)
  bm = temp
  likelihood = temp$hmm$likelihood
  bestpoint = startPoint
  return(bm)
}

simHMM = function(hmm, nSessions)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  states   = c()
  emission = c()
  agent_obs = c()
  lengths = c()
  nAgentResponse = hmm$nAgentResponse
  for (n in 1:nSessions) {
    length = round(rnorm(1,12,2))
    lengths = c(lengths, length)
    states   = c(states, sample(hmm$States,1,prob=hmm$startProbs))
    for(i in 1:length)
    {
      emi      = get_emi(sample(hmm$Symbols, 1, prob=hmm$emissionProbs[states[i],]))
      emission = c(emission, emi[1])
      agent_obs = c(agent_obs, emi[2])
      if (i<length) {
        state  = sample(hmm$States, 1, prob=hmm$transProbs[((states[i]-1)*nAgentResponse +1) : (states[i]*nAgentResponse),])
        states = c(states, state)
      }
      
    }
  }
  
  return(list(states=states,observation=emission, agent_obs = agent_obs, lengths = lengths))
}

simHMMwithSeed = function(hmm, nSessions,seed)
{
  set.seed(seed)
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  states   = c()
  emission = c()
  agent_obs = c()
  lengths = c()
  nAgentResponse = hmm$nAgentResponse
  for (n in 1:nSessions) {
    length = round(rnorm(1,6,4))
    if (length<=0) {
      length = 2
    }
    lengths = c(lengths, length)
    states   = c(states, sample(hmm$States,1,prob=hmm$startProbs))
    for(i in 1:length)
    {
      if (length(hmm$emissionProbs[states[i],])!=length(hmm$Symbols)) {
        print(i)
      }
      emi      = get_emi(sample(hmm$Symbols, 1, prob=hmm$emissionProbs[states[i],]),hmm)
      emission = c(emission, emi[1])
      agent_obs = c(agent_obs, emi[2])
      if (i<length) {
        state  = sample(hmm$States, 1, prob=hmm$transProbs[get_state(states[i], agent_obs[i],nAgentResponse),])
        states = c(states, state)
      }
      
    }
  }
  
  return(list(states=states,observation=emission, agent_obs = agent_obs, lengths = lengths))
}

simHMMwithSeedandAction = function(hmm, nSessions,seed,action=-1,sent=-1)
{
  set.seed(seed)
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  states   = c()
  emission = c()
  agent_obs = c()
  lengths = c()
  nAgentResponse = hmm$nAgentResponse
  for (n in 1:nSessions) {
    length = 10
    lengths = c(lengths, length)
    states   = c(states, sample(hmm$States,1,prob=hmm$startProbs))
    for(i in 1:length)
    {
      emi      = get_emi(sample(hmm$Symbols, 1, prob=hmm$emissionProbs[states[i],]),hmm)
      if (sent!=-1) {
        emission = c(emission,sent)
      }
      else{
        emission = c(emission, emi[1])
      }
      if (action!=-1) {
        agent_obs = c(agent_obs, action)
      }
      else{
        agent_obs = c(agent_obs, emi[2])
      }
      if (i<length) {
        state  = sample(hmm$States, 1, prob=hmm$transProbs[get_state(states[i], agent_obs[i],nAgentResponse),])
        states = c(states, state)
      }
      
    }
  }
  
  return(list(states=states,observation=emission, agent_obs = agent_obs, lengths = lengths))
}

get_state = function(state, agent, nAgentResponse) {
  return ((((state-1)*nAgentResponse+1):(state*nAgentResponse))[agent])
}

get_emi = function(emi,hmm){
  nAgentResponse = hmm$nAgentResponse
  nSents = length(hmm$Symbols)/nAgentResponse
  for (sentiment in 1:nSents) {
    if (emi <= (sentiment*nAgentResponse) && emi >=((sentiment-1)*nAgentResponse+1)){
      if (emi == (sentiment*nAgentResponse)) {
        return (c(sentiment,nAgentResponse))
      }
      else{
        return (c(sentiment, emi%%nAgentResponse)) 
      }
    }
  }
}

viterbiMaster = function(hmm, obs, agent_res, lengths=NULL){
  if (length(lengths)!=0) {
    indices = get_indeces(lengths)
    paths = c()
    for (col in 1:(dim(indices)[2])) {
      if(length(obs[indices[,col][1]:indices[,col][2]])>1){
        paths = c(paths, viterbi(hmm,obs[indices[,col][1]:indices[,col][2]],agent_res[indices[,col][1]:indices[,col][2]]))
      }
    }
    
  }
  return(paths)
}

viterbi = function(hmm, observation, agent_response)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations  = length(observation)
  nStates    = length(hmm$States)
  nAgentResponse = hmm$nAgentResponse
  v          = array(NA,c(nStates,nObservations))
  prev          = array(NA,c(nStates,nObservations))
  dimnames(v)= list(states=hmm$States,index=1:nObservations)
  # Init
  for(state in hmm$States)
  {
    for (k in 1:nAgentResponse) {
      if (agent_response[1] == k) {
        v[state,1] = log(hmm$startProbs[state] * hmm$emissionProbs[state,(((observation[1]-1)*nAgentResponse+1) : (observation[1]*nAgentResponse))[k]])
        prev[state,1] = NaN      
      }
    }
  }
  # Iteration
  for(k in 2:nObservations)
  {
    for(state in hmm$States)
    {
      for (kj in 1:nAgentResponse) {
        maxi = -Inf
        best_prev = NaN
        if (agent_response[k]==kj) {
          for(previousState in hmm$States)
          {
            for (ki in 1:nAgentResponse) {
              if (agent_response[k-1]==ki) {
                temp = v[previousState,k-1] + log(hmm$transProbs[(((previousState-1)*nAgentResponse+1) : (previousState*nAgentResponse))[ki],state]) 
                if (temp>maxi){
                  maxi = temp
                  best_prev = previousState
                }
              }
            }
          }
          v[state,k] = log(hmm$emissionProbs[state,(((observation[k]-1)*nAgentResponse+1) : (observation[k]*nAgentResponse))[kj]]) + maxi
          prev[state,k]=best_prev
        }
      }
    }
  }
  # Traceback
  print(v)
  print(prev)
  backback = FALSE
  viterbiPath = rep(NA,nObservations)
  for(state in hmm$States)
  {
    for (kFin in 1:nAgentResponse) {
      if (agent_response[nObservations]==kFin) {
        if(max(v[,nObservations])==v[state,nObservations])
        {
          viterbiPath[nObservations] = state
          previous = state
          backback = TRUE
          break
        }  
      }
    }
    if (backback) {
      break
    }
  }
  
  for(k in (nObservations-1):1)
  {
    backback = FALSE
    for(state in hmm$States)
    {
      if(prev[previous,k+1]==state){ 
        viterbiPath[k] = state
        previous = state
        backback = TRUE
        break
      } 
    }
  }
  return(viterbiPath)
}

get_matrix_rows_from_agent_response = function(hmm,k){
  rows = c()
  for (state in hmm$States) {
    rows = c(rows,(((state-1)*hmm$nAgentResponse+1) : (state*hmm$nAgentResponse))[k])
  }
  return(rows)
}

forward = function(hmm, observation, agent_obs)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations  = length(observation)
  nStates    = length(hmm$States)
  if (is.na(observation[length(observation)])) {
    print(observation)
  }
  nAgentResponse = hmm$nAgentResponse
  f          = array(NA,c(nStates,nObservations,nAgentResponse))
  dimnames(f)= list(states=hmm$States,index=1:nObservations, agentResponse = 1:nAgentResponse)
  # Init
  for(state in hmm$States)
  {
    for(k in 1:nAgentResponse)
    {
      if(agent_obs[1]==k)
      {
        f[state,1,k] = log(hmm$startProbs[state] * hmm$emissionProbs[state,(((observation[1]-1)*nAgentResponse+1) : (observation[1]*nAgentResponse))[k]])
      }
      else
      {
        f[state,1,k] = -Inf
      }
    }
  }
  # Iteration
  for(k in 2:nObservations)
  {
    for(state in hmm$States)
    {
      for(kj in 1:nAgentResponse)
      {
        if(agent_obs[k]==kj)
        {
          logsum = -Inf
          for(previousState in hmm$States)
          {
            for(ki in 1:nAgentResponse)
            {
              if(agent_obs[k-1]==ki)
              {
                temp   = f[previousState,k-1,ki] + log(hmm$transProbs[(((previousState-1)*nAgentResponse+1) : (previousState*nAgentResponse))[ki],state])
                logsum = log(exp(temp) + exp(logsum))
              }
            }
          }
          f[state,k,kj] = log(hmm$emissionProbs[state,(((observation[k]-1)*nAgentResponse+1) : (observation[k]*nAgentResponse))[kj]]) + logsum
        }
        else
        {
          f[state,k,kj] = -Inf
        }
      }
    }
  }
  return(f)
}

convert_i_1 = function(state)
{
  if(state == 1)
  {
    return(2)
  }
  if(state == 2)
  {
    return(1)
  }
  if(state == 3)
  {
    return(0)
  }
}

backward = function(hmm, observation, agent_obs)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations  = length(observation)
  nStates    = length(hmm$States)
  nAgentResponse = hmm$nAgentResponse
  b          = array(NA,c(nStates,nObservations,nAgentResponse))
  dimnames(b)= list(states=hmm$States,index=1:nObservations, agentResponse = 1:nAgentResponse)
  # Init
  for(state in hmm$States)
  {
    for(k in 1:nAgentResponse)
    {
      if(agent_obs[nObservations]==k)
      {
        b[state,nObservations,k] = log(1)
      }
      else
      {
        b[state,nObservations,k] = -Inf
      }
    }
  }
  # Iteration
  for(k in (nObservations-1):1)
  {
    for(state in hmm$States)
    {
      for(ki in 1:nAgentResponse)
      {
        if(agent_obs[k]==ki)
        {
          logsum = -Inf
          for(nextState in hmm$States)
          {
            for (kj in 1:nAgentResponse) {
              if (agent_obs[k+1]==kj) {
                temp   = b[nextState,k+1,kj] + log(hmm$transProbs[(((state-1)*nAgentResponse+1) : (state*nAgentResponse))[ki],nextState]*hmm$emissionProbs[nextState,(((observation[k+1]-1)*nAgentResponse+1) : (observation[k+1]*nAgentResponse))[kj]])
                logsum = log(exp(temp) + exp(logsum))
              }
            }
          }
          b[state,k,ki] = logsum    
        }
        else
        {
          b[state,k,ki] = -Inf
        }
      }
      
    }
  }
  return(b)
}


posterior = function(hmm, observation)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  f = forward(hmm, observation)
  b = backward(hmm, observation)
  probObservations = f[1,length(observation)]
  for(i in 2:length(hmm$States))
  {
    j = f[i,length(observation)]
    if(j > - Inf)
    {
      probObservations = j + log(1+exp(probObservations-j))
    }
  }
  posteriorProb = exp((f+b)-probObservations)
  return(posteriorProb)
}

baumWelch = function(hmm, observation, agent_obs, maxIterations=100, delta=1E-9, pseudoCount=0, lengths=c(),fromObs=1,toObs=-1,realData=FALSE)
{
  tempHmm = hmm
  tempHmm$transProbs[is.na(hmm$transProbs)]       = 0
  tempHmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  tempHmm$likelihood = 0
  diff = c()
  indices = get_indeces(lengths)
  if (toObs==-1) {
    toObs = (dim(indices)[2])
  }
  for(i in 1:maxIterations)
  {
    T    = tempHmm$transProbs
    T[,] = 0
    E      = tempHmm$emissionProbs
    E[,]   = 0
    S = tempHmm$startProbs
    S[] = 0
    L = 0
    if (length(lengths)!=0) {
      
      iteration = 1
      for (col in fromObs:toObs) {
        if(length(observation[indices[,col][1]:indices[,col][2]])>1){
          bw = baumWelchRecursion(tempHmm, observation[indices[,col][1]:indices[,col][2]], agent_obs[indices[,col][1]:indices[,col][2]],i, iteration)
          T  = T + bw$TransitionMatrix
          E  = E + bw$EmissionMatrix
          S  = S + bw$StartProb
          L = L + log(bw$likelihood)
          iteration = iteration+1
        }
      }
      T = T/length(lengths)
      E = E/length(lengths)
      S = S/length(lengths)
      for (row in 1:(dim(T)[1])){
        zero_vector = compare.list(T[row,], double(dim(T)[2]))
        if (zero_vector && zero_vector){
          T[row,] = integer(dim(T)[2]) + 1 
        }
      }
      T = (T/apply(T,1,sum))
      E[E==0] = 0.0000001
      E[E==Inf] = 9e+303
      E = (E/apply(E,1,sum))
      S = (S/sum(S))
      d = sqrt(sum((tempHmm$transProbs-T)^2)) + sqrt(sum((tempHmm$emissionProbs-E)^2))+ sqrt(sum((tempHmm$startProbs-S)^2))
      diff = c(diff, d)
      tempHmm$transProbs    = T
      tempHmm$emissionProbs = E
      tempHmm$startProbs = S
      tempHmm$likelihood = L
      if(d < delta)
      {
        break
      }
    }
    if(d < delta)
    {
      break
    }
  }
  tempHmm$transProbs[is.na(hmm$transProbs)]       = NA
  tempHmm$emissionProbs[is.na(hmm$emissionProbs)] = NA
  return(list(hmm=tempHmm,difference=diff))
}

get_indeces = function(lengths){
  indices = matrix(c(1,lengths[1]),1)
  for (l in lengths[-1]) {
    indices = matrix(c(indices, matrix(c(indices[length(indices)]+1,indices[length(indices)]+l),1)),2)
  }
  return(indices)
}