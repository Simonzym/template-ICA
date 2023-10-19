library(combinat)
library(gtools)
library(Matrix)
#dimension
#voxels
V = 1000
#ICs
Q = 10
#subjects
N = 100
T = 20
#observation
y = list()
for(v in 1:V){
  y[[v]] = matrix(1, nrow = N*T, ncol = 1)
}


#hyper-parameter
#variance of first-level model (ICA)
sigma2 = 0.001
#number of mixture normals
m = 2

#initilize parameters
#mixtures weights for S
pis = matrix(c(0.4, 0.6), nrow = Q, ncol = m, byrow = T)
mus = matrix(1:2, nrow = Q, ncol = m, byrow = T)
sigma2s = matrix(c(0.1, 0.2), nrow = Q, ncol = m, byrow = T)
#mixing matrix
M = matrix(1, nrow = T*N, ncol = Q)
#covairance of ICA model
Sig = kronecker(diag(1,N), sigma2*diag(1, T))
Sig_inv = kronecker(diag(1,N), sigma2^(-1)*diag(1, T))

#mixture weights for conditional of s,z|y
pis_c = list()
mus_c = list()
sigma2s_c = list()

for(v in 1:V){
  pis_c[[v]] = matrix(c(0.4, 0.6), nrow = Q, ncol = m, byrow = T)
  mus_c[[v]] = matrix(1:2, nrow = Q, ncol = m, byrow = T)
  sigma2s_c[[v]] = matrix(c(0.1, 0.2), nrow = Q, ncol = m, byrow = T)
}


#update parameters
expect_s = list()
expect_ss = list()
res = permutations(n=m, r=Q, v = seq(1:m), repeats.allowed = T)
for(v in 1:V){
  for(i in 1:nrow(res)){
    pos = matrix(c(1:Q, res[i,]), nrow = Q)
    pi_i = pis_c[[v]][pos]
    mu_i = mus_c[[v]][pos]
    sigma2s_i = sigma2s_c[[v]][pos]
    expect_s[[v]] = expect_s[[v]] + prod(pi_i)*mu_i
    expect_ss[[v]] = expect_ss[[v]] + prod(pi_i)*(diag(sigma2s_i)+mu_i%*%t(mu_i))
}
}
M1 = 0
M2 = 0
pi_sum = 0
pi_mu = 0
pi_sigma = 0
for(v in 1:V){
  M1 = M1 + y[[v]]%*%expect_s[[v]]
  M2 = M2 + expect_ss[[v]]
  pi_sum = pi_sum + pis_c[[v]]
  pi_mu = pi_mu + pis_c[[v]]*mus_c[[v]]
  pi_sigma = pi_sigma + pis_c[[v]] * (sigma2s_c[[v]] + mus_c[[v]]^2 - 2*mus*mus_c[[v]] + mus^2)
}
M = M1%*%solve(M2)
pis = pi_sum/V
mus = pi_mu/pi_sum
sigma2s = pi_sigma/pi_sum

#update parameters of s,z|y
gamma_c = list()
for(v in 1:V){
  gamma_c[[v]] = pis*sqrt(sigma2s_c[[v]]/sigma2s)*
    exp(0.5*(mus_c[[v]]^2/sigma2s_c[[v]] - (mus^2+sigma2s_c[[v]])/sigma2s))
}

for(v in 1:V){
  pis_c[[v]] = gamma_c[[v]]/apply(gamma_c[[v]],1,sum)
}

B = t(M)%*%Sig_inv%*%M
for(v in 1:V){
  sigma2s_c[[v]] = (diag(B) + 1/sigma2s)^(-1)
}

b = list()
for(v in 1:V){
  b[[v]] = t(M)%*%Sig_inv%*%y[[v]]
}


mat_b = list()
mu_sigma = c(t(mus))/c(t(sigma2s))
for(v in 1:V){
  mul = matrix(0, nrow = m*Q, ncol = m*Q)
  for(q in 1:Q){

      bs = B[1,]*sigma2s_c[[v]]
      bs = c(t(bs))
      blocks = matrix(rep(bs, m), nrow = m, byrow = T)
      mul[((q-1)*m+1):(q*m), ] = blocks
      mul[((q-1)*m+1):(q*m), ((q-1)*m+1):(q*m)] = diag(sigma2s[q,])
      
      
  }
  mat_b[[v]] = mul
  bb = rep(b[[v]], each = m)
  mu_b = solve(mul)%*%(bb+mu_sigma)
  mus_c[[v]] = matrix(mu_b, nrow = Q, byrow = T)
}


#stop criteria




