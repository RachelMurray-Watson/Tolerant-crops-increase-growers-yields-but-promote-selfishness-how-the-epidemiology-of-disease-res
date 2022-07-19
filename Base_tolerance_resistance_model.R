par(las = 1, mfrow = c(1,2), mar = c(5.1, 5.1, 3.1, 4.1), family = "Times New Roman")
library(deSolve)
library(fields)
library(colorspace)
pal <- sequential_hcl(20,"BurgYl")
pal <- rev(pal)
Resistance_tolerance_model <-  function(time, y, params){
  S_C <- y[1]
  E_C <- y[2]
  I_C <- y[3]
  S_U <- y[4]
  E_U <- y[5]
  I_U <- y[6]
  
  N <- params$N
  
  beta <- params$beta
  beta_C <- params$delta_beta * params$beta
  gamma <- params$gamma 
  delta_sigma <- params$delta_sigma
  mu <- params$mu
  mu_C <- params$mu_C
  epsilon <- params$epsilon
  epsilon_C <- params$delta_epsilon *params$epsilon
  YU <- params$Y
  YC <- params$delta_Y * params$Y
  phi_C <- params$phiC
  phi_R <- params$phiR
  eta <- params$eta
  LU <- params$L
  LC <- params$delta_L * params$L

  P_SU <- P_EU <- YU
  P_SC <- P_EC <-  YC - phi_C
  
  P_IUR <- YU - phi_R*LU 
  P_ICR <- YC - phi_C - phi_R*LC
  P_IUH <- YU - LU 
  P_ICH <- YC - phi_C - LC 
  
   q_U <-beta*(delta_sigma*I_C + I_U)/(beta*(delta_sigma*I_C + I_U) + gamma)
  q_C <- beta_C*(delta_sigma*I_C + I_U)/(beta_C*(delta_sigma*I_C + I_U) + gamma)
  
  P_U <- YU - q_U*(epsilon/(epsilon + gamma))*LU*(gamma/( mu + gamma) + mu/(gamma + mu)*phi_R)
  P_C <- YC - phi_C - q_C*(epsilon_C/(epsilon_C + gamma))*LC*(gamma/(mu_C + gamma) + mu_C/(mu_C + gamma)*phi_R)
  z_SC <- max(0, 1 - exp(-eta*(P_U - P_SC)))
  z_EC <- max(0, 1 - exp(-eta*(P_U - P_EC)))
  z_ICR <- max(0, 1 - exp(-eta*(P_U - P_ICR)))
  z_ICH <- max(0, 1 - exp(-eta*(P_U - P_ICH)))
  
  z_SU <- 0
  z_EU <- max(0, 1 - exp(-eta*(P_C - P_EU)))
  z_IUR <- max(0, 1 - exp(-eta*(P_C - P_IUR)))
  z_IUH <- max(0, 1 - exp(-eta*(P_C - P_IUH)))
  
  
  dS_C <- gamma*(S_C*(1- z_SC) + I_C*(1 - z_ICH) + E_C*(1-z_EC) + E_U*z_EU + I_U*z_IUH)  + mu_C*I_C*(1-z_ICR) + I_U*z_IUR*mu - beta_C*(I_U + delta_sigma*I_C)*S_C - gamma*S_C 
  dE_C <-  beta_C*(I_U + delta_sigma*I_C)*S_C - gamma*E_C - epsilon_C*E_C
  dI_C <- epsilon_C*E_C - gamma*I_C - mu_C*I_C
  dS_U <- gamma*(S_U + E_U*(1-z_EU) + I_U*(1-z_IUH) + S_C*(z_SC)  + I_C*(z_ICH) + E_C*(z_EC))+ I_U*(1-z_IUR)*mu + mu_C*I_C*z_ICR - beta*(I_U + delta_sigma*I_C)*S_U - gamma*S_U 
  dE_U <-  beta*(I_U + delta_sigma*I_C)*S_U - gamma*E_U  - epsilon*E_U
  dI_U <- epsilon*E_U - gamma*I_U - mu*I_U 
  
  
  
  return(list(c(dS_C, dE_C, dI_C, dS_U, dE_U,dI_U)))
  
}
readParams <- function(N,beta,delta_beta,gamma,Y,delta_Y,phiC,phiR,L,delta_L,nu,delta_nu,Delta = Delta,eta = eta,delta_sigma = delta_sigma,epsilon = epsilon,delta_epsilon = delta_epsilon)
{
  retval <- list(N = N, 
                 beta = beta/N, 
                 delta_beta = delta_beta,
                 gamma = gamma, 
                 Y = Y,
                 delta_Y = delta_Y, 
                 phiC= phiC,
                 phiR = phiR, 
                 L = L, 
                 delta_L = delta_L, 
                 nu = nu, 
                 delta_nu = delta_nu,
                 Delta = Delta,
                 mu = 1/(((1/nu) - 1/2)*Delta), 
                 mu_C = 1/(((1/(delta_nu*nu)) - 1/2)*Delta),
                 eta = eta,
                 delta_sigma = delta_sigma,
                 epsilon = epsilon,
                 delta_epsilon = delta_epsilon)
  return(retval)
}

#Default tolerance parameters
parameters <- readParams(beta = .055, delta_beta = 1,
                         gamma = 1/120,N = 1, Y = 1, delta_Y = 1, delta_L = 0.1,
                         phiR = 0.7, phiC = 0.1, 
                         L = 0.6, eta = 10, nu = 1, delta_nu = 0.1, delta_sigma = 1, epsilon = 1/40, delta_epsilon = 1, Delta = 120) 



#Default resistance parameters
parameters <- readParams(beta = .055, delta_beta = 0.5,
                         gamma = 1/120,N = 1, Y = 1, delta_Y = 1, delta_L = 1,
                         phiR = 0.7, phiC = 0.1, 
                         L = 0.6, eta = 10, nu = 1, delta_nu =1, delta_sigma = 0.5, epsilon = 1/40, delta_epsilon = 2, Delta = 120) 


# Initial conditions - S_C, E_C, I_C, S_U, E_U, I_U
yini <- c(0.1,0,0,(1-0.1 - 0.01),0,0.01)
# 100 seasons
times <- seq(0,(100*120),10)

out <- ode(func =Resistance_tolerance_model ,  y = yini, time = times, parms = parameters)

plot(out[,1], out[,2], lwd = 4, ty = "l", ylim = c(0,1), col = "coral", main = "", ylab = expression("Proportion of growers"), xlab = expression("Time (seasons)"), xaxt = "n"   )
axis(side = 1, at = seq(0, max(times ), max(times /5)),labels = seq(0, max(times /300), max(times /300)/5), cex.axis = 1.6)
lines(out[,1], out[,3], lwd = 4, col = "darkgoldenrod1")
lines(out[,1], out[,4], lwd = 4, col = "firebrick")
lines(out[,1], out[,5], lwd = 4, col = "dodgerblue3")
lines(out[,1], out[,6], lwd = 4, lty = 1, col = "blue")
lines(out[,1], out[,7], lwd = 4, lty = 1, col = "light blue")
legend("topright", c(expression("S"[U]*""),expression("E"[U]*""),expression("I"[U]*""), expression("S"[C]*""),expression("E"[C]*""),expression("I"[C]*"")), ncol = 2, lty = 1, lwd = 3, col = c("dodgerblue3", "blue", "lightblue","coral", "orange", "firebrick"), cex = 1.35)

