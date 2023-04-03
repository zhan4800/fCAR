function delta=MH_stepsize1(theta_sd,theta_mle,multiple)
% This function compute the step size for Metropolis-Hastings update(proposal distribution) 
% based on MLE estimations of Theta and their standard deviation estimates.
% Input:
%       theta_sd:  standard error.
%       theta_mle: MLE estimate.
%       multiple: the multiple of the standard deviation to control step
%       size.
% Output:
%       delta: the step size used for log normal proposal. 
%       log sigma^2_new=log sigma^2_old+eps, where eps~N(0, delta^2), for
%       all j,k. 
%delta=sqrt(log(0.5*(1+sqrt(1+4*theta_sd.^2*multiple^2./theta_mle.^2))));
delta=multiple*sqrt(log(0.5*(1+sqrt(1+4*theta_sd.^2./theta_mle.^2))));
