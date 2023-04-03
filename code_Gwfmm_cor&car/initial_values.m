function [tau,nu2,PI,alpha,a_nu2,b_nu2,a_pi,b_pi]=initial_values(betans,Vbetans,model,wavespecs,MCMCspecs,meanop)
%% This function compute initial values for tau, nu2 (scale mixture models para) and pi.
% Input:
%       Theta_mle0: MLE estimates of Theta. 
%       Vbetans0: MLE estimates of Beta.
%       model:
% Output:
%       tau: p x K variance/scaling parameter of beta.
%       nu2: p x J the lasso parameters, depend on i and j, but not k.
%       pi: p x J the probability parameters for beta equal to zero, depend
%       on i and j, but not k.
%       a_nu2, b_nu2: hyperparameters for the gamma prior on nu2
%       a_pi, b_pi: hyperparameters for the beta prior on pi
%
%
%  Functions needed: repvec.m : same as Splus function rep applied to vectors
%
%%#zhu#% This code is written by assuming that \pi_{ij} depend on
%%only i and j, not k. \tau_{ijk}.
% see my detailed derivation. T_{ijk}=tau_{ijk}/V_{ijk}.

% (1) Find initial values of tau, nu2, PI.
tau = Vbetans;

%% Compute initial values of the Lasso parameters from MLE estimates.
max_nu2=1; % If nu2 are too large, then the Expontial prior variance of lam,phi,psi will be too small, which will result in Lindley's paradox.
% From another point of view, after integrating the scales, we go DE
% prior with scaling parameter 1/nu, nu has to be small enough to make sure
% a heavy tail disttribution.
extra_scale=10; % to make a heavier tail, I can add a extra scale to scale the nu2.
meanop=uneqkron(wavespecs.Kj)*diag(1./wavespecs.Kj);
nu2=min(2./(Vbetans*meanop)/extra_scale,max_nu2);


J_temp=1:wavespecs.J;
pstart0=rep(repmat(0.8*2.^(-(J_temp-1)),model.p,1),wavespecs.Kj);
Zeta=betans./sqrt(Vbetans);
tau0=Vbetans.*max((Zeta).^2-1,1e-6);
BigTau=tau0./Vbetans; 
O=min(pstart0./(1-pstart0).*(1+BigTau).^(-.5).*exp(Zeta.^2/2.*(BigTau./(BigTau+1))),MCMCspecs.maxO);
alpha=O./(O+1);
%meanop=uneqkron(wavespecs.Kj)*diag(1./wavespecs.Kj);
PI=alpha*meanop; 

%% Compute the initial values for pi_{ij} and Hyperprior of nu2_(ij)
[a_nu2,b_nu2,]=Gamma_param_vec(nu2,repmat(MCMCspecs.nu_var,size(nu2,1),size(nu2,2)),0);
[a_pi,b_pi]=Beta_para(PI,MCMCspecs.pi_var,0,0);
