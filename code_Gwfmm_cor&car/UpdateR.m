function [flag,auto,Vbetans,Wv, W, Riv]=UpdateR(auto,beta,theta,Vbetans,Wv,W, D,model,wavespecs,MCMCspecs, Riv)
%%%%    UpdateR(model,beta,theta,L1,L2,W,a,b,propsd_theta): Updates the
%%%%    parameters for auto correlations in E_{jk}*, the autho correlations
%%%%    are assumed to be common across (j,k).
%%%%
%%%%    Input:  model = structure containing model information.           
%                   X = (n x p) design matrix for fixed effects functions;
%                   Z = cell array of H matrices, each containing (n x m_h) design
%                       matrix for set of random effects functions; 
%               beta: p x K matrix containing values of fixed effects from
%                       current MCMC iteration
%               theta: H+c x * matrix containing current values of variance
%                       components
%                           (*=1 if model.covQ=0, J if model.covQ=1, K if model.covQ=2)
%               betans: p x K matrix of non-shrinkage estimates of fixed
%                           effect coefficients
%               Vbetans: p x K matrix containing estimated variances of each of the
%                       non-shrinkage fixed effects estimates
%               L1: 1 x * matrix containing values of log|Sigma_jk|
%               L2: 1 x * matrix containing values of 
%                   (d_jk-X beta_jk)'Sigma(Theta_jk)^{-1} (d_jk-X beta_jk)
%                   (summed over k if covQ=1)
%
%               W = structure containing:
%                   XtX=X'X
%                   XtZ=X'Z
%                   XtD=X'D
%                   ZtD=Z'D 
%              auto: a structure that contains current parameters of rho
%              and all other info.
%                   auto.rho = a cell of size number of blocks by 1, each cell
%                     contains the parameters for the particular block.
%                   auto.rhosd = same structure as rho, but specify the proposal
%               standard deviation of the correponding rho parameter. 
%              auto.Ti, auto.num_of_blocks.
%              auto.cd = matrix of length(rho) by 2 containing the Beta priors Beta(c,d) for
              
%       Output:
%               theta = new value of theta
%               newtheta = 1 x * matrix containing indicator if "new"
%                           theta accepted by Metropolis
%               XvX, Xvd, dvd, XvZ, ZvZ, Zvd = updated values of
%                           X'(V^(-1))X, etc
%               betans, Vbetans = updated values of betans, Vbetans
%               L1, L2: new values of L1 and L2
%
%%%%    Functions needed: 
%%%%      GetGCP.m, Get_L2.m, rep.m                  
Wv_old=Wv;

% propose a new value of rho.
rho_new = propose_rho(auto.rho, auto.rhosd);
auto_new=auto;
auto_new.rho=rho_new;
Rivnew=autocor(auto_new);

% update W and Wv to the new values using the proposed rho.
Wnew = GetW_cor1(model, D, Rivnew);
[~,Vbetans_new,Wv_new]=GetGCP_cor1(theta, Wnew, Rivnew, model,wavespecs,1,D,MCMCspecs);
Wv_new.L2=Get_L2(beta,Wv_new,wavespecs);

A1=sum(.5*(Wv_old.L1-Wv_new.L1+Wv_old.L2-Wv_new.L2)); % log likelihood ratio
A2=sum((auto.prior(:,1)-1)'.*(log(rho_new)-log(auto.rho))+(auto.prior(:,2)-1)'.*(log(1-rho_new)-log(1-auto.rho)));
Prop_ratio=sum(log(rho_new)-log(auto.rho)+log(1-rho_new)-log(1-auto.rho)); 
log_ratio=A1+A2+Prop_ratio;
if any(isnan(log_ratio)==1)||any(isinf(log_ratio)) 
    error('Update Theta M-H logratio contains Inf or NaN');
end
flag=log(rand(1,1))<log_ratio; % with probability ratio, accept the new theta.

if (log(rand(1,1))<log_ratio)
    flag = 1; %accepted.
    auto = auto_new;
    Vbetans = Vbetans_new;
    W = Wnew;
    Wv = Wv_new;
    Riv = Rivnew;
end
