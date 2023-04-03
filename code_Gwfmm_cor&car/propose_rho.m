function rho_new = propose_rho(rho, rhosd)
% This function generate rho_new using random
% walk proposal, using logit transform.
% i.e., logit(rho_new)=N(logit(rho),rho_sd.^2)

logit_rhonew=log(rho)-log(1-rho)+randn(1,length(rho)).*rhosd;
rho_new=exp(logit_rhonew)./(1+exp(logit_rhonew));

