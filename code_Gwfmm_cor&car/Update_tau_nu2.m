function [tau,nu2]=Update_tau_nu2(beta,gamma,tau,nu2,a_nu2,b_nu2,minVC,wavespecs)
% This function update tau_{ij} from inverse gaussian distribution and nu2
% from gamma distribution.

K=wavespecs.K;
nu2_1=repcol(nu2,wavespecs.Kj);

gam=logical(gamma);
tau(gam)=max(1./invgaussrnd(sqrt(nu2_1(gam)./beta(gam).^2),nu2_1(gam)),minVC);
% constrain tau(or psi) that it is greater than minVC=1e-10, which prevents
% it from being exact 0.


sumop=uneqkron(wavespecs.Kj);
nu2=gamrnd(gamma*sumop+a_nu2,1./((tau.*gamma)*sumop/2+b_nu2)); 

% the following code puts an upper bound on nu2 to lessen the laplace
% shrinkage. 
% meanop = uneqkron(wavespecs.Kj)*diag(1./wavespecs.Kj);
% a = repmat(2./(theta(1,:)*meanop),5,1);
% nu2=min(gamrnd(gamma*sumop+a_nu2,1./((tau.*gamma)*sumop/2+b_nu2)),a); 

