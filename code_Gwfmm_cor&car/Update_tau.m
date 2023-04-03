function tau=Update_tau(beta,gamma,a_tau,b_tau,meanop,Vbetans,minT)
% This function update tau_{ij} from inverse gamma distribution.

tau=max(1./gamrnd(a_tau+gamma*meanop/2,1./(b_tau+(gamma.*beta.^2)*meanop/2)),Vbetans*meanop*minT);