function [auto,Vbetans,Wv, W, Riv,nreject]=UpdateR_car(auto,beta,theta,Vbetans,Wv,W, D,model,wavespecs,MCMCspecs, Riv)
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
[n,K] = size(D);

switch auto.type
    case 'independent'
        nreject = 0;
        
    case 'separable' % then there is a single rho to update
        auto_new=auto;
        % calculate the current neg-loglikelihood
        LL = -sum(.5*(Wv_old.L1+Wv_old.L2))-n*log(2*pi)*K/2;
        % sample an auxilary variable u from Unif(0,f(rho))
        logU = log(rand(1))+LL;
        if any(isnan(logU)==1)||any(isinf(logU)) 
            error('Update rho loglikelihood contains Inf or NaN');
        end
        accept = 0; nreject = 0;   % accept: binary indicating if the proposed value is acceptted
                                  % nreject: record the number of rejection
        
        while(accept == 0)
            
        % propose a new value
        rho_prop = betarnd(auto.prior(1),auto.prior(2));   % draw a sample from beta prior
        auto_new.rho = rho_prop;
        Rivnew=autocor_car(auto_new,wavespecs);
        
        % update W and Wv to the new values using the proposed rho.
        Wnew = GetW_car(model, D, Rivnew);
        [~,Vbetans_new,Wv_new]=GetGCP_car(theta, Wnew, Rivnew, model,wavespecs,1,D,MCMCspecs);
        Wv_new.L2=Get_L2(beta,Wv_new,wavespecs);
        LL_new = -sum(.5*(Wv_new.L1+Wv_new.L2))-n*log(2*pi)*K/2;
        if (LL_new > logU)
            accept = 1;
            auto = auto_new;
            Vbetans = Vbetans_new;
            W = Wnew;
            Wv = Wv_new;
            Riv = Rivnew;
        else
            nreject = nreject+1;
        end
        
        end  % end of while statement

        
    case 'change.over.jk'
        K=wavespecs.K;
        nreject = 0;  % nreject: record the number of rejection
        accept = zeros(K,1); % accept: binary indicating if the proposed values are acceptted for all wavelets
         
        auto_new=auto;
        % calculate the current neg-loglikelihood
        LL = -.5*(Wv_old.L1+Wv_old.L2)-n*log(2*pi);
        % sample an auxilary variable u from Unif(0,f(rho)) for ALL wavelets
        logU = log(rand(1,K))+LL; % or generate logu = log(rand(1,K))-LL
            
        while any(accept==0)  % check if all wavelets have been updated
            
            whichtoupdate = find(accept==0); % focus only on those not yet updated wavelets
            
            % draw samples from beta prior for those not yet updated wavelets
            auto_new.rho(whichtoupdate) = betarnd(auto.prior(1),auto.prior(2),length(whichtoupdate),1);   
            Rivnew=autocor_car(auto_new,wavespecs);
            
            % update W and Wv to the new values using the proposed rho.
            Wnew = GetW_car(model, D, Rivnew);
            [~,Vbetans_new,Wv_new]=GetGCP_car(theta, Wnew, Rivnew, model,wavespecs,1,D,MCMCspecs);
            Wv_new.L2=Get_L2(beta,Wv_new,wavespecs);
            LL_new = -.5*(Wv_new.L1+Wv_new.L2)-n*log(2*pi);
                        
            for j = whichtoupdate
                if (LL_new(j) > logU(j))
                    accept(j) = 1;
                    auto.rho(j) = auto_new.rho(j);
                    Vbetans(j) = Vbetans_new(j);

                    W.XtX(:,:,j) = Wnew.XtX(:,:,j);
                    W.XtD(:,j) = Wnew.XtD(:,j);
                    W.Dtd(j,j) = Wnew.Dtd(j,j);
                    W.DtD(j) = Wnew.DtD(j);
                    W.ZtD(:,j) = Wnew.ZtD(:,j);
                    W.XtZ(:,:,j) = Wnew.XtZ(:,:,j);
                    W.ZtZ(:,:,j) = Wnew.ZtZ(:,:,j);
                
                    Wv.XvX(:,:,j) = Wv_new.XvX(:,:,j);
                    Wv.Xvd(:,j) = Wv_new.Xvd(:,j);
                    Wv.XvZ(:,:,j) = Wv_new.XvZ(:,:,j);
                    Wv.ZvZ(:,:,j) = Wv_new.ZvZ(:,:,j);
                    Wv.dvd(j) = Wv_new.dvd(j);
                    Wv.Zvd(:,j) = Wv_new.Zvd(:,j);
                    Wv.L1(j) = Wv_new.L1(j);
                    Wv.L2(j) = Wv_new.L2(j);

                    Riv(:,:,j) = Rivnew(:,:,j);
                end        
            end  % end of j wavelet loop 
            
            nreject=nreject+1; 
        
        end  % end of while loop
        
        
    case'change.over.j'
        J=wavespecs.J;
        sumblk_matrix = uneqkron(wavespecs.Kj);  % KxJ matrix for sumation
        accept = zeros(J,1); % accept: binary indicating if the proposed values are acceptted for all wavelets
         
        auto_new=auto;
        % calculate the current neg-loglikelihood
        LL = (-.5*(Wv_old.L1+Wv_old.L2)-n*log(2*pi))*sumblk_matrix;
        % sample an auxilary variable u from Unif(0,f(rho)) for ALL wavelets
        logU = log(rand(1,J))+LL; % or generate logu = log(rand(1,K))-LL
        
        while any(accept==0)  % check if all wavelets have been updated
            
            whichtoupdate = find(accept==0); % focus only on those not yet updated wavelets
            
            % draw samples from beta prior for those not yet updated wavelets
            auto_new.rho(whichtoupdate) = betarnd(auto.prior(1),auto.prior(2),length(whichtoupdate),1);   
            Rivnew=autocor_car(auto_new,wavespecs);
            
            % update W and Wv to the new values using the proposed rho.
            Wnew = GetW_car(model, D, Rivnew);
            [~,Vbetans_new,Wv_new]=GetGCP_car(theta, Wnew, Rivnew, model,wavespecs,1,D,MCMCspecs);
            Wv_new.L2=Get_L2(beta,Wv_new,wavespecs);
            LL_new = (-.5*(Wv_new.L1+Wv_new.L2)-n*log(2*pi))*sumblk_matrix;
                        
            for j = whichtoupdate
                if (LL_new(j) > logU(j))
                    accept(j) = 1;
                    auto.rho(j) = auto_new.rho(j);
                    kj = find(sumblk_matrix(:,j)==1);
                    Vbetans(kj) = Vbetans_new(kj);
                
                    W.XtX(:,:,kj) = Wnew.XtX(:,:,kj);
                    W.XtD(:,kj) = Wnew.XtD(:,kj);
                	W.Dtd(kj,kj) = Wnew.Dtd(kj,kj);
                	W.DtD(kj) = Wnew.DtD(kj);
                	W.ZtD(:,kj) = Wnew.ZtD(:,kj);
                	W.XtZ(:,:,kj) = Wnew.XtZ(:,:,kj);
                    W.ZtZ(:,:,kj) = Wnew.ZtZ(:,:,kj);
                
                    Wv.XvX(:,:,kj) = Wv_new.XvX(:,:,kj);
                    Wv.Xvd(:,kj) = Wv_new.Xvd(:,kj);
                    Wv.XvZ(:,:,kj) = Wv_new.XvZ(:,:,kj);
                    Wv.ZvZ(:,:,kj) = Wv_new.ZvZ(:,:,kj);
                    Wv.dvd(kj) = Wv_new.dvd(kj);
                	Wv.Zvd(:,kj) = Wv_new.Zvd(:,kj);
                    Wv.L1(kj) = Wv_new.L1(kj);
                    Wv.L2(kj) = Wv_new.L2(kj);

                    Riv(:,:,kj) = Rivnew(:,:,kj);
                end        
            end  % end of j wavelet loop 
            
            nreject=nreject+1; 
        
        end  % end of while loop
        
        
end    % end of switch/case statement


