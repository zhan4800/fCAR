function [auto,Vbetans,Wv, W, Riv,nreject]=UpdateR_MH_car(auto,beta,theta,Vbetans,Wv,W, D,model,wavespecs,MCMCspecs, Riv, sampleU)
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

switch auto.type
    case 'independent'
        auto.rho = 0;
        nreject=0;
        
    case 'separable' % then there is a single rho to update
        nreject=0;
        
        auto_new=auto;
        rho_prop = propose_rho(auto.rho, auto.rhosd);
        auto_new.rho = rho_prop;
        Rivnew=autocor_car(auto_new,wavespecs);
        
        % update W and Wv to the new values using the proposed rho.
        Wnew = GetW_car(model, D, Rivnew);
        [~,Vbetans_new,Wv_new]=GetGCP_car(theta, Wnew, Rivnew, model,wavespecs,1,D,MCMCspecs);
        Wv_new.L2=Get_L2(beta,Wv_new,wavespecs);

        A1=sum(.5*(Wv_old.L1-Wv_new.L1+Wv_old.L2-Wv_new.L2)); % log likelihood ratio
        A2=(auto.prior(1)-1)'.*(log(rho_prop)-log(auto.rho))+(auto.prior(2)-1)'.*(log(1-rho_prop)-log(1-auto.rho));
        Prop_ratio=log(rho_prop)-log(auto.rho)+log(1-rho_prop)-log(1-auto.rho); 
        log_ratio=A1+A2+Prop_ratio;
        if any(isnan(log_ratio)==1)||any(isinf(log_ratio)) 
            error('Update Theta M-H logratio contains Inf or NaN');
        end

        if (log(rand(1,1))<log_ratio)
            nreject = 1; %accepted.
            auto = auto_new;
            Vbetans = Vbetans_new;
            W = Wnew;
            Wv = Wv_new;
            Riv = Rivnew;
        end

        
    case 'change.over.jk'
        K=wavespecs.K;
        nreject = zeros(1,K);  % nreject: record the number of rejected j (wavelet)
        
        auto_new=auto;
        rho_prop = propose_rho(auto.rho, auto.rhosd);
        auto_new.rho = rho_prop;
        Rivnew=autocor_car(auto_new,wavespecs);
        
        % update W and Wv to the new values using the proposed rho.
        Wnew = GetW_car(model, D, Rivnew);
        [~,Vbetans_new,Wv_new]=GetGCP_car(theta, Wnew, Rivnew, model,wavespecs,1,D,MCMCspecs);
        Wv_new.L2=Get_L2(beta,Wv_new,wavespecs);

        A1=.5*(Wv_old.L1-Wv_new.L1+Wv_old.L2-Wv_new.L2); % log likelihood ratio
        A2=(auto.prior(1)-1)'.*(log(rho_prop)-log(auto.rho))+(auto.prior(2)-1)'.*(log(1-rho_prop)-log(1-auto.rho));
        Prop_ratio=log(rho_prop)-log(auto.rho)+log(1-rho_prop)-log(1-auto.rho); 
        log_ratio=A1+A2+Prop_ratio;
        
        for j=1:K
            if (log(rand(1,1))<log_ratio(j))
                nreject(j) = 1;
                auto.rho(j) = auto_new.rho(j);
                Vbetans(:,j) = Vbetans_new(:,j);
                
                W.XtX(:,:,j) = Wnew.XtX(:,:,j);
                W.XtD(:,j) = Wnew.XtD(:,j);
                %W.Dtd(j,j) = Wnew.Dtd(j,j);
                W.DtD(j) = Wnew.DtD(j);
                
                Wv.XvX(:,:,j) = Wv_new.XvX(:,:,j);
                Wv.Xvd(:,j) = Wv_new.Xvd(:,j);
                Wv.dvd(j) = Wv_new.dvd(j);
                Wv.L1(j) = Wv_new.L1(j);
                Wv.L2(j) = Wv_new.L2(j);

                Riv(:,:,j) = Rivnew(:,:,j);
                
                if (model.H>0)&&(sampleU==1) 
                    W.ZtD(:,j) = Wnew.ZtD(:,j);
                    W.XtZ(:,:,j) = Wnew.XtZ(:,:,j);
                    W.ZtZ(:,:,j) = Wnew.ZtZ(:,:,j);    
                    Wv.XvZ(:,:,j) = Wv_new.XvZ(:,:,j);
                    Wv.ZvZ(:,:,j) = Wv_new.ZvZ(:,:,j);
                    Wv.Zvd(:,j) = Wv_new.Zvd(:,j);
                end

            end        
        end 
        
        
    case'change.over.j'
        J=wavespecs.J;
        sumblk_matrix = uneqkron(wavespecs.Kj);  % KxJ matrix for sumation
        nreject = zeros(1,J); % accept: binary indicating if the proposed values are acceptted for all wavelets
         
        auto_new=auto;
        rho_prop = propose_rho(auto.rho, auto.rhosd);
        auto_new.rho = rho_prop;
        Rivnew=autocor_car(auto_new,wavespecs);
        
        % update W and Wv to the new values using the proposed rho.
        Wnew = GetW_car(model, D, Rivnew);
        [~,Vbetans_new,Wv_new]=GetGCP_car(theta, Wnew, Rivnew, model,wavespecs,1,D,MCMCspecs);
        Wv_new.L2=Get_L2(beta,Wv_new,wavespecs);

        A1=.5*(Wv_old.L1-Wv_new.L1+Wv_old.L2-Wv_new.L2)*sumblk_matrix; % log likelihood ratio
        A2=(auto.prior(1)-1)'.*(log(rho_prop)-log(auto.rho))+(auto.prior(2)-1)'.*(log(1-rho_prop)-log(1-auto.rho));
        Prop_ratio=log(rho_prop)-log(auto.rho)+log(1-rho_prop)-log(1-auto.rho); 
        log_ratio=A1+A2+Prop_ratio;
        
        for j=1:J
            if (log(rand(1,1))<log_ratio(j))
                nreject(j) = 1;
                auto.rho(j) = auto_new.rho(j);
                kj = find(sumblk_matrix(:,j)==1);
                Vbetans(:,kj) = Vbetans_new(:,kj);
                
                W.XtX(:,:,kj) = Wnew.XtX(:,:,kj);
                W.XtD(:,kj) = Wnew.XtD(:,kj);
                %W.Dtd(kj,kj) = Wnew.Dtd(kj,kj);
                W.DtD(kj) = Wnew.DtD(kj);
                
                Wv.XvX(:,:,kj) = Wv_new.XvX(:,:,kj);
                Wv.Xvd(:,kj) = Wv_new.Xvd(:,kj);
                Wv.dvd(kj) = Wv_new.dvd(kj);
                Wv.L1(kj) = Wv_new.L1(kj);
                Wv.L2(kj) = Wv_new.L2(kj);

                Riv(:,:,kj) = Rivnew(:,:,kj);
                
                if (model.H>0)&&(sampleU==1) 
                    W.ZtD(:,kj) = Wnew.ZtD(:,kj);
                    W.XtZ(:,:,kj) = Wnew.XtZ(:,:,kj);
                    W.ZtZ(:,:,kj) = Wnew.ZtZ(:,:,kj);
                    Wv.XvZ(:,:,kj) = Wv_new.XvZ(:,:,kj);
                    Wv.ZvZ(:,:,kj) = Wv_new.ZvZ(:,:,kj);
                    Wv.Zvd(:,kj) = Wv_new.Zvd(:,kj);
                end

            end        
        end
        
        
end    % end of switch/case statement


