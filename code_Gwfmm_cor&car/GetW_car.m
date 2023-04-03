function W=GetW_car(model,D,Rivmat)
% This function get the W matrix for the case that correlation of E_{jk}^*
% is R.
%%%%% W=GetW(model,D)
%%%%%   Compute cross-products matrices
%%%%%           X'RinvX       X'RinvZ       X'RinvD
%%%%%           Z'RinvX       Z'RinvZ       Z'RinvD
%%%%%           D'RinvX       D'RinvZ       D'RinvD
% Input:    
%
%           model = structure with integer elements:
%               X = (n x p) design matrix for fixed effects functions;
%               Z = cell array containg H matrices (each n x m_h); design
%                       matrices for random effects functions; 
%           D = (n x K) matrix of wavelet coefficients
%           Rinv is a 3D array of CAR correlation marices 
%
% Output:   W = structure with elements:
%           XtX,XtZ,XtD,ZtD,DtD

H=model.H;
    if (H>0)
        Z=model.Z{1};
        if (H>1)
            for h=2:H
                Z=[Z,model.Z{h}];    %#ok<*AGROW> %% concatenate columns of model.Z{h}, h=1,...,H to form single Z matrix.
            end
        end
    else
        Z=0;       %#zhu#% may be make Z to be zeros(n,1) will make matrix product straight forward. 
    end
    
% if size(Rinv,3)==1; % then it is a seperable functinal CAR model
%     fprintf('\n covariance matrix R is constant over j and k. \n');
%     
%     W.XtX=model.X'*Rinv*model.X;
%     W.XtD=model.X'*Rinv*D;
%     W.DtD=D'*Rinv*D;
%     
%     W.ZtZ=Z'*Rinv*Z;
% 
%     if Z==0
%         W.ZtD=zeros(1,size(D,2)); 
%         W.XtZ=zeros(size(model.X,2),1);  
%     else
%         W.ZtD=Z'*Rinv*D;
%         W.XtZ=model.X'*Rinv*Z;
%     end
% 
% else % when covariance matrix R changes over j and/or k
    p = size(model.X,2);
    K = size(D,2);
    
    W.XtX = NaN(p,p,K);
    W.XtD = NaN(p,K);
    W.DtD = NaN(1,K);
    
    for j=1:K
        W.XtX(:,:,j) = model.X'*Rivmat(:,:,j)*model.X;
        W.XtD(:,j) = model.X'*Rivmat(:,:,j)*D(:,j);
        W.DtD(j) = D(:,j)'*Rivmat(:,:,j)*D(:,j);
    end
    
    if (H==0)
        W.ZtD = zeros(1,K);
        W.XtZ = zeros(p,1,K);
        W.ZtZ = zeros(1,1,K); 
    else
        h = size(Z,2);
        W.ZtD = NaN(h,K);
        W.XtZ = NaN(p,h,K);
        W.ZtZ = NaN(h,h,K);
        for j=1:K
            W.ZtD(:,j) = Z'*Rivmat(:,:,j)*D(:,j);
            W.XtZ(:,:,j) = model.X'*Rivmat(:,:,j)*Z;
            W.ZtZ(:,:,j) = Z'*Rivmat(:,:,j)*Z;
        end
    end
    
% end

end
        
    
    
        

