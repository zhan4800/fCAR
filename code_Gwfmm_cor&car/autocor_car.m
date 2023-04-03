function Rivmat=autocor_car(auto,wavespecs)
% This function transforms the parameters of CAR models to the autocorrelation matrix R.
% Input: auto a structure that contains the following:
%        auto.rho: scalor or vector containing the parameter of a proper CAR model at each wavelet jk.
%        auto.type: the type of the functional CAR model, 'separable',
%        'change.over.j', or 'change.over.jk'
%        auto.V: spatial/neighborhood matrix
%        auto.Dv: diagonal matrix of sum of neighbors
%        auto.prior: (a,b) values for the beta prior distribution of rho
%        wavespecs: structure giving the wavelet transform specification
% Output:
% Rivmat: array of CAR correlation marices corresponding to each (or all) wavelet jk.

J=length(auto.rho); % the number of independent CAR models
N=size(auto.V,2); 
K=wavespecs.K;
rho=auto.rho;

switch auto.type
    case 'independent'
        Rivmat = repmat(eye(N),[1,1,K]);
    case 'separable'
        temp = auto.Dv - rho*auto.V;
        Rivmat = repmat(temp,[1,1,K]);
    case 'change.over.jk'
        Rivmat=NaN(N,N,K);
        for j = 1:K
            Rivmat(:,:,j) = auto.Dv - rho(j)*auto.V;
        end
    case'change.over.j'
        Rivmat=NaN(N,N,K);
        waveblk=[0;cumsum(wavespecs.Kj)];
        for j = 1:J
            temp = auto.Dv - rho(j)*auto.V;
            Rivmat(:,:,(waveblk(j)+1):waveblk(j+1)) = repmat(temp,[1,1,wavespecs.Kj(j)]);
        end
end            
        
end

        

% 
% switch autoc.type
%     
%     case 'AR1'
%         autoc.rho
%         
%     case 'Mixed1'
%         err=laprnd(n,wavespecs.K,zeros(n,wavespecs.K), repmat(sqrt(S0)',n,1));
%         u=laprnd(M,wavespecs.K,zeros(M,wavespecs.K),repmat(sqrt(Q0)',M,1));
%     case 'Mixedp'
%         nu=3;
%         err=trnd(nu,n,wavespecs.K).*repmat(sqrt(S0)'/sqrt(nu/(nu-2)),n,1);
%         u=trnd(nu,M,wavespecs.K).*repmat(sqrt(Q0)'/sqrt(nu/(nu-2)),M,1);
%     case 'DampExpo'
%         nu=2;
%         err=trnd(nu,n,wavespecs.K).*repmat(sqrt(S0)'/2,n,1);
%         u=trnd(nu,M,wavespecs.K).*repmat(sqrt(Q0)'/2,M,1);    
% end


