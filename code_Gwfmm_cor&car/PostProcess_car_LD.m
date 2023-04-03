function [betahat,beta_025CI,beta_975CI,alphahat,accept_rate_theta,ghat,ghatns,Q025_ghat,Q975_ghat,thetahat,rhohat,accept_rate_rho,theta_025CI,theta_975CI,tauhat,pihat,...
    Sigma,uhat,u_025CI,u_975CI,g_Ut,g_Ut025,g_Ut975]=...
    PostProcess_car_LD(paramroute,betans,theta,model,wavespecs,B,get_sigma,sampleU)

p=model.p;
T=wavespecs.T;
K=wavespecs.K;
H=model.H;

% step (1) find posterior mean and quantiles for beta. and that of
% gBeta(t).
fid = fopen(strcat(paramroute,'/MCMC_beta_rec.txt'),'r');
N = 1000;

for iter=1:(B/N)
    MCMC_beta_row=textscan(fid,repmat('%f',1,p*K),N,'Delimiter',',');
    MCMC_beta_row=[MCMC_beta_row{:}];
    MCMC_g_row=NaN(N,p*T);
    for i=1:p     
        gsamples=idwt_rows(MCMC_beta_row(:,((i-1)*K+1):(i*K)),wavespecs);  %%% 18 seconds for B=100
        start_col=(i-1)*T+1;
        end_col=i*T;
        MCMC_g_row(:,start_col:end_col)=gsamples;  % Assemble MCMC samples of data-space betas   
    end
    dlmwrite(strcat(paramroute,'/MCMC_gBetat_rec.txt'),MCMC_g_row,'-append','precision','%12.6e');
    fprintf('%d \n',iter);
end
clear MCMC_g_row;

betatemp = NaN(3,p*K);
gtemp = NaN(3,p*T);
N = 20000;

for iter = 1:floor(K*p/N)
    fid = fopen(strcat(paramroute,'/MCMC_beta_rec.txt'),'r');
    fid2 = fopen(strcat(paramroute,'/MCMC_gBetat_rec.txt'),'r');
    
    start_col = N*(iter-1)+1;
    end_col = N*iter;

    MCMC_beta_col=textscan(fid,[repmat('%*f',1,start_col-1),repmat('%f',1,N),'%*[^\n]'],'delimiter',',');
    MCMC_beta_col=[MCMC_beta_col{:}];
    betatemp(1,start_col:end_col)=mean(MCMC_beta_col);
    betatemp(2:3,start_col:end_col)=quantile(MCMC_beta_col,[0.025,0.975]);

    MCMC_g_col = textscan(fid2,[repmat('%*f',1,start_col-1),repmat('%f',1,N),'%*[^\n]'],'delimiter',',');
    MCMC_g_col=[MCMC_g_col{:}];
    gtemp(1,start_col:end_col)=mean(MCMC_g_col);
    gtemp(2:3,start_col:end_col)=quantile(MCMC_g_col,[0.025,0.975]);
end

if mod(K*p,N)~=0
    fid = fopen(strcat(paramroute,'/MCMC_beta_rec.txt'),'r');
    fid2 = fopen(strcat(paramroute,'/MCMC_gBetat_rec.txt'),'r');
    
    start_col = N*iter+1;
    end_col = K*p;
    
    MCMC_beta_col=textscan(fid,[repmat('%*f',1,start_col-1),repmat('%f',1,(end_col-start_col+1))],'delimiter',',');
    MCMC_beta_col=[MCMC_beta_col{:}];
    betatemp(1,start_col:end_col)=mean(MCMC_beta_col);
    betatemp(2:3,start_col:end_col)=quantile(MCMC_beta_col,[0.025,0.975]);

    end_col = T*p;
    MCMC_g_col = textscan(fid2,[repmat('%*f',1,start_col-1),repmat('%f',1,(end_col-start_col+1))],'delimiter',',');
    MCMC_g_col=[MCMC_g_col{:}];
    gtemp(1,start_col:end_col)=mean(MCMC_g_col);
    gtemp(2:3,start_col:end_col)=quantile(MCMC_g_col,[0.025,0.975]);
end

betahat = reshape(betatemp(1,:),K,p)';
beta_025CI = reshape(betatemp(2,:),K,p)';
beta_975CI = reshape(betatemp(3,:),K,p)';

ghat = reshape(gtemp(1,:),T,p)';
Q025_ghat = reshape(gtemp(2,:),T,p)';
Q975_ghat = reshape(gtemp(3,:),T,p)';

fprintf('\n Done with function coefficient \n');
clear MCMC_beta_col MCMC_g_col;
fclose(fid); fclose(fid2);

ghatns=idwt_rows(betans,wavespecs);%#zhu#% note that betans is not changed during MCMC, it is only an initial values.


% step (2) find posterior mean and quantiles for alpha,acpt_theta,tau,pi. 
fid = fopen(strcat(paramroute,'/MCMC_alpha_rec.txt'));
N = 100;
alphahat_temp = NaN(B/N,K*p);

for iter=1:(B/N)
    MCMC_alpha_row=textscan(fid,repmat('%f',1,p*K),N,'Delimiter',',');
    MCMC_alpha_row=[MCMC_alpha_row{:}];
    alphahat_temp(iter,:)=mean(MCMC_alpha_row);
end
alphahat=reshape(mean(alphahat_temp)',K,p)';
clear MCMC_alpha_row alphahat_temp;

fid = fopen(strcat(paramroute,'/MCMC_flag_theta_rec.txt'));
accept_rate_theta_temp = NaN(B/N,K*size(theta,1));
for iter=1:(B/N)
    MCMC_flag_theta_row=textscan(fid,repmat('%f',1,size(theta,1)*K),N,'Delimiter',',');
    MCMC_flag_theta_row=[MCMC_flag_theta_row{:}];
    accept_rate_theta_temp(iter,:)=mean(MCMC_flag_theta_row);
end
accept_rate_theta=reshape(mean(accept_rate_theta_temp),size(theta,1),K)';
clear MCMC_flag_theta_row accept_rate_theta_temp;

% MCMC_flag_R=dlmread(strcat(paramroute,'/MCMC_flag_R_rec.txt'));
% accept_rate_rho=mean(MCMC_flag_R);
% clear MCMC_flag_R;

MCMC_tau=dlmread(strcat(paramroute,'/MCMC_tau_rec.txt'));
tauhat=reshape(mean(MCMC_tau),p,size(MCMC_tau,2)/p);
clear MCMC_tau;

MCMC_pi=dlmread(strcat(paramroute,'/MCMC_pi_rec.txt'));
pihat=reshape(mean(MCMC_pi)',p,wavespecs.J);
clear MCMC_pi;

% step (3) get covariance Sigma of U_h(t) or E_c(t) in time domain, if needed. 
if sampleU==1
    M=model.M;
    %----------------------------------------------------------------------    
    fid = fopen(strcat(paramroute,'/MCMC_U_rec.txt'));
    N=1000;
    %----------------------------------------------------------------------
    
    for iter=1:(B/N)
        MCMC_U_row=textscan(fid,repmat('%f',1,M*K),N,'Delimiter',',');
        MCMC_U_row=[MCMC_U_row{:}];
        MCMC_gUt_row=NaN(N,M*T);
        for j=1:M     
            gUtsamples=idwt_rows(MCMC_U_row(:,((j-1)*K+1):(j*K)),wavespecs);
            start_col=(j-1)*T+1;
            end_col=j*T;
            MCMC_gUt_row(:,start_col:end_col)=gUtsamples;  
        end
        dlmwrite(strcat(paramroute,'/MCMC_gUt_rec.txt'),MCMC_gUt_row,'-append','precision','%12.6e');
    end
    clear MCMC_U_row MCMC_gUt_row;
    
    Utemp = NaN(3,M*K);
    gUttemp = NaN(3,M*K);    
    N = 10000;

    for iter = 1:floor(K*M/N)
        fid = fopen(strcat(paramroute,'/MCMC_U_rec.txt'));
        fid2 = fopen(strcat(paramroute,'/MCMC_gUt_rec.txt'),'r');
        
        start_col = N*(iter-1)+1;
        end_col = N*iter;

        MCMC_U_col=textscan(fid,[repmat('%*f',1,start_col-1),repmat('%f',1,N),'%*[^\n]'],'delimiter',',');
        MCMC_U_col=[MCMC_U_col{:}];
        Utemp(1,start_col:end_col)=mean(MCMC_U_col);
        Utemp(2:3,start_col:end_col)=quantile(MCMC_U_col,[0.025,0.975]);

        MCMC_gUt_col = textscan(fid2,[repmat('%*f',1,start_col-1),repmat('%f',1,N),'%*[^\n]'],'delimiter',',');
        MCMC_gUt_col=[MCMC_gUt_col{:}];
        gUttemp(1,start_col:end_col)=mean(MCMC_gUt_col);
        gUttemp(2:3,start_col:end_col)=quantile(MCMC_gUt_col,[0.025,0.975]);
    end

    if mod(K*M,N)~=0
        fid = fopen(strcat(paramroute,'/MCMC_U_rec.txt'));
        fid2 = fopen(strcat(paramroute,'/MCMC_gUt_rec.txt'),'r');

        start_col = N*iter+1;
        end_col = K*M;
    
        MCMC_U_col=textscan(fid,[repmat('%*f',1,start_col-1),repmat('%f',1,end_col-start_col+1)],'delimiter',',');
        MCMC_U_col=[MCMC_U_col{:}];
        Utemp(1,start_col:end_col)=mean(MCMC_U_col);
        Utemp(2:3,start_col:end_col)=quantile(MCMC_U_col,[0.025,0.975]);

        MCMC_gUt_col = textscan(fid2,[repmat('%*f',1,start_col-1),repmat('%f',1,end_col-start_col+1)],'delimiter',',');
        MCMC_gUt_col=[MCMC_gUt_col{:}];
        gUttemp(1,start_col:end_col)=mean(MCMC_gUt_col);
        gUttemp(2:3,start_col:end_col)=quantile(MCMC_gUt_col,[0.025,0.975]);
    end

    uhat = reshape(Utemp(1,:),K,M)';
    u_025CI = reshape(Utemp(2,:),K,M)';
    u_975CI = reshape(Utemp(3,:),K,M)';

    g_Ut = reshape(gUttemp(1,:),T,M)';
    g_Ut025 = reshape(gUttemp(2,:),T,M)';
    g_Ut975 = reshape(gUttemp(3,:),T,M)';

    clear MCMC_U_col MCMC_gUt_col;
else
    uhat=[];
    u_025CI=[];
    u_975CI=[];
    g_Ut=[];
    g_Ut025=[];
    g_Ut975=[];
end

fid = fopen(strcat(paramroute,'/MCMC_rho_rec.txt'));
N = 100;
rhohat_temp = NaN(B/N,K);
for iter=1:(B/N)
    MCMC_rho_row=textscan(fid,repmat('%f',1,K),N,'Delimiter',',');
    MCMC_rho_row=[MCMC_rho_row{:}];
    rhohat_temp(iter,:)=mean(MCMC_rho_row);
end
rhohat = mean(rhohat_temp);
clear MCMC_rho_row rhohat_temp;

fid = fopen(strcat(paramroute,'/MCMC_nreject_R_rec.txt'));
accept_rate_rho_temp = NaN(B/N,K);
for iter=1:(B/N)
    MCMC_flag_R_row=textscan(fid,repmat('%f',1,K),N,'Delimiter',',');
    MCMC_flag_R_row=[MCMC_flag_R_row{:}];
    accept_rate_rho_temp(iter,:)=mean(MCMC_flag_R_row);
end
accept_rate_rho = mean(accept_rate_rho_temp);
clear MCMC_flag_R_row accept_rate_rho_temp;


thetatemp = NaN(3,size(theta,1)*size(theta,2));
N = 20000;

for iter = 1:floor(size(theta,1)*size(theta,2)/N)
    fid = fopen(strcat(paramroute,'/MCMC_theta_rec.txt'));
    start_col = N*(iter-1)+1;
    end_col = N*iter;

    MCMC_theta_col=textscan(fid,[repmat('%*f',1,start_col-1),repmat('%f',1,N),'%*[^\n]'],'delimiter',',');
    MCMC_theta_col=[MCMC_theta_col{:}];
    thetatemp(1,start_col:end_col)=mean(MCMC_theta_col);
    thetatemp(2:3,start_col:end_col)=quantile(MCMC_theta_col,[0.025,0.975]);
end

if mod(size(theta,1)*size(theta,2),N)~=0
    fid = fopen(strcat(paramroute,'/MCMC_theta_rec.txt'));
    start_col = N*iter+1;
    end_col = size(theta,1)*size(theta,2);
    
    MCMC_theta_col=textscan(fid,[repmat('%*f',1,start_col-1),repmat('%f',1,end_col-start_col+1)],'delimiter',',');
    MCMC_theta_col=[MCMC_theta_col{:}];
    thetatemp(1,start_col:end_col)=mean(MCMC_theta_col);
    thetatemp(2:3,start_col:end_col)=quantile(MCMC_theta_col,[0.025,0.975]);
end

thetahat=reshape(thetatemp(1,:),size(theta,1),size(theta,2));
theta_025CI=reshape(thetatemp(2,:),size(theta,1),size(theta,2));
theta_975CI=reshape(thetatemp(3,:),size(theta,1),size(theta,2));
clear MCMC_theta_col thetatemp;

if (get_sigma==1)
    Sigma.td=NaN(T,T,H+model.c); 
    Sigma.diagtd=NaN(T,H+model.c);
    Sigma.Rhotd=NaN(T,T,H+model.c);    
    for h=1:model.H
        [Sigma.td(:,:,h),Sigma.diagtd(:,h),Sigma.Rhotd(:,:,h)]=GetSigma(thetahat(h,:),wavespecs);
    end
    for cc=1:model.c
    [Sigma.td(:,:,model.H+cc),Sigma.diagtd(:,model.H+cc),Sigma.Rhotd(:,:,model.H+cc)]=GetSigma(thetahat(model.H+cc,:),wavespecs);
    end
else
    Sigma=[];    
end

