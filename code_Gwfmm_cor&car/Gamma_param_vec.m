function [alpha1,beta1,mean1]=Gamma_param_vec(mode1,var,plt)
% according to the mode and variance, compute the two parameters of gamma.
% both mode1 and var are vectors of same size.
alpha1=1+mode1.^2./(2*var)+mode1./(2*var).*sqrt(4*var+mode1.^2);
beta1=(alpha1-1)./mode1;
mean1=beta1./alpha1;

if plt==1
    figure()
    fprintf('\n Now start plotting gamma densities.\n')
    
    for i=1:100    
        idx=input('\n Please input the index of the matrix that you want to plot gamma density with, in form of n\n');
        arg=linspace(min(1e-10,mode1(idx(1),idx(2))),mode1(idx(1),idx(2))+10,200);
        dens=NaN(length(arg),1);
        for k=1:length(arg)
            dens(k)=gampdf(arg(k),alpha1(idx(1),idx(2)),1/beta1(idx(1),idx(2)));
        end   
        plot(arg,dens)
        hold on
        plot([mean1(idx(1),idx(2)),mean1(idx(1),idx(2))],[0,max(dens)*1.1],'-k');
        plot([mode1(idx(1),idx(2)),mode1(idx(1),idx(2))],[0,max(dens)*1.1],'--k');
        hold off
    
        st=input('\n Do you want to stop plotting densities? [y/n]\n','s');
        
        if strmatch(st,'y')==1
            break
        end        
    end
end