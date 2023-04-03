% this script test for autocorrelation matrix.
auto.type='AR1';
auto.rho=0.6;
auto.Ti=cell(3,1);
auto.Ti{1}=linspace(0,1,3)';
auto.Ti{2}=linspace(0,1,4)';
auto.Ti{3}=linspace(0,1,5)';
%[Rmat,Lmat,R1,L1]=autocor(auto);
[Rmat,Rivmat, Lmat,R1,L1,Rinv]=autocor(auto);
% Mixed.p
auto.type='Mixed.p';
auto.rho=cell(3,1);
auto.rho{1}=[0.6,0.5];
auto.rho{2}=0.3;
auto.rho{3}=[0.4,0.3,0.2];
auto.Ti=cell(3,1);
auto.Ti{1}=linspace(0,1,3)';
auto.Ti{2}=linspace(0,1,4)';
auto.Ti{3}=linspace(0,1,5)';
%[Rmat,Lmat,R1,L1]=autocor(auto);
[Rmat,Rivmat, Lmat,R1,L1,Rinv]=autocor(auto);

% DampedExponential
auto.type='DampedExponential';
auto.ts=0.5;
auto.rho=0.6;
auto.Ti=cell(3,1);
auto.Ti{1}=linspace(0,1,3)';
auto.Ti{2}=linspace(0,1,4)';
auto.Ti{3}=linspace(0,1,5)';
%[Rmat,Lmat,R1,L1]=autocor(auto);
[Rmat,Rivmat, Lmat,R1,L1,Rinv]=autocor(auto);

