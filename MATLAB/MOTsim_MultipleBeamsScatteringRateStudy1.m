clear
clc

s0 = 160;
N = 40;
delta = 20;

Rs0 = s0./(1 + s0 + 4*delta.^2);
sf = s0/N;
Rsf = sf./(1 + sf + 4*delta.^2);


tau0 = -log(1-rand(100000,1))/Rs0;
subplot(1,2,1);
hist(tau0,50);
Rs0_exp = 1./mean(tau0);
fprintf('\n\n    Rs0 = %f\n',Rs0);
fprintf('Rs0_exp = %f\n',Rs0_exp);

tauF_each = -log(1-rand(100000,N))/Rsf;
tauF = min(tauF_each,[],2);
tauF = min(tauF_each,[],2);
subplot(1,2,2);
hist(tauF,50);
Rsf_exp = 1./mean(tauF);
fprintf('    Rsf = %f\n',Rsf);
fprintf('Rsf_exp = %f\n\n\n',Rsf_exp);

