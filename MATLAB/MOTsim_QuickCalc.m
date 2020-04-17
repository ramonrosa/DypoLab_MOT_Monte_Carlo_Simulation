clear
clc
close all

%% Constantes

h = 6.62607004e-34;
mub = 9.274009994e-24;
e = 1.60217662e-19;
kb = 1.38064852e-23;
c = 2.99792458e8;

%% Atomo

m = (1.660539040e-27)*163.9291748;
Gamma = 136e3;
lambda = 626e-9;
Jgnd = 8;
Jexc = 9;
gjg = 1.24;
gje = 1.29;

%%

T=10e-6;


%% Import data

[FILE,PATH] = uigetfile('*.dat','Select file','G:\Meu Drive\PosDoc\2018\Simulacoes\MOT_Simulation_v1\Results\');
FILENAME = strcat(PATH,FILE);
% FILENAME = 'G:\Meu Drive\PosDoc\2018\Simulacoes\MOT_Simulation_v1\Results\Test_v1.2_TrapDepth_Dy6B_20180611_133347\TrapDepth.dat';
D = importdata(FILENAME);

%%



vb = D(:,1);
vc = D(:,2);

F = @(a) a(1)*((vb/a(2)).^2).*exp(-(vb/a(2)).^2/2);
E = @(a) sum((vc - F(a)).^2);
options = optimset('TolX',1e-10,'TolFun',1e-10);
[af,~,exitflag,output] = fminsearch(E,[max(vc),1],options);

TEMPERATURE = af(2)^2*m/kb;
[TEMPNUM,TEMPPREF] = Num2Sci(TEMPERATURE);
TEMPSTR = sprintf('T_{M-B} = %.2f %sK',TEMPNUM,TEMPPREF);

bar(vb,vc,1,'EdgeColor',[0 0 0],'FaceColor',[0.5 0.5 1.0],'LineWidth',1);
hold on;
plot(vb,F(af),'color',[1 0 0],'linewidth',2);
title(TEMPSTR);
xlim([0 vb(end)]);
ylim([0 max([F(af);vc])*1.05]);
hold off;
xlabel('Velocity (m/s)');
ylabel('Normalized frequency');
set(gca,'Ytick',[]);