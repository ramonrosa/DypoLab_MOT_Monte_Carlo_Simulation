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
DATA = importdata(FILENAME);

%%

vb = DATA(2:end,1);
tw = DATA(1,2:end);
D = DATA(2:end,2:end);

FM = zeros(size(D));
TEMPERATURE = zeros(1,size(D,2));

for tw_count=1:numel(tw)
    vc = D(:,tw_count);

    F = @(a) a(1)*((vb/a(2)).^2).*exp(-(vb/a(2)).^2/2);
    E = @(a) sum((vc - F(a)).^2);
    options = optimset('TolX',1e-10,'TolFun',1e-10,'MaxIter',1e5,'MaxFunEvals',1e6);
    [af,~,exitflag,output] = fminsearch(E,[max(vc),1],options);

    TEMPERATURE(tw_count) = (af(2)^2)*m/kb;
    FM(:,tw_count) = F(af);
end

%%

h = figure();
for tw_count=1:numel(tw)
    figure(h);
    bar(vb,D(:,tw_count),1,'EdgeColor',[0 0 0],'FaceColor',[0.5 0.5 1.0],'LineWidth',1);
    hold on;
%     plot(vb,FM(:,tw_count),'color',[1 0 0],'linewidth',2);
    
    [TEMPNUM,TEMPPREF] = Num2Sci(TEMPERATURE(tw_count));
    TEMPSTR = sprintf('T_{M-B} = %.2f %sK %e',TEMPNUM,TEMPPREF);
    
    title(TEMPSTR);
    
    %xlim([0 vb(end)]);
    xlim([0 0.5]);
    ylim([0 max([FM(:,tw_count);D(:,tw_count)])*1.05]);
    hold off;
    xlabel('Velocity (m/s)');
    ylabel('Normalized frequency');
    set(gca,'Ytick',[]);
    drawnow();
    pause(1);
end