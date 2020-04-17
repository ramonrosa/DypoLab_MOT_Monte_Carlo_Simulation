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

%%

h = figure();
subplot(1,2,1);
% for tw_count=1:1:(numel(tw)-1)
for tw_count=1:10
%     bar(vb,D(:,tw_count),1,'EdgeColor',[0 0 0],'FaceColor',[0.5 0.5 1.0],'LineWidth',1);
    y = D(:,tw_count);
    mask = ones(11,1);
    mask = mask/sum(mask(:));
    yc = conv(y,mask,'same')
%     yc = yc/sum(yc(:));
    yc = yc - (tw_count-1)*0.008;
    plot(vb,yc,'color',[0 0 0],'linewidth',2);
    hold on;
    text(0.3,-(tw_count-1)*0.008,sprintf(' %.0f ms',tw(tw_count)*1000),'HorizontalAlignment','left','VerticalAlignment','middle');
end
hold off;
xlabel('Velocity (m/s)');
ylabel('Normalized frequency');
set(gca,'Ytick',[]);
axis tight;
xlim([0 0.29]);
