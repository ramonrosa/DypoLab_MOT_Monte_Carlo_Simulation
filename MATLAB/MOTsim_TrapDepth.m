clear
clc
% close all

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


%% Read file

[FILE,PATH] = uigetfile('*.dat','Select file','G:\Meu Drive\PosDoc\2018\Simulacoes\MOT_Simulation_v1\Results\');
FILENAME = strcat(PATH,FILE);
% FILENAME = 'G:\Meu Drive\PosDoc\2018\Simulacoes\MOT_Simulation_v1\Results\Test_v1.2_TrapDepth_Dy6B_20180611_133347\TrapDepth.dat';
DATA = importdata(FILENAME);

%%

detuning = DATA(2:end,1);
v = DATA(1,2:end);
D = DATA(2:end,2:end)';
T = linspace(0,0.05,1000);

TE = zeros(numel(T),numel(detuning));
TE(1,:) = D(1,:);

for TEc = 2:numel(T)
    f = ((v/sqrt(kb*T(TEc)/m)).^2).*exp( -((v/sqrt(kb*T(TEc)/m)).^2)/2 );
    f = (f/sum(f))';
    TE(TEc,:) = sum(bsxfun(@times,f,D),1);
end

%%

DTrap = zeros(size(detuning));
for detc = 1:numel(detuning)
    idtrap = find(TE(:,detc)>exp(-1),1,'last');
    if(numel(idtrap)==0)
        idtrap = 1;
    end
    DTrap(detc) = T(idtrap);
end


%%
figure('position',[388,367,1233,420]);
subplot(1,2,1)
    imagesc(detuning,v,D);
    set(gca,'Ydir','Normal');
    set(gca,'Xdir','Reverse');
    xlabel('Detuning (units of \Gamma)');
    ylabel('Initial velocity (m/s)');
    cbh = colorbar;
    ylabel(cbh,'Trapped atoms ratio');
subplot(1,2,2);
    imagesc(detuning,T,TE);
    set(gca,'Ydir','Normal');
    set(gca,'Xdir','Reverse');
    xlabel('Detuning (units of \Gamma)');
    ylabel('Initial temperature (K)');
    cbh = colorbar;
    ylabel(cbh,'Trapped atoms ratio');
    hold on;
    plot(detuning,DTrap,'Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],'Color',[1 0 0],'Linewidth',2);
    hold off;