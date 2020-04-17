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

ANGLE = 2.5:2.5:60.0;

%%
for k=1:numel(ANGLE)
    PATH = uigetdir('G:\Meu Drive\PosDoc\2018\Simulacoes\MOT_Simulation_v1\Results\ResultsToKeepForGood\Geometry_Optimization_3B',sprintf('Select folder: %f deg',ANGLE(k)));
    FILENAME = strcat(PATH,'\TrapDepth.dat');
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
    cDETUNING{k} = detuning;
    cTRAPDEPTH{k} = DTrap;
end

% for k=1:5
%     cTRAPDEPTH{k} = [cTRAPDEPTH{k};zeros(10,1)];
% end

%%

MAT = cell2mat(cTRAPDEPTH);

size(detuning);
imagesc(detuning,90-ANGLE,1000*MAT');
set(gca,'XDir','Reverse');
set(gca,'YDir','Normal');
xlabel ('Detuning (units of \Gamma)');
ylabel ('Angle of beams with {\it xy}-plane (degrees)');
colorbar;
title('Trap depth (mK)');