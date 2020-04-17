clear
clc
close all

load 'G:\Meu Drive\PosDoc\2018\Simulacoes\MOT_Simulation_v1\Results\Geometry_Optimization_3B\results_mat.mat'

iMAT = interp2(MAT,6)';
idet = linspace(detuning(1),detuning(end),size(iMAT,2))';
iang = linspace(ANGLE(1),ANGLE(end),size(iMAT,1));

imagesc(idet,90-iang,1000*iMAT);
set(gca,'XDir','Reverse');
set(gca,'YDir','Normal');
xlabel ('Detuning (units of \Gamma)');
ylabel ('Angle of beams in respect to z-axis (degrees)');
colorbar;
title('Trap depth (mK)');
hold on;
% contour(idet,90-iang,iMAT*1000,2:2:12,'color',[0 0 0],'ShowText','on','LabelSpacing',1e6);
CH = contour(idet,90-iang,iMAT*1000,3:3:12,'color',[0 0 0]);
clabel(CH,'FontSize',15,'FontWeight','bold','Color','black');
hold off;