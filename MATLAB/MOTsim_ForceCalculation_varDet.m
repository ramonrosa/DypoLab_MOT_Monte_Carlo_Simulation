clear
clc
close all


detlist = linspace(0,-150,31)';

ravg = zeros(size(detlist,1),3);
rstd = zeros(size(detlist,1),3);

fprintf('\n\n');
for k=1:numel(detlist);
    [ravg(k,:),rstd(k,:)] = MOToptim_ForceCalculation_2(detlist(k));
    fprintf('%d/%d\n',k,numel(detlist));
end

%%

xavg = ravg(:,1)*1000;
yavg = ravg(:,2)*1000;
zavg = ravg(:,3)*1000;

xstd = rstd(:,1)*1000;
ystd = rstd(:,2)*1000;
zstd = rstd(:,3)*1000;

% subplot(1,3,1);
% errorbar(detlist,xavg,xstd,'Marker','o','MarkerFaceColor',[1 0.3 0.3],'MarkerEdgeColor',[0 0 0],'color',[0 0 0],'linestyle','-');
% set(gca,'XDir','Reverse');
% xlabel('Detuning (units of Gamma)');
% ylabel('<x> (mm)');
% 
% subplot(1,3,2);
% errorbar(detlist,yavg,ystd,'Marker','o','MarkerFaceColor',[1 0.3 0.3],'MarkerEdgeColor',[0 0 0],'color',[0 0 0],'linestyle','-');
% set(gca,'XDir','Reverse');
% xlabel('Detuning (units of Gamma)');
% ylabel('<y> (mm)');
% 
% subplot(1,3,3);
errorbar(detlist,zavg,zstd,'Marker','o','MarkerFaceColor',[1 0.3 0.3],'MarkerEdgeColor',[0 0 0],'color',[0 0 0],'linestyle','-');
set(gca,'XDir','Reverse');
xlabel('Detuning (units of Gamma)');
ylabel('<z> (mm)');