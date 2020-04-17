clear
clc
close all

% x = v/a;
% a = sqrt(kb*T/m);

% f = Maxwell-Boltzmann distribution
f = @(x) sqrt(2/pi)*(x.^2).*exp(-(x.^2)/2);

%y = CDF(f);
y = @(x) erf(x/sqrt(2)) - sqrt(2/pi)*(x.*exp(-x.^2/2));

options = optimset('TolFun',1e-16,'TolX',1e-16,'MaxIter',1e7,'MaxFunEvals',1e9,'Display','off');

y0set = [0,0.0010,0.0100,0.0500,0.1000,0.2000,0.3000,0.4000,0.5000,0.6000,0.7000,0.8000,0.9000,0.9900,0.9990,0.9999];
x0set = zeros(size(y0set));

for k=2:numel(y0set)
    y0 = y0set(k);
    yr = @(x) y(x) - y0;
    x0set(k) = fsolve(yr,1,options);
end
fprintf('\n\n');
fprintf('%.4f,',y0set);
fprintf('\n\n');
fprintf('%.16f,',x0set);
fprintf('\n\n');
fprintf('%.16f,',f(x0set));
fprintf('\n\n');