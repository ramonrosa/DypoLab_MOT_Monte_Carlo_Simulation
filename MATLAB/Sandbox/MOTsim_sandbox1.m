clear
clc
close all

n01 = (rand(3,1)-0.5)*2;
m01 = rand(1)*10;
n02 = (rand(3,1)-0.5)*2;
m02 = rand(1)*10;

nf1 = (rand(3,1)-0.5)*2;



n01 = n01/norm(n01);
n02 = n02/norm(n02);
nf1 = nf1/norm(nf1);


v01 = m01*n01;
v02 = m02*n01;

b = -(nf1'*v01 + nf1'*v02);
c = v01'*v02;

x1 = (-b-sqrt(b^2 - 4*c))/2;
x2 = (-b+sqrt(b^2 - 4*c))/2;

mf1 = x1;

vf1 = mf1*nf1;
vf2 = v01 + v02 - vf1;

fprintf('%f %f\n\n',x1,x2);
fprintf('%f %f %f\n',v01 + v02);
fprintf('%f %f %f\n\n',vf1 + vf2);
fprintf('%f\n',v01'*v01 + v02'*v02);
fprintf('%f\n\n\n\n',vf1'*vf1 + vf2'*vf2);

mf1 = x2;

vf1 = mf1*nf1;
vf2 = v01 + v02 - vf1;