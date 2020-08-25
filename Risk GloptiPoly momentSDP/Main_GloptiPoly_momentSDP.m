%% Risk = probability (g(x)>=0)  where x: random vector
% Eaxample:  Probability ( -x(1)^4+0.5*(x(1)^2-x(2)^2)+0.1 >=0 ) 
% where x1 hase Uniform probability distribution over [0,1]
%       x2 has Beta(1.3,3) probability distribution over [0,1]
% Lecture 10: Probabilistic Nonlinear Safety Verification, rarnop.mit.edu 
%% Ashkan Jasour, Research Scientist, MIT 2020
% jasour.mit.edu  rarnop.mit.edu
%%
clc;clear all;close all
%%
nx=2; % number of uncertain parameters
d=10; % relaxation order: SDP uses 2*d number of the moments of uncertainties to calculate the Risk.

%% Safety Constraint : Risk= Prob(g(x)>=0)  x:random vector
mpol('x',1,nx); 
g=-x(1)^4+0.5*(x(1)^2-x(2)^2)+0.1;

%% Moments Information

%moments of Lebesgue Measure over [-1,1]^2 to calculate the integral
u=1;l=-1; yL=[2];for i=1:2*d ;yL(i+1,1)= ( u^(i+1) - l^(i+1) )/(i+1);end 
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(nx,k)]; end; 
yL=prod(yL(vpow+1),2);

%moments of Uniform probability distribution over [0,1] for uncertain variable x1
u=1;l=0;yx1=[1];for i=1:2*d ;yx1(i+1,1)=(1/(u-l))*((u^(i+1) - l^(i+1))/(i+1));end 

%moments of Beta(aB,bB) probability distribution over [0,1] for uncertain variable x2
aB=1.3;bB=3; yx2=[1];for k=1:2*d; yx2=[yx2;(aB+k-1)/(aB+bB+k-1)*yx2(end) ]; end;

% moments of joint distribution of x1 and x2
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(nx,k)]; end; 
yx1x2=yx1(vpow(:,1)+1).*yx2(vpow(:,2)+1);

%% GloptiPoly
% Gloptipoly solves the moment SDP.
% We use the solution of the dual SDP (i.e., c_mom) to calculate the risk.
% c_mom is the coefficients of the polynomial indicator function.
[c_mom]=func_Glopti(nx,g,d,yL);
Risk_mom_dual=sum(yx1x2.*c_mom);

% polynomial indicator function obtianed from the solution of Dual SDP
z=sym('x',[1 nx]); P=sum(prod(z.^vpow,2).*c_mom); 
[x1,x2] = meshgrid(-0.9:0.01:0.9,-0.9:0.01:0.9);P=eval(P);
close all;surfc(x1,x2,P,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.7);
camlight; lighting gouraud; hold on;grid on;set(gca,'fontsize',25)
title('polynomial indicator function');


%% Monte Carlo
x=random('uniform',l,u,1,10^6); x=[x;random('beta',aB,bB,1,10^6)];    
g_monte=-x(1,:).^4+0.5*(x(1,:).^2-x(2,:).^2)+0.1;
Risk_monte=size(find(g_monte>=0),2)/10^6;

%% Compare
clc; 
display('   GloptiPoly_dual,     monte carlo')
display([Risk_mom_dual, Risk_monte])
