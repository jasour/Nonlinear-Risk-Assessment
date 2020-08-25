function [A,C,b,Mind]=func_Standard_SDP_Gen(nx,d,g,yL)

%% polynomial approximation of Indicator function of safety set g(x)>=0
% Find polynomial Indicator function P(x):  
% min_{p(x)} Integral_[-1,1]^nx P(x) dx    :obj
% s.t. P(x)>=1  on g(x)>=0          :cons(1)
%      P(x)>=0                      :cons(2)

% SOS Optimization:
% Find polynomial Indicator function P(x):  
% min_{P(x),sigma_0(x),sigma_1(x)} Integral_[-1,1]^nx P(x) dx    :obj
% s.t. P(x)-1= sigma_0(x)+sigma_1(x)g(x)        :cons(1)  ----> Coeffs of right and left habd sides match
%      P(x)>=0 , sigma_0(x)>=0 , sigma_1(x)>=0  :cons(2)  ----> Gram matrix of P(x) and sigma_0 , sigma_1 is PSD

% SDP Optimization:
% Find polynomial Indicator function P(x):  
% min_{Q,Q0,Q1} Integral_[-1,1]^nx P(x) dx    :obj
% s.t. <Ap_i,Q>-<As0_0,Q0>-<As1g_i,Q1> =0 for i>0  : cons(1)
%      <Ap_i,Q>-<As0_0,Q0>-<As1g_i,Q1> =1 for i=0       
%      Q >=0 , Q0 >=0, Q1>=0                       :cons(2)
% where 
%Q: Gram matrix of P(x), Q0:Gram matrix of sigma_0(x), and Q1: Gram matrix of sigma_1(x)

%Lecture 5: Duality of SOS and Moment based Semidefinite Programs(SDPs), Appendix II
% https://rarnop.mit.edu/Lectures-Codes
%%
clc;
nvar=nx; 
Ny= round(factorial(nvar+2*d)/(factorial(2*d)*factorial(nvar))); % Number of coefficients of P(x): Poly of order 2*d in x  

%% 1: P(x)=Sum p_i*x^i = x'Qx,  Q: Gram matrix
% coeffs of P(x): <Ap_i,Q> where Q is the Gram matrix of P(x) polynomial
% integral of P: int P(x)dx = Sum p_i*E[x^i] with respect to Lebesgue Measure
%= Sum p_i*moments of order "i" of Lebesgue Measure 
%=   Mind*Q, where Mind is the moment matirx of order 2d of Lebesgue Measure over [-1,1]^2 
[Ap,Mind]=func_standard_moment(nvar,d);

%% 2: sigma_1(x)*g(x) 
% coeffs of sigma_1(x)*g(x): <As1g_i,Q1> where Q1 is the Gram matrix of sigma_1(x) polynomial
[As1g]=func_standard_localizing(g,nvar,d);

%% Generate A, C, b in Standard SDP min <C,X> s.t. <A,X>=b
% where X=blkdiag{Q,Q0,Q1} 

% <Ap_i,Q>-<As0_0,Q0>-<As1g_i,Q1>
As0=cellfun(@(x) x*(-1),Ap,'un',0); 
AAs1g=cellfun(@(x) x*(-1),As1g,'un',0);
for i=1:Ny; A{i}=   blkdiag(Ap{i},As0{i},AAs1g{i}); end

%  Mind*Q
C=blkdiag(yL(Mind),zeros(size(As0{1})),zeros(size(AAs1g{1})));

m=size(Ap,2); 
b=zeros(m,1); b(1)=1;


end

 
