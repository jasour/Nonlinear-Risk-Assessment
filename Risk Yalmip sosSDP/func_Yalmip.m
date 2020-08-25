function [c_sos]=func_Yalmip(nx,x,g,d,yL)
%% sum-of-squares formulation:
% Find polynomial Indicator function P(x):  
% min_{p(x)} Integral_[-1,1]^nx P(x) dx    :obj
% s.t. P(x)>=1  on g(x)>=0          :cons(1)
%      P(x)>=0                      :cons(2)

% SOS Optimization:
% Find polynomial Indicator function P(x):  
% min_{P(x),sigma_0(x),sigma_1(x)} Integral_[-1,1]^nx P(x) dx    :obj
% s.t. P(x)-1= sigma_0(x)+sigma_1(x)g(x)        :cons(1)--->(P(x)-1-sigma_1(x)g(x)) is SOS polynomial
%      P(x)>=0 , sigma_0(x)>=0 , sigma_1(x)>=0  :cons(2)--->P(x),sigma_0(x),sigma_1(x) are SOS polynomials

% Lecture 10: Probabilistic Nonlinear Safety Verification, rarnop.mit.edu 
%%
clc
% SDP Solver
ops = sdpsettings('solver','mosek');

% polynomial Indicator function P(x) 
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(nx,k)]; end % monomials
cp=sdpvar(size(vpow,1),1); %coefficients 
X=1;for i=1:nx; X=X.*x(i).^vpow(:,i); end
P=cp'*X; % polynomial Indicator function P(x) 
P_Int=cp'*yL;% Integral of P(x) with respect to Lebesgue Measure over [-1,1]^2

% sigma_1(x)
[s1,c1] = polynomial(x,2*d);

% SOS Conditions
F = [sos(P-1-[s1]*g), sos(P), sos(s1)];
% Solve SOS based SDP
[sol,v,Q]=solvesos(F,P_Int,ops,[c1;cp]);

% Obtained coefficients of polynomial Indicator function 
c_sos=double(cp);

end