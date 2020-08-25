function [dual]=func_Glopti(nx,g,d,yL)

%% measure-moment formulation:
% max_{mu} vol(mu) 
%   s.t mu supported on g(x)>=0        : support constraint
%       mu <= given Lebesgue Measure   : measure constraint
% where mu is unknown measure defined on safety set g(x)>=0.
% Risk = Expected value of {Polynomial Indicator function of the set {x:g(x)>=0} } with respect to the probability distribution of random vector x.
% Solution of the dual SDP is the coefficients of the polynomial indicator function.
% Lecture 10: Probabilistic Nonlinear Safety Verification, rarnop.mit.edu 
%%
clc

% SDP solvers
mset clear; warning('off','YALMIP:strict')
mset('yalmip',true);mset(sdpsettings('solver','mosek')); % SDP sovers: mosek, sedumi, sdpt3,...

% Lebesgue Measure over [-1,1]^2
mpol('x',1,nx); mu = meas(x); 
y = mom(mmon(x,2*d)); % unknown moments of mu 

% slack measure
mpol('xs',1,nx); mu_s = meas(xs);
ys = mom(mmon(xs,2*d)); % unknown moments

% moment SDP: msdp( max(vol(mu)), measure constraint, support constraint);
% measure constraint: mu <= given Lebesgue Measure ----> mu + mu_s = given Lebesgue Measure  
% -----> measure constraint in terms of the moments:  y+ys=yL
P = msdp(max(mass(mu)), y+ys==yL, g>=0); 

% solve moment SDP relaxation
[stat,obj,mm,dual]=msol(P); 


end