# Nonlinear-Risk-Assessment

# Risk = probability (g(x)>=0)  where x is a random vector and g(x)>=0 is a nonlinear safety constraint

Eaxample:  Probability ( -x(1)^4+0.5*(x(1)^2-x(2)^2)+0.1 >=0 ) 
where, x1 hase Uniform probability distribution over [0,1] and x2 has Beta(1.3,3) probability distribution over [0,1].

More information: rarnop.mit.edu, Lecture 10: Probabilistic Nonlinear Safety Verification.

1) Nonlinear-Risk-Assessment in CVX (Standard SDP Formulation)

2) Nonlinear-Risk-Assessment in GloptiPoly (moment SDP Formulation)

3) Nonlinear-Risk-Assessment in Yalmip (sum-of-squares (SOS) SDP Formulation)
