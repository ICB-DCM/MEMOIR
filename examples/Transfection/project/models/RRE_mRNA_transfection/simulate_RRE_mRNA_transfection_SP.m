function [y,sy] = simulate_RRE_mRNA_transfection_SP(tout,phi,kappa)

% Simulate model
[status,t,x,y,sx,sy] = simulate_RRE_mRNA_transfection(tout,phi,kappa);