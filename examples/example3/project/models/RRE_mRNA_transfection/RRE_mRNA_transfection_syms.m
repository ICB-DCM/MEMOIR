function [model] = RRE_mRNA_transfection_syms()

%% CVODES OPTIONS

% absolute tolerance
model.atol = 1e-8;
% relative tolerance
model.rtol = 1e-8;
% maximal number of steps
model.maxsteps = 1e4;

%% STATES

syms MRNA GFP

x = [MRNA,GFP];

%% PARAMETERS

syms t0 kTL_m0 pbeta delta

p = [t0,kTL_m0,pbeta,delta];

%% INPUT 

u = sym.empty(0,0);

%% SYSTEM EQUATIONS

xdot = sym(zeros(size(x)));

% RRE-MEAN
xdot(1) = -delta*pbeta*MRNA;
xdot(2) = kTL_m0*MRNA-pbeta*GFP;

%% INITIAL CONDITIONS

x0 = sym(zeros(2,1));

x0(1) = 1;

%% OBSERVABLES

y = sym(zeros(1,1));

y(1) = GFP;


%% SYSTEM STRUCT

model.sym.x = x;
model.sym.u = u;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;
end