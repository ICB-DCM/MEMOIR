function out = MEMsigma_noise_decay_diag_diag_1(in1)
%MEMSIGMA_NOISE_DECAY_DIAG_DIAG_1
%    out = MEMSIGMA_NOISE_DECAY_DIAG_DIAG_1(IN1)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    05-Dec-2016 13:01:38

p_m0= in1(1);
p_delta= in1(2);
out = zeros(1,1);
out(1) = 1/100;
out = reshape(out,[1  1]);