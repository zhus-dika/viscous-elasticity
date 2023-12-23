function [alpha, beta, gamma, h, tau, N, M, a, b, Q1, Q2] = input_values()
  a = [1 1];
  b = [1 1];
  Q1 = 1 / sqrt(2) * [0 -1 0 0 0 1; -1 0 0 0 1 0; 0 a(2) -2 * b(2) 0 0 a(2); 0 b(2) a(2) 1 0 b(2); 0 b(2) a(2) -1 0 b(2); 1 0 0 0 1 0];
  Q2 = 1 / sqrt(2) * [-1 0 0 0 1 0; 0 1 0 0 0 1; 0 -b(2) a(2) 1 0 b(2); 0 -a(2) -2 * b(2) 0 0 a(2); 0 -b(2) a(2) -1 0 b(2); 1 0 0 0 1 0];
  tau = 0.1;
  gamma = 1.1;
  N = 100;
  M = 100;
  h = [0.01 0.01];
  alpha = [tau/h(1) tau/h(2)];
  beta = alpha / gamma;
 end
