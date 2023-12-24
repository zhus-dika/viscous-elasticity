function [gamma, h, tau, N, M, d, valpha] = input_values()
  gamma = 2.04;
  a = [gamma / 3 * (sqrt(2) + 1 / sqrt(3 * gamma^2 - 4)), 1 / (3 * gamma) * (2 * sqrt(2) + sqrt(3 * gamma^2 - 4))];
  b = [gamma / 3 * (-sqrt(2) / 2 + 1 / sqrt(3 * gamma^2 - 4)), 1 / (3 * gamma) * (-sqrt(2) + sqrt(3 * gamma^2 - 4))];
  Q1 = 1 / sqrt(2) * [0 -1 0 0 0 1; -1 0 0 0 1 0; 0 a(2) -2 * b(2) 0 0 a(2); 0 b(2) a(2) 1 0 b(2); 0 b(2) a(2) -1 0 b(2); 1 0 0 0 1 0];
  Q2 = 1 / sqrt(2) * [-1 0 0 0 1 0; 0 1 0 0 0 1; 0 -b(2) a(2) 1 0 b(2); 0 -a(2) -2 * b(2) 0 0 a(2); 0 -b(2) a(2) -1 0 b(2); 1 0 0 0 1 0];
  tau = 0.01665;
  N = 100;
  M = 100;
  h = [0.025 0.025];
  alpha = [tau/h(1) tau/h(2)];
  beta = alpha / gamma;
  d = 0.4;
  valpha = 0.1;
  theta = d / (valpha * gamma^2);
  g = [valpha/2, 2 * valpha / (3 * gamma^2), valpha / (3 * gamma^2) * sqrt(2 * (3 * gamma ^ 2 - 4)), valpha / (3 * gamma^2) * (3 * gamma ^ 2 - 4)];
 end
