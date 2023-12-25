global tau
global N
global M
global h
global gamma
global d
global valpha
tau = 0.01665;
N = 100;
M = 100;
h = [0.025 0.025];
gamma = 2.04;
d = 0.4;
valpha = 0.1;
global a
global b
global Q1
global Q2
global alpha
global beta
global g
[a,b,Q1,Q2,alpha,beta,g] = input_values(gamma, h, tau, N, M, d, valpha);

global sigma11
global sigma12
global sigma22
global F1
global F2
global F3
global F4
global F5
global F6
sigma11 = zeros(N, M);
sigma12 = zeros(N, M);
sigma22 = zeros(N, M);
F1 = zeros(N, M);
F2 = zeros(N, M);
F3 = zeros(N, M);
F4 = zeros(N, M);
F5 = zeros(N, M);
F6 = zeros(N, M);
%initial condition
w_n = zeros(6, N, M);
%boundary conditions
%D={(x,y): 0<=x<=1, 0<=y<=1}
w_n_next = next_layer_calc(w_n, Q1, Q2);
disp(w_n_next(:,99,99));

