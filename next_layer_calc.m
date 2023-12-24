function w_n_next = next_layer_calc(w_n, Q1, Q2)
  %w_n - 6 x N x M - dimensional tensor
  %>> global b
  %>> b=7
  global alpha, beta, gamma, tau, g, sigma11, sigma12, F1, F2, F3, F4, F5, F6;
  N = size(w_n, 2);
  M = size(w_n, 3);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %(r ^ (1)) ^ n layer calculation
  r_1_n = w_n;
  for l=1:N
    for k=1:M
      r_1_n(:,l,k) = transpose(Q1) * w_n(:,l,k);
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %(r ^ (1)) ^ (n + 1/2) layer calculation
  r_1_n_half = r_1_n;
  for l=2:N
    for k=1:M
      r_1_n_half(1,l,k) = (1-beta(1)) * r_1_n(1,l,k) + beta(1) * r_1_n(1,l-1,k);
      r_1_n_half(2,l,k) = (1-alpha(1)) * r_1_n(2,l,k) + alpha(1) * r_1_n(2,l-1,k);
    end
  end
  r_1_n_half(3,l,k) = r_1_n(3,l,k);
  r_1_n_half(4,l,k) = r_1_n(4,l,k);
  for l=1:N-1
    for k=1:M
      r_1_n_half(5,l,k) = (1-beta(1)) * r_1_n(5,l,k) + beta(1) * r_1_n(5,l+1,k);
      r_1_n_half(6,l,k) = (1-alpha(1)) * r_1_n(6,l,k) + alpha(1) * r_1_n(6,l+1,k);
    end
  end
  %boundary conditions
  r_1_n_half(1,1,:) = - r_1_n_half(5,1,:) + gamma * sqrt(2) * sigma12(1,:);
  r_1_n_half(2,1,:) = - r_1_n_half(6,1,:) + sqrt(2) * sigma11(1,:);
  r_1_n_half(5,N,:) = - r_1_n_half(1,N,:) + gamma * sqrt(2) * sigma12(N,:);
  r_1_n_half(6,N,:) = - r_1_n_half(2,N,:) + sqrt(2) * sigma11(N,:);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %(r ^ (2)) ^ (n + 1/2) layer calculation
  r_2_n_half = w_n;
  for l=1:N
    for k=1:M
      r_2_n_half(:,l,k) = transpose(Q2) * Q1 * r_1_n_half(:,l,k);
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %(r ^ (2)) ^ (n + 1) layer calculation
  r_2_n_next = r_2_n_half;
  for l=1:N
    for k=2:M
      r_2_n_next(1,l,k) = (1-beta(2)-tau*g(1)) * r_2_n_half(1,l,k) + beta(2) * r_2_n_half(1,l,k-1) - tau * g(1) * r_2_n_half(5,l,k) - tau * F1(l,k);
      r_2_n_next(2,l,k) = (1-alpha(2)-tau*g(2)) * r_2_n_half(2,l,k) + alpha(2) * r_2_n_half(2,l,k-1)
      -tau*g(3)*r_2_n_half(3,l,k)+tau*g(2)*r_2_n_half(6,l,k)-tau*F2(l,k);
    end
  end
  r_2_n_next(3,:,:) = (1-tau*g(4))*r_2_n_half(3,:,:)-tau*g(3)*r_2_n_half(2,:,:)+tau*g(3)*r_2_n_half(6,:,:)-tau*F3;
  r_2_n_next(4,:,:) = (1-2*tau*g(1))*r_2_n_half(4,:,:)-tau*F4;
   for l=1:N
    for k=1:M-1
      r_2_n_next(5,l,k) = (1-beta(2)-tau*g(1)) * r_2_n_half(5,l,k) + beta(2) * r_2_n_half(5,l,k+1) - tau * g(1) * r_2_n_half(1,l,k) - tau * F5(l,k);
      r_2_n_next(6,l,k) = (1-alpha(2)-tau*g(2)) * r_2_n_half(6,l,k) + alpha(2) * r_2_n_half(6,l,k+1)
      +tau*g(2)*r_2_n_half(2,l,k)+tau*g(3)*r_2_n_half(3,l,k)-tau*F6(l,k);
    end
  end
  %boundary conditions
  r_2_n_next(1,:,1) = - r_2_n_next(5,:,1) + gamma * sqrt(2) * sigma12(:,1);
  r_2_n_next(2,:,1) = - r_2_n_next(6,:,1) - sqrt(2) * sigma22(:,1);
  r_2_n_next(5,:,M) = - r_2_n_next(1,:,M) + gamma * sqrt(2) * sigma12(:,M);
  r_2_n_next(6,:,M) = - r_2_n_next(2,:,M) + sqrt(2) * sigma22(:,M);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %w ^ (n + 1) layer calculation
  w_n_next = w_n;
  for l=1:N
    for k=1:M
      w_n_next(:,l,k) = transpose(Q2) * r_2_n_next(:,l,k);
    end
  end

 end
