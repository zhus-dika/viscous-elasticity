function w_n_next = next_layer_calc(w_n, Q1, Q2)
  %w_n - 6 x N x M - dimensional tensor
  global tau
  global N
  global M
  global h
  global gamma
  global d
  global valpha
  global a
  global b
  global Q1
  global Q2
  global alpha
  global beta
  global g
  global sigma11
  global sigma12
  global sigma22
  global F1
  global F2
  global F3
  global F4
  global F5
  global F6
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
  r_1_n_half(3,:,:) = r_1_n(3,:,:);
  r_1_n_half(4,:,:) = r_1_n(4,:,:);
  for l=1:N-1
    for k=1:M
      r_1_n_half(5,l,k) = (1-beta(1)) * r_1_n(5,l,k) + beta(1) * r_1_n(5,l+1,k);
      r_1_n_half(6,l,k) = (1-alpha(1)) * r_1_n(6,l,k) + alpha(1) * r_1_n(6,l+1,k);
    end
  end
  %boundary conditions
  for k=1:M
    r_1_n_half(1,1,k) = - r_1_n_half(5,1,k) + gamma * sqrt(2) * sigma12(1,k);
    r_1_n_half(2,1,k) = - r_1_n_half(6,1,k) + sqrt(2) * sigma11(1,k);
    r_1_n_half(5,N,k) = - r_1_n_half(1,N,k) + gamma * sqrt(2) * sigma12(N,k);
    r_1_n_half(6,N,k) = - r_1_n_half(2,N,k) + sqrt(2) * sigma11(N,k);
  end
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
      r_2_n_next(2,l,k) = (1-alpha(2)-tau*g(2)) * r_2_n_half(2,l,k) + alpha(2) * r_2_n_half(2,l,k-1)-tau*g(3)*r_2_n_half(3,l,k)+tau*g(2)*r_2_n_half(6,l,k)-tau*F2(l,k);
    end
  end
  for l=1:N
    for k=1:M
      r_2_n_next(3,l,k) = (1-tau*g(4))*r_2_n_half(3,l,k)-tau*g(3)*r_2_n_half(2,l,k)+tau*g(3)*r_2_n_half(6,l,k)-tau*F3(l,k);
      r_2_n_next(4,l,k) = (1-2*tau*g(1))*r_2_n_half(4,l,k)-tau*F4(l,k);
    end
  end
  for l=1:N
    for k=1:M-1
      r_2_n_next(5,l,k) = (1-beta(2)-tau*g(1)) * r_2_n_half(5,l,k) + beta(2) * r_2_n_half(5,l,k+1) - tau * g(1) * r_2_n_half(1,l,k) - tau * F5(l,k);
      r_2_n_next(6,l,k) = (1-alpha(2)-tau*g(2)) * r_2_n_half(6,l,k) + alpha(2) * r_2_n_half(6,l,k+1)+tau*g(2)*r_2_n_half(2,l,k)+tau*g(3)*r_2_n_half(3,l,k)-tau*F6(l,k);
    end
  end
  %boundary conditions
  for l=1:N
    r_2_n_next(1,l,1) = - r_2_n_next(5,l,1) + gamma * sqrt(2) * sigma12(l,1);
    r_2_n_next(2,l,1) = - r_2_n_next(6,l,1) - sqrt(2) * sigma22(l,1);
    r_2_n_next(5,l,M) = - r_2_n_next(1,l,M) + gamma * sqrt(2) * sigma12(l,M);
    r_2_n_next(6,l,M) = - r_2_n_next(2,l,M) + sqrt(2) * sigma22(l,M);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %w ^ (n + 1) layer calculation
  w_n_next = w_n;
  for l=1:N
    for k=1:M
      w_n_next(:,l,k) = transpose(Q2) * r_2_n_next(:,l,k);
    end
  end

 end
