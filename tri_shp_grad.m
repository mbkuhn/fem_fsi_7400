function N = tri_shp_grad(n,J_mat)
% Shape function gradient for linear triangular element
N = zeros(2,1);
N(1) = max(2-n,0)*1 ...      % N_1 shape fcn grad wrt ksi
      +(1-abs(2-n))*0 ...    % N_2 shape fcn grad wrt ksi
      +max(n-2,0)*-1;        % N_3 shape fcn grad wrt ksi
N(2) = max(2-n,0)*0 ...      % N_1 shape fcn grad wrt eta
      +(1-abs(2-n))*1 ...    % N_2 shape fcn grad wrt eta
      +max(n-2,0)*-1;        % N_3 shape fcn grad wrt eta
N = J_mat \ N;  % Multiply by inverse of the Jacobian matrix
end