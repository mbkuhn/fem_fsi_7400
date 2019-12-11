function N = tri_shp(n,ksi,eta)
% Shape function for linear triangular element
N = max(2-n,0)*ksi ...      % N_1 shape fcn
   +(1-abs(2-n))*eta ...    % N_2 shape fcn
   +max(n-2,0)*(1-ksi-eta); % N_3 shape fcn
end