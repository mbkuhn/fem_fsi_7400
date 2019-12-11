function [J,J_mat,G] = tri_Jac(x,y)
% Computes the Jacobian for a linear triangular element
% inputs must have 3 points, they are assumed to be in the convention
J_mat = zeros(2);
J_mat(1,1) = x(1)-x(3); % dx/dksi
J_mat(1,2) = y(1)-y(3); % dy/dksi
J_mat(2,1) = x(2)-x(3); % dx/deta
J_mat(2,2) = y(2)-y(3); % dy/deta

J = det(J_mat);

% element metric
G_mat = inv(J_mat);
G = zeros(2);
for i = 1:2
    for j = 1:2
        for k = 1:2
            G(i,j) = G(i,j) + G_mat(i,k)*G_mat(j,k);
        end
    end
end

end