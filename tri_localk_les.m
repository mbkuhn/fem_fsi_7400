function kl = tri_localk_les(J_mat,shpa,shpb,lambda,mu)
% Compute local stiffness matrix for linear triangle element
% les = Linear ElastoStatic
% Inputs: J_mat= Jacobian matrix
%         shpa = functional form of shape function with index a
%         shpb = functional form of shape function with index b
% Output: kl = local stiffness matrix

% Construct D, which is 3x3 in 2D
D = zeros(3);
D(1,1) = lambda+2*mu; D(2,2) = D(1,1);
D(1,2) = lambda     ; D(2,1) = D(1,2);
D(3,3) = mu         ;

% No need to sum over Gauss pts, shape fcn gradients are constant
kl= zeros(2*3);
% Loop shape functions
for nna=1:3
    for nnb=1:3
        % Construct Ba, Bb
        Ba = zeros(3,2); Bb = Ba;
        tmp = shpa(nna,J_mat); Ba(1,1)=tmp(1);
        tmp = shpa(nna,J_mat); Ba(2,2)=tmp(2);
        Ba(3,1) = Ba(2,2); Ba(3,2) = Ba(1,1);
        tmp = shpb(nnb,J_mat); Bb(1,1)=tmp(1);
        tmp = shpb(nnb,J_mat); Bb(2,2)=tmp(2);
        Bb(3,1) = Bb(2,2); Bb(3,2) = Bb(1,1);
        % Calculate B^T_A D B_B
        kl(2*nna-1:2*nna,2*nnb-1:2*nnb) = transpose(Ba)*D*Bb;
    end
end  