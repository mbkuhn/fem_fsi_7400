function [gradG,resG] = tri_localtanres_LE(gauss,shp_mat,shpg_mat,...
    yddotl,yl,c)
% Compute local stiffness matrix for linear triangle element
% les = Linear ElastoStatic
% Inputs: J_mat= Jacobian matrix
%         shpa = functional form of shape function with index a
%         shpb = functional form of shape function with index b
% Output: kl = local stiffness matrix

ng = gauss.ng;
W_G = gauss.W_G;

% Construct D, which is 3x3 in 2D
D = zeros(3);
D(1,1) = c.lambda_s+2*c.mu_s; D(2,2) = D(1,1);
D(1,2) = c.lambda_s         ; D(2,1) = D(1,2);
D(3,3) = c.mu_s;

BaT = zeros(6,3);
Bb = zeros(3,6);

yl_vec = [yl(1,:) yl(2,:) yl(3,:)]';

% Construct Ba^T
BaT(1:2:5,1)=shpg_mat(1,:)';
BaT(2:2:6,2)=shpg_mat(2,:)';
BaT(1:2:5,3) = BaT(2:2:6,2); 
BaT(2:2:6,3) = BaT(1:2:5,1);
% Calculate B^T_A D
tmp0 = BaT*D;
Bb(1,1:2:5)=shpg_mat(1,:);
Bb(2,2:2:6)=shpg_mat(2,:);
Bb(3,1:2:5) = Bb(2,2:2:6);
Bb(3,2:2:6) = Bb(1,1:2:5);
tmp1 = tmp0*Bb;
% Calculate terms
gradG= tmp1*c.afbdt2;
resG = tmp1*yl_vec;

% Loop shape functions
for nna=1:3
    for nnb=1:3
        rinds = 2*nna-1:2*nna; cinds = 2*nnb-1:2*nnb; 
        
        % Sum over Gauss pts for acceleration term (and force)
        for ig=1:ng
            tmp2=shp_mat(nna,ig)*shp_mat(nnb,ig)*c.rho_s*W_G(ig);
            % Terms are only on the diagonal
            gradG(rinds(1),cinds(1)) = gradG(rinds(1),cinds(1))+tmp2*c.am;
            gradG(rinds(2),cinds(2)) = gradG(rinds(2),cinds(2))+tmp2*c.am;
            resG(rinds) = resG(rinds) + tmp2*(yddotl(nnb,:)'-c.f);
        end
    end
end  