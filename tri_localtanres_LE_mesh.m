function [gradG,resG] = tri_localtanres_LE_mesh(shpg_mat,Jc,yl,c)
% Compute local tangent matrix and residual vector for
% linear elastostatic in RBVMS framework to model mesh motion

Jc = Jc^(-c.ma);

% Construct D, which is 3x3 in 2D
D = zeros(3);
D(1,1) = c.lambda_m*Jc+2*c.mu_m*Jc; D(2,2) = D(1,1);
D(1,2) = c.lambda_m*Jc            ; D(2,1) = D(1,2);
D(3,3) = c.mu_m*Jc;

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
