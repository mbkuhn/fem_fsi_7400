function resG = tri_localres_heatdiff(gauss,J_mat,shp,shpg,coeff,udl,ul)
% Compute local stiffness matrix for linear triangle element
% Inputs: ng = number of Gauss points
%         ksi_G = ksi coords of Gauss points
%         eta_G = eta coords of Gauss points
%         W_G = weights of Gauss points
%         shpa = functional form of shape function with index a
%         shpb = functional form of shape function with index b
% Output: kl = local stiffness matrix

ng = gauss.ng;
ksi_G = gauss.ksi_G;
eta_G = gauss.eta_G;
W_G = gauss.W_G;

resG= zeros(3,1);
% Loop shape functions
for nna=1:3
    for nnb=1:3
        % Loop Gauss points
        for ig=1:ng
            % Sum over Gauss pts of Na*Nb*wt
            resG(nna) = resG(nna) ...
                        + shp(nna,ksi_G(ig),eta_G(ig)) ...
                        * shp(nnb,ksi_G(ig),eta_G(ig)) ...
                        * W_G(ig)*coeff(1)*udl(nnb) ...
                        + sum(shpg(nna,J_mat).*shpg(nnb,J_mat))...
                        * W_G(ig)*coeff(2)*ul(nnb);
        end
    end
end  