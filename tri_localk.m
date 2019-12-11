function kl = tri_localk(ng,ksi_G,eta_G,W_G,shpa,shpb)
% Compute local stiffness matrix for linear triangle element
% Inputs: ng = number of Gauss points
%         ksi_G = ksi coords of Gauss points
%         eta_G = eta coords of Gauss points
%         W_G = weights of Gauss points
%         shpa = functional form of shape function with index a
%         shpb = functional form of shape function with index b
% Output: kl = local stiffness matrix

kl= zeros(3);
% Loop shape functions
for nna=1:3
    for nnb=1:3
        % Loop Gauss points
        for ig=1:ng
            % Sum over Gauss pts of Na*Nb*wt
            kl(nna,nnb) = kl(nna,nnb) ...
                        + shpa(nna,ksi_G(ig),eta_G(ig)) ...
                        * shpb(nnb,ksi_G(ig),eta_G(ig))*W_G(ig);
        end
    end
end  