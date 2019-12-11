function fl = tri_localf(ng,ksi_G,eta_G,W_G,shp,Ff)
% Compute local force vector for linear triangle element
% Inputs: ng = number of Gauss points
%         ksi_G = ksi coords of Gauss points
%         eta_G = eta coords of Gauss points
%         W_G = weights of Gauss points
%         shp = functional form of shape function
%         Ff = functional form of forcing function
% Output: fl = local stiffness matrix

fl= zeros(2*3,1);
% Loop shape functions
for nna=1:3
    % Loop Gauss points
    for ig=1:ng
        fl(2*nna-1)=fl(2*nna-1)+shp(nna,ksi_G(ig),eta_G(ig))*Ff{1}*W_G(ig);
        fl(2*nna  )=fl(2*nna  )+shp(nna,ksi_G(ig),eta_G(ig))*Ff{2}*W_G(ig);
    end
end

end