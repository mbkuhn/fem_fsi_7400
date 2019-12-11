function resG = tri_localres_NS(gauss,J_mat,shp,shpg,coeff,udl,ul,pl)
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

% Get residuals at gauss points for stability terms
Rm = zeros(2,ng);
Rc = zeros(1,ng);
for nnb=1:3
    shpgB = shpg(nnb,J_mat);
    for ig=1:ng
        shpB = shp(nnb,ksi_G(ig),eta_G(ig));
        for ii=1:2
            Rm(ii,ig) = Rm(ii,ig) ...
                + shpB*udl(nnb,ii)*coeff(1) ...
                + shpB*sum(ul(nnb,:)'.*shpgB)*ul(nnb,ii)*coeff(2) ...
                + shpB*coeff(3) ...
                + sum(shpgB.*(shpgB*ul(nnb,ii)))*coeff(4) ...
                + shpgB(ii)*pl(nnb)*coeff(5);
        end
        Rc(ig) = Rc(ig) ...
            + sum(shpgB.*ul(nnb,:)')*coeff(10);
    end
end

% Calculate prime quantities
upl = coeff(12)*Rm;
ppl = coeff(13)*Rc;

% 3 shape functions, 2 directions of udot, 1 direction of p (scalar)
resG= zeros(3*3,1);
% Loop shape functions
for nna=1:3
    for nnb=1:3
        % Store shape function gradients
        shpgA = shpg(nna,J_mat);
        shpgB = shpg(nnb,J_mat);
        % Loop Gauss points
        for ig=1:ng
            % Store shape functions
            shpA = shp(nna,ksi_G(ig),eta_G(ig));
            shpB = shp(nnb,ksi_G(ig),eta_G(ig));
            % Loop vector indices
            for ii=1:2
                % Row indices
                rind = 3*nna-2+ii-1;
                
                % Term 1 in Rm
                resG(rind) = resG(rind) ...
                    + shpA*shpB*udl(nnb,ii)*W_G(ig)*coeff(1);
                % Term 2 in Rm
                resG(rind) = resG(rind) ...
                    + shpA*shpB*sum(ul(nnb,:)'.*shpgB)*ul(nnb,ii)...
                    * W_G(ig)*coeff(2);
                % Term 3 in Rm
                resG(rind) = resG(rind) ...
                    + shpA*W_G(ig)*coeff(3);
                % Term 4 in Rm
                resG(rind) = resG(rind) ...
                    + sum(shpgA.*(shpgB*ul(nnb,ii)+shpgB(ii)*ul(nnb,:)'))...
                    * W_G(ig)*coeff(4);
                % Term 5 in Rm
                resG(rind) = resG(rind) ...
                    + shpgA(ii)*pl(nnb)*W_G(ig)*coeff(5);
                % Term 6 in Rm
                resG(rind) = resG(rind) ...
                    + shpA*sum(upl(:,ig).*shpgB)*ul(nnb,ii)*W_G(ig)*coeff(6);
                % Term 7 in Rm
                resG(rind) = resG(rind) ...
                    + sum(shpgA*upl(ii,ig).*upl(:,ig))*W_G(ig)*coeff(7);
                % Term 8 in Rm
                resG(rind) = resG(rind) ...
                    + sum(shpgA*shpB.*ul(nnb,:)')*upl(ii,ig)*W_G(ig)*coeff(8);
                % Term 9 in Rm
                resG(rind) = resG(rind) ...
                    + shpgA(ii)*ppl(ig)*W_G(ig)*coeff(9);
                
            end
            % Term 10 in Rc
            resG(3*nna) = resG(3*nna) ...
                + shpA*sum(shpgB.*ul(nnb,:)')*W_G(ig)*coeff(10);
            % Term 11 in Rc
            resG(3*nna) = resG(3*nna) ...
                + shpA*sum(shpgB.*upl(:,ig))*W_G(ig)*coeff(11);
        end
    end
end  