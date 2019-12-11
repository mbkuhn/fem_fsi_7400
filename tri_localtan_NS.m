function gradG = tri_localtan_NS(gauss,J_mat,shp,shpg,ul,coeff)
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

% 3 shape functions, 2 directions of udot, 1 direction of p (scalar)
gradG= zeros(3*3);
% Loop shape functions
for nna=1:3
    for nnb=1:3
        % Store shape function gradients
        shpgA = shpg(nna,J_mat); shpgB = shpg(nnb,J_mat);
        % Loop Gauss points
        for ig=1:ng
            % Store shape functions
            shpA = shp(nna,ksi_G(ig),eta_G(ig));
            shpB = shp(nnb,ksi_G(ig),eta_G(ig));
            % Loop vector indices
            for ii=1:2
                for jj=1:2
                    % Row and column indices
                    rind = 3*nna-2+ii-1; cind = 3*nnb-2+jj-1;
                    
                    % Term 1 in dRm/dudot
                    gradG(rind,cind) = gradG(rind,cind) ...
                        + shpA*(ii==jj)*shpB ...
                        * W_G(ig)*coeff(1);
                    % Term 2 in dRm/dudot
                    gradG(rind,cind) = gradG(rind,cind) ...
                        + sum(shpgB.*ul(nnb,:)')*(ii==jj)...
                        * W_G(ig)*coeff(2);
                    % Term 3 in dRm/dudot
                    gradG(rind,cind) = gradG(rind,cind) ...
                        +(sum(shpgA.*shpgB)*(ii==jj)+shpgA(jj)*shpgB(ii))...
                        * W_G(ig)*coeff(3);
                    % Term 4 in dRm/dudot
                    gradG(rind,cind) = gradG(rind,cind) ...
                        + sum(shpgA.*ul(nna,:)')*(ii==jj)*shpB...
                        * W_G(ig)*coeff(4);
                    % Term 5 in dRm/dudot
                    gradG(rind,cind) = gradG(rind,cind) ...
                        + sum(shpgA.*ul(nna,:)')*(ii==jj)...
                        * sum(shpgB.*ul(nnb,:)')...
                        * W_G(ig)*coeff(5);
                    % Term 6 in dRm/dudot
                    gradG(rind,cind) = gradG(rind,cind) ...
                        + shpgA(ii)*shpgB(jj)...
                        * W_G(ig)*coeff(6); 
                end
                % Row and column indices
                rind = 3*nna-2+ii-1; cind = 3*nnb;
                
                % Term 7 in dRm/dp
                gradG(rind,cind) = gradG(rind,cind) ...
                    + shpgA(ii)*shpB ...
                    * W_G(ig)*coeff(7);
                % Term 8 in dRm/dp
                gradG(rind,cind) = gradG(rind,cind) ...
                    + sum(shpgA.*ul(nna,:)')*shpgB(ii)...
                    * W_G(ig)*coeff(8);
            end

            for jj=1:2
                % Row and column indices
                rind = 3*nna; cind = 3*nnb-2+jj-1;
                
                % Term 9 in dRc/dudot
                gradG(rind,cind) = gradG(rind,cind) ...
                    + shpA*shpgB(jj)...
                    * W_G(ig)*coeff(9);
                % Term 10 in dRc/dudot
                gradG(rind,cind) = gradG(rind,cind) ...
                    + shpgA(jj)*shpB...
                    * W_G(ig)*coeff(10);
                % Term 11 in dRc/dudot
                gradG(rind,cind) = gradG(rind,cind) ...
                    + shpgA(jj)*sum(shpgB.*ul(nnb,:)')...
                    * W_G(ig)*coeff(11);
            end
            
            % Term 12 in dRc/dp
            gradG(3*nna,3*nnb) = gradG(3*nna,3*nnb) ...
                + sum(shpgA.*shpgB)...
                * W_G(ig)*coeff(12);
        end
    end
end  