function [gradG,resG] = tri_localtanres_NS(gauss,shp_mat,shpg_vec,...
    ul,udl,pl,ydl,c,eG)
% Compute local tangent matrix and residual vector for
% an element according to ALE NS equations

ng = gauss.ng;
W_G = gauss.W_G;

% 3 shape functions, 2 directions of udot, 1 direction of p (scalar)
gradG= zeros(3*3);
resG= zeros(3*3,1);

% Calculate things that do not depend on gauss points
shpgshpg = shpg_vec'*shpg_vec;
% gug= grad u_gauss, gpg = grad p_gauss
gpg = (shpg_vec*pl)';
gug = shpg_vec*ul;

% Loop gauss points
for ig=1:ng
    % Initialize quantities at gauss points
    pg = 0.0;
    ug = zeros(1,2); 
    udg = ug; umg = ug; % udot_g, umesh_g
    % Calculate values at gauss points
    for nnb=1:3
        % Store shape function
        shpB = shp_mat(nnb,ig);
        % Calculate u, udot, p, um at gauss points
        ug  = ug +ul(nnb,:)*shpB;
        udg = udg+udl(nnb,:)*shpB;
        pg  = pg +pl(nnb)*shpB;
        umg = umg+(ul(nnb,:)-ydl(nnb,:))*shpB;
    end
    % Calculate parameters
    tau_m = (umg*(eG*umg')+...
        c.c1*(c.mu_f/c.rho_f)^2*sum(sum(eG.*eG))+4/c.dt^2)^(-0.5);
    nu_c = 1/(trace(eG)*tau_m);
    % Calculate pointwise residuals
    Rm = c.rho_f*udg + c.rho_f*umg*gug + gpg;
    Rc = sum(diag(gug));
    % Calculate primed quantities
    upg = -tau_m/c.rho_f*Rm;
    ppg = -c.rho_f*nu_c*Rc;
    
    % Loop shape functions
    for nna=1:3
        % Store shape function and gradient
        shpgA = shpg_vec(:,nna);
        shpA = shp_mat(nna,ig);
        
        % Calculate repetitively used terms
        umgshpgA = umg*shpgA;
        
        % Calculate residual vector components
        % Loop vector indices
        for ii=1:2
            % Row indices
            rind = 3*nna-2+ii-1;
            
            % Term 1 in Rm
            resG(rind) = resG(rind) ...
                + shpA*udg(ii)*W_G(ig)*c.rho_f;
            % Term 2 in Rm
            resG(rind) = resG(rind) ...
                + shpA*umg*gug(:,ii)*W_G(ig)*c.rho_f;
            % Term 3 in Rm (ignore)
            % Term 4 in Rm
            resG(rind) = resG(rind) ...
                + (shpgA(1)*(gug(1,ii)+gug(ii,1)) ...
                +  shpgA(2)*(gug(2,ii)+gug(ii,2)))...
                * W_G(ig)*c.mu_f/2;
            % Term 5 in Rm
            resG(rind) = resG(rind) ...
                + shpgA(ii)*pg*W_G(ig)*(-1);
            % Term 6 in Rm
            resG(rind) = resG(rind) ...
                +shpA*upg*gug(:,ii)*W_G(ig)*c.rho_f;
            % Term 7 in Rm
            resG(rind) = resG(rind) ...
                +umgshpgA*upg(ii)*W_G(ig)*c.rho_f*(-1);
            % Term 8 in Rm
            resG(rind) = resG(rind) ...
                +umgshpgA*upg(ii)*W_G(ig)*c.rho_f*(-1);
            % Term 9 in Rm
            resG(rind) = resG(rind) ...
                +shpgA(ii)*ppg*W_G(ig)*c.rho_f*(-1);
            % Term 10 in Rc
            resG(3*nna) = resG(3*nna) ...
                + shpA*gug(ii,ii)*W_G(ig);
            % Term 11 in Rc
            resG(3*nna) = resG(3*nna) ...
                + shpgA(ii)*upg(ii)*W_G(ig)*(-1);
        end
        
        % Loop shape functions again for tangent matrix
        for nnb=1:3
            % Store shape function and gradient
            shpgB = shpg_vec(:,nnb);
            shpB = shp_mat(nnb,ig);
            % Calculate repetitively used terms
            shpgA_shpgB = shpgshpg(nna,nnb);
            umgshpgB   = umg*shpgB;
            % Row and column indices (ii==jj) for these terms
            % Term 1 in dRm/dudot
            term1 = shpA*shpB*c.rho_f*c.am;
            % Term 2 in dRm/dudot
            term2 = umgshpgB*c.rho_f*c.afgdt;
            % Term 3 in dRm/dudot, part 1
            term3 = shpgA_shpgB*c.mu_f*c.afgdt/2;
            % Term 4 in dRm/dudot
            term4 = umgshpgA*shpB*tau_m*c.rho_f*c.am;
            % Term 5 in dRm/dudot
            term5 = umgshpgA*umgshpgB*tau_m*c.rho_f*c.afgdt;
            % Combine
            term12345 = (term1+term2+term3+term4+term5)*W_G(ig);
            % Loop vector indices
            for ii=1:2
                % Incorporate diagonal terms
                rind = 3*nna-2+ii-1; cind = 3*nnb-2+ii-1;
                % Term 1 in dRm/dudot
                gradG(rind,cind) = gradG(rind,cind) ...
                    + term12345;
                for jj=1:2
                    % Row and column indices
                    rind = 3*nna-2+ii-1; cind = 3*nnb-2+jj-1;
                    
                    % Term 3 in dRm/dudot, part 2
                    gradG(rind,cind) = gradG(rind,cind) ...
                        + shpgA(jj)*shpgB(ii)...
                        * W_G(ig)*c.mu_f*c.afgdt/2;

                    % Term 6 in dRm/dudot
                    gradG(rind,cind) = gradG(rind,cind) ...
                        + shpgA(ii)*shpgB(jj)...
                        * W_G(ig)*nu_c*c.rho_f*c.afgdt;
                end
                % Row and column indices
                rind = 3*nna-2+ii-1; cind = 3*nnb;
                
                % Term 7 in dRm/dp
                gradG(rind,cind) = gradG(rind,cind) ...
                    + shpgA(ii)*shpB ...
                    * W_G(ig)*(-1);
                % Term 8 in dRm/dp
                gradG(rind,cind) = gradG(rind,cind) ...
                    + umgshpgA*shpgB(ii)...
                    * W_G(ig)*tau_m*(-1);
            end
            
            for jj=1:2
                % Row and column indices
                rind = 3*nna; cind = 3*nnb-2+jj-1;
                
                % Term 9 in dRc/dudot
                gradG(rind,cind) = gradG(rind,cind) ...
                    + shpA*shpgB(jj)...
                    * W_G(ig)*c.afgdt;
                % Term 10 in dRc/dudot
                gradG(rind,cind) = gradG(rind,cind) ...
                    + shpgA(jj)*shpB...
                    * W_G(ig)*tau_m*c.am;
                % Term 11 in dRc/dudots
                gradG(rind,cind) = gradG(rind,cind) ...
                    + shpgA(jj)*umgshpgB...
                    * W_G(ig)*tau_m*c.afgdt;
            end
            
            % Term 12 in dRc/dp
            gradG(3*nna,3*nnb) = gradG(3*nna,3*nnb) ...
                + shpgA_shpgB...
                * W_G(ig)*tau_m/c.rho_f;
        end
    end
end
