function [dudot,dp,dyddot,err] = ...
    main_fem_step(gradG,resG,udot_am,u_af,p,yddot_am,ydot_af,y_af,solver)
% Input : mesh_file = file name where mesh data is stored
%         Ff = force vector function/expression
%         BC = cell array containing BC information
% Output: u = solution vector
%         f = force vector
%         x = x coords
%         y = y coords
%         conn = connectivity array

% Zero the tangent matrix and residual
gradG(:,:) = 0.0; resG(:) = 0.0;

% Unpack variables from structures
conn = solver.mesh.conn;
X = solver.mesh.x;
Y = solver.mesh.y;
nelem = solver.mesh.nelem;
ng = solver.gauss.ng;
ksi_G = solver.gauss.ksi_G;
eta_G = solver.gauss.eta_G;

% Set up function handle for shape functions for parent element
% Could be moved out of the time loop !!
shpg = @(n,J_mat) tri_shp_grad(n,J_mat);
shp  = @(n,ksi,eta) tri_shp(n,ksi,eta);

% Get shape functions ahead of time
shp_mat  = zeros(3);
for nnb=1:3
    for ig=1:ng
        shp_mat(nnb,ig) = shp(nnb,ksi_G(ig),eta_G(ig));
    end
end

% Loop elements
for ne=1:nelem
    % Order points to correspond with convention, get node locations
    % Could this be taken out of the time loop?
    [conn(ne,:),xl,yl] = tri_reorder(X,Y,conn(ne,:));
    
    % Get local arrays of variables in residual
    udot_aml = [udot_am(2*conn(ne,:)-1), udot_am(2*conn(ne,:))];
    u_afl = [u_af(2*conn(ne,:)-1), u_af(2*conn(ne,:))];
    pl = p(conn(ne,:));
    yddot_aml = [yddot_am(2*conn(ne,:)-1), yddot_am(2*conn(ne,:))];
    ydot_afl  = [ydot_af(2*conn(ne,:)-1),  ydot_af(2*conn(ne,:)) ];
    y_afl     = [y_af(2*conn(ne,:)-1),     y_af(2*conn(ne,:))    ];

    % Calculate Jacobian and shape gradients based on reference config
    [Jr,J_matr,~] = tri_Jac(xl,yl);
    % Get shape gradients (depend on Jacobian matrix)
    shpg_matr = zeros(2,3);
    for nnb=1:3
        shpg_matr(:,nnb) = shpg(nnb,J_matr);
    end
    
    % Depending on fluid or solid element:
    if (solver.mesh.tag(ne)) % solid
        % Get local tangent matrix and residual, then assemble
        [gradGl,resGl] = tri_localtanres_LE(solver.gauss,shp_mat,...
            shpg_matr,yddot_aml,y_afl,solver.coeff);
        [gradG,resG] = global_assemble_s(gradG,resG,gradGl,resGl,...
            ne,conn(ne,:),Jr,solver.BC);
    else % (solver.mesh.tag(ne)==0) fluid
        % Calculate Jacobian for fluid, where mesh can change
        [Jf,J_matf,eG] = tri_Jac(xl+y_afl(:,1),yl+y_afl(:,2));
        % Get shape gradients (depend on Jacobian matrix)
        shpg_matf = zeros(2,3);
        for nnb=1:3
            shpg_matf(:,nnb) = shpg(nnb,J_matf);
        end
        % Get local tangent matrix and residual from fluid eq and mesh eq
        [gradGlf,resGlf] = tri_localtanres_NS(solver.gauss,shp_mat,...
            shpg_matf,u_afl,udot_aml,pl,ydot_afl,solver.coeff,eG);
        [gradGlm,resGlm] = tri_localtanres_LE_mesh(shpg_matr,Jf,...
            y_afl,solver.coeff);
        % Assemble
        [gradG,resG] = global_assemble_f(gradG,resG,gradGlf,gradGlm,...
            resGlf,resGlm,ne,conn(ne,:),Jr,Jf,solver.BC);
    end
end

% Mandate that all pressure adds up to 0
for nn = 1:solver.mesh.nnode
    gradG(5*solver.mesh.nnode+1,5*nn-2) = 1.0;
end

% Implement dirichlet conditions
[gradG,resG] = global_BCs(gradG,resG,solver.BC);

% Solve matrix equation
sgradG = sparse(gradG);
sresG  = sparse(resG);
du = sgradG\(-sresG);

% Compute residual error
err = norm(resG);

% Isolate dudot, dp, and dyddot
dudot = zeros(2*solver.mesh.nnode,1);
dp = zeros(solver.mesh.nnode,1);
dyddot = zeros(2*solver.mesh.nnode,1);
for nn = 1:solver.mesh.nnode
    dudot(2*nn-1:2*nn)  = du(5*nn-4:5*nn-3);
    dp(nn)              = du(5*nn-2);
    dyddot(2*nn-1:2*nn) = du(5*nn-1:5*nn);
end