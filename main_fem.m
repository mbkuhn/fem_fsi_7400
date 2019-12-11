function [u,f,x,y,conn] = main_fem(mesh_file,Ff,BC,sprop,ndim_out)
% Input : mesh_file = file name where mesh data is stored
%         Ff = force vector function/expression
%         BC = cell array containing BC information
% Output: u = solution vector
%         f = force vector
%         x = x coords
%         y = y coords
%         conn = connectivity array

% Get constants from solid properties structure
lambda = sprop.lambda;
mu = sprop.mu;

% Read mesh files
[x,y,conn,nx,ny,nelem,npe] = read_mesh(mesh_file);
Lx = max(x)-min(x); Ly = max(y)-min(y);

% Initialize global arrays (2D)
k = zeros(ndim_out*nx*ny); f = zeros(ndim_out*nx*ny,1); 

% Set up Gauss pts (linear triangle, ng = 3)
ng = 3;
[ksi_G,eta_G,W_G] = tri_Gpts_3(); 

% Loop elements
for ne=1:nelem
    % Order points correspond with convention, get node locations
    [conn(ne,:),xl,yl] = tri_reorder(x,y,conn(ne,:));
        
    % Calculate Jacobian for element
    [J,J_mat] = tri_Jac(xl,yl);
    % J_mat = [dx/dksi, dx/deta]
    %         [dy/dksi, dy/deta]
    
    % Set up function handle for shape functions for parent element
    shpg = @(n,J_mat) tri_shp_grad(n,J_mat);
    shp  = @(n,ksi,eta) tri_shp(n,ksi,eta);
    
    % Compute local stiffness matrix
    kl = tri_localk_les(J_mat,shpg,shpg,lambda,mu);
    
    % Compute local force vector
    fl = tri_localf(ng,ksi_G,eta_G,W_G,shp,Ff);
    
    % Assemble local stiffness matrix and force vector into global
    [k,f] = global_assemble(k,f,kl,fl,x,y,conn(ne,:),J,BC);
end

% Solve matrix equation
u = k\f;