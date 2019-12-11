clear; clc; clf;

%% Memory structures (Modules)
mesh = struct();  % mesh-related parameters
fprop = struct(); % fluid properties
sprop = struct(); % solid properties
time  = struct(); % time-related parameters
iter  = struct(); % iteration-related parameters
gauss = struct(); % Gaussian quadrature
plots =  struct(); % plotting-related values
coeff = struct();  % coefficients needed in local grad and res
solver = struct(); % encasing structure

% Mesh option
mesh_opt = 'meshgen';

%% Inputs
% Dimensions and mesh specified within mesh_gen function
if (strcmp(mesh_opt,'uniform'))
    mesh.nx = 17;
    mesh.ny = 17;
    mesh.nnode = mesh.nx*mesh.ny;
    mesh.Lx = 1.0;
    mesh.Ly = 1.0;
    mesh.filename = 'mesh';
elseif (strcmp(mesh_opt,'bar'))
    mesh.nx = 17;
    mesh.ny = 5;
    mesh.nnode = mesh.nx*mesh.ny;
    mesh.Lx = 10.0;
    mesh.Ly = 1.0;
    mesh.filename = 'mesh';
elseif (strcmp(mesh_opt,'cavity'))
    mesh.nx = 17;
    mesh.ny = 17;
    mesh.nnode = mesh.nx*mesh.ny;
    mesh.Lx = 1.0;
    mesh.Ly = 1.0;
    mesh.filename = 'mesh';
elseif (strcmp(mesh_opt,'poiseuille'))
    mesh.nx = 17;
    mesh.ny = 17;
    mesh.nnode = mesh.nx*mesh.ny;
    mesh.Lx = 1.0;
    mesh.Ly = 1.0;
    mesh.filename = 'mesh';
end
% Solid properties
sprop.rho = 0.1; %1; %0.1;
sprop.nu = 0.35; %0.3; %0.35;
sprop.E  = 2.5e6; %1e8; %2.5e6;
% Fluid properties
fprop.mu = 1.82e-4;
fprop.rho = 1.18e-3;
% Mesh properties
mprop.nu = 0.35;
mprop.alpha = 3; % related to Jacobian-dependent stiffness

% Problem-designated parameters
ndim_out = 2;
temp_sim = 1;

% User-designated parameters
init_var = zeros(1,4); % u(:), y(:)
init_var(1) = 51.3;
init_var(2) = 10; %10; % perturb the y velocity
time.dt = 1e-2;
time.endtime = 3;
time.max_tstep = 5e6;
time.rho = 0.2;
time.outputint = 2e-2;
iter.tol = 1e-2;
iter.max_iter = 10;
plots.plotnow = 0;
c1 = 4;

% Output filename
filename  = 'mesh1_i10_dt01';
text_file = [filename '.txt'];
data_file = [filename '.dat'];
meta_file = [filename '.met'];

%% Body forces
fxmag = 0; % 1e8;
fymag = 0; %-1e4;
phasem = 0; phaseoff = 0.5;
f1 = @(t) fxmag*cos(phasem*pi()*t+phaseoff);
f2 = @(t) fymag*sin(phasem*pi()*t+phaseoff);
% Set up force vector function, 2D
Ff = {f1, f2};  % functions
f = zeros(2,1); % values

%% Constants and Parameters Calculated from Inputs
sprop.lambda = sprop.nu*sprop.E/((1+sprop.nu)*(1-2*sprop.nu));
sprop.mu = sprop.E/(2+2*sprop.nu);

mprop.lambda = mprop.nu/((1+mprop.nu)*(1-2*mprop.nu));
mprop.mu = 1/(2+2*mprop.nu);

time.outputint = time.outputint-eps;
time.alpha_m = (3-time.rho)/(2*(1+time.rho));
time.alpha_f = 1/(1+time.rho);
time.gamma = 0.5+time.alpha_m-time.alpha_f;
time.beta = 1/4*((1+time.alpha_m-time.alpha_f)^2);
plots.disp_nodes = [];

%% Mesh
% Get initial mesh and 
% Designate which nodes are solid or fluid (fluid = 0, solid = 1)
switch (mesh_opt)
    case('meshgen')
        mesh = meshgen_nofile(mesh);
        % mesh.tag(:) = 0; % All fluid
        % mesh.tag(:) = 1; % All solid
    case('uniform')
        mesh = uniform_tri(mesh);
        mesh.tag = zeros(size(mesh.conn,1),1); % All fluid
        mesh.tag([314,42,318,444,307,317,289,76]) = 1; % Obstacle
        mesh.block_elem = [314,42,318,444,307,317,289,76];
        mesh.block_base = [144 145 146];
        mesh.block_intf = [144 146 161 163 178 179 180];
        % mesh.tag = ones(size(mesh.conn,1),1); % All solid
%         % Putting a defect in the mesh
%         mesh.y((3*mesh.nx+1)/2:mesh.nx:((mesh.nx+1)/2+mesh.nx*(mesh.ny-2)))=...
%             mesh.y((3*mesh.nx+1)/2:mesh.nx:((mesh.nx+1)/2+mesh.nx*(mesh.ny-2)))...
%             -mesh.Ly/(mesh.ny-1)*0.5;
    case('bar')
        mesh = uniform_tri(mesh);
        mesh.tag = ones(size(mesh.conn,1),1); % All solid
        mesh.block_base = [1 18 35 52 69];
    case('cavity')
        mesh = uniform_tri(mesh);
        mesh.tag = zeros(size(mesh.conn,1),1); % All fluid
    case('poiseuille')
        mesh = uniform_tri(mesh);
        mesh.tag = zeros(size(mesh.conn,1),1); % All fluid
end

%% Set up Gauss pts (linear triangle, ng = 3)
gauss.ng = 3;
[gauss.ksi_G,gauss.eta_G,gauss.W_G] = tri_Gpts_3(); 

%% Boundary Conditions
% Set boundary conditions
% Locations where mesh is fixed should also be where solid is fixed
% At interface nodes, I should have a pressure BC where nothing is done,
% and at the remaining solid nodes, I'll set the pressure diagonal to 1
switch (mesh_opt)
    case('meshgen')
        nBC = 19;
        BC = cell(nBC,3);
        % cell(number of BCs, BC parameters)
        % BC parameters
        % 1: component of solution designated by BC
        % 2: value of solution at BC
        % 3: node list for BC location
        iBC = 1;
        % Fluid BCs
        % Inlet
        BC{iBC,1} = 'dxf'; BC{iBC,2} = init_var(1); BC{iBC,3} = mesh.inlet; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0   ; BC{iBC,3} = mesh.inlet; iBC=iBC+1;
        % Obstacle
        BC{iBC,1} = 'dxf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.wall; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.wall; iBC=iBC+1;
        % Top and bottom
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.top; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bottom; iBC=iBC+1;
        
        % Fix mesh at edges of domain and at obstacle
        % Inlet
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.inlet; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.inlet; iBC=iBC+1;
        % Obstacle
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.wall; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.wall; iBC=iBC+1;
        % Top
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.top; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.top; iBC=iBC+1;
        % Bottom
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bottom; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bottom; iBC=iBC+1;
        % Outlet
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.outlet; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.outlet; iBC=iBC+1;
        
        % Fix solid at base of bar
        BC{iBC,1} = 'dxs'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bar_base; iBC=iBC+1;
        BC{iBC,1} = 'dys'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bar_base; iBC=iBC+1;
        % Matching BC at bar
        BC{iBC,1} = 'm';   BC{iBC,2} = mesh.bar_intf_elem;
        BC{iBC,3} = mesh.bar_intf; iBC=iBC+1;
        
    case('cavity')
        nBC = 16;
        BC = cell(nBC,3);
        iBC = 1;
        % Left
        BC{iBC,1} = 'dxf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.lnodes; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.lnodes; iBC=iBC+1;
        % Top and bottom
        BC{iBC,1} = 'dxf'; BC{iBC,2} = 100*fprop.mu/fprop.rho; BC{iBC,3} = mesh.tnodes; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.tnodes; iBC=iBC+1;
        BC{iBC,1} = 'dxf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        % Right
        BC{iBC,1} = 'dxf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.rnodes; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.rnodes; iBC=iBC+1;
        
        % Fix mesh at boundaries
        % Lefts
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.lnodes; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.lnodes; iBC=iBC+1;
        % Top and bottom
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.tnodes; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.tnodes; iBC=iBC+1;
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        % Right
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;

    case('poiseuille')
        nBC = 14;
        BC = cell(nBC,3);
        iBC = 1;
        % Left
        BC{iBC,1} = 'dxf'; BC{iBC,2} = 100*fprop.mu/fprop.rho; BC{iBC,3} = mesh.lnodes; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.lnodes; iBC=iBC+1;
        % Top and bottom
        BC{iBC,1} = 'dxf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.tnodes; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.tnodes; iBC=iBC+1;
        BC{iBC,1} = 'dxf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        
        % Fix mesh at boundaries
        % Lefts
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.lnodes; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.lnodes; iBC=iBC+1;
        % Top and bottom
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.tnodes; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.tnodes; iBC=iBC+1;
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        % Right
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        
    case('uniform')
        nBC = 17;
        BC = cell(nBC,3);
        iBC = 1;
        % Inlet
        BC{iBC,1} = 'dxf'; BC{iBC,2} = init_var(1); BC{iBC,3} = mesh.lnodes; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.lnodes; iBC=iBC+1;
        % Top and bottom
        BC{iBC,1} = 'dxf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.tnodes; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.tnodes; iBC=iBC+1;
        BC{iBC,1} = 'dxf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        % Base of obstacle
        BC{iBC,1} = 'dxf'; BC{iBC,2} = 0; BC{iBC,3} = [144 161 178]; iBC=iBC+1;
        BC{iBC,1} = 'dyf'; BC{iBC,2} = 0; BC{iBC,3} = [144 161 178]; iBC=iBC+1;
        
        % Fix mesh at boundaries, done by setting 'solid' BCs
        % Inlet
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.lnodes; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.lnodes; iBC=iBC+1;
        % Top and bottom
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.tnodes; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.tnodes; iBC=iBC+1;
        BC{iBC,1} = 'dxm'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        BC{iBC,1} = 'dym'; BC{iBC,2} = 0; BC{iBC,3} = mesh.bnodes; iBC=iBC+1;
        % Fix mesh at base of obstacle
        BC{iBC,1} = 'dxs'; BC{iBC,2} = 0; BC{iBC,3} = mesh.block_base; iBC=iBC+1;
        BC{iBC,1} = 'dys'; BC{iBC,2} = 0; BC{iBC,3} = mesh.block_base; iBC=iBC+1;
        % Obstacle - match (interface nodes only)
        BC{iBC,1} = 'm'; BC{iBC,2} = mesh.block_elem;
        BC{iBC,3} = mesh.block_intf; iBC=iBC+1;
    case('bar')
        nBC = 2;
        BC = cell(nBC,3);
        iBC = 1;
        BC{iBC,1} = 'dxs'; BC{iBC,2} = 0; BC{iBC,3} = mesh.block_base; iBC=iBC+1;
        BC{iBC,1} = 'dys'; BC{iBC,2} = 0; BC{iBC,3} = mesh.block_base; iBC=iBC+1;
        plots.disp_nodes = zeros(1,mesh.ny);
        for nn = 1:mesh.ny
            plots.disp_nodes(nn) = mesh.nx+(nn-1)*mesh.ny;
        end
        
end

%% Initialization routines
plots = plots_init(mesh,plots,ndim_out);

%% Coefficients for local routines
coeff.rho_f = fprop.rho;
coeff.mu_f = fprop.mu;
coeff.rho_s = sprop.rho;
coeff.mu_s  = sprop.mu;
coeff.lambda_s = sprop.lambda;
coeff.mu_m  = mprop.mu;
coeff.ma = mprop.alpha;
coeff.lambda_m = mprop.lambda;
coeff.am = time.alpha_m;
coeff.afgdt = time.alpha_f*time.gamma*time.dt;
coeff.afbdt2= time.alpha_f*time.beta*time.dt^2;
coeff.c1 = c1;
coeff.dt = time.dt;
coeff.f  = f;

%% Run Solver
% Initialize arrays
if (temp_sim)
    % If simulation is temporal, use temporal loop
    % Put variables in structure 'solver'
    solver.BC = BC;
    solver.Ff = Ff;
    solver.sprop = sprop;
    solver.fprop = fprop;
    solver.gauss = gauss;
    solver.iter = iter;
    solver.time = time;
    solver.ndim_out = ndim_out;
    solver.plots = plots;
    solver.coeff = coeff;
    solver.mesh = mesh;
    solver.ftxt = text_file;
    solver.fdat = data_file;
    solver.fmet = meta_file;
    % Initial guess/starting point
    [u0,udot0,p0,yddot0,ydot0,y0] = set_initial_state(mesh,BC,init_var);
    % Initialize global arrays
    gradG0 = zeros(5*mesh.nnode+1,5*mesh.nnode); 
    resG0 = zeros(5*mesh.nnode+1,1);
    % Solve time evolving problem
    [u,udot,p,yddot,ydot,y] = time_loop(gradG0,resG0,...
        u0,udot0,p0,yddot0,ydot0,y0,solver);
else
    % Call solver function
    [u,f,x,y,conn] = main_fem(mesh.filename,Ff,BC,sprop,ndim_out);
    plots_current(mesh,plots,ndim_out,u);
end