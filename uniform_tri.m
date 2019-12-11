function mesh = uniform_tri(mesh)
% Input : (nx,ny) = number of nodes in (x,y)
%         (Lx,Ly) = length of domain in each direction
% Output: conn = connectivity of points
%         x    = vector of x positions of points
%         y    = vector of y positions of points

% Get variables from structure
nx = mesh.nx; ny = mesh.ny;
Lx = mesh.Lx; Ly = mesh.Ly;
mesh_file = mesh.filename;

x = zeros(nx*ny,1);
ytemp = 0:Ly/(ny-1):Ly;
y = zeros(nx*ny,1);
for j=1:ny
    x(1+(j-1)*nx:nx+(j-1)*nx) = 0:Lx/(nx-1):Lx;
    y(1+(j-1)*nx:nx+(j-1)*nx) = ytemp(j);
end
conn = delaunay(x,y);
% Get number of elements and number of nodes
[nelem,npe] = size(conn);

% Write information to file
fid = fopen(mesh_file,'w');
fwrite(fid,nx,'integer*4');
fwrite(fid,ny,'integer*4');
fwrite(fid,nelem,'integer*4');
fwrite(fid,npe,'integer*4');
fwrite(fid,x,'double');
fwrite(fid,y,'double');
fwrite(fid,conn,'double');
fclose(fid);

% Make node lists for BCs
mesh.lnodes = 1:nx:(1+nx*(ny-1));
mesh.rnodes = mesh.lnodes+(nx-1);
mesh.bnodes = 1:nx;
mesh.tnodes = mesh.bnodes+nx*(ny-1);

% Keep some variables in structure
mesh.nelem = nelem;
mesh.npe = npe;
mesh.conn = conn;
mesh.x = x;
mesh.y = y;
end