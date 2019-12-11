function [x,y,conn,nx,ny,nelem,npe] = read_mesh(mesh_file)
% Read information from file
fid = fopen(mesh_file,'r');
nx = fread(fid,1,'integer*4');
ny = fread(fid,1,'integer*4');
nelem = fread(fid,1,'integer*4');
npe = fread(fid,1,'integer*4');
x = fread(fid,[nx*ny,1],'double');
y = fread(fid,[nx*ny,1],'double');
conn = fread(fid,[nelem,npe],'double');
fclose(fid);
end