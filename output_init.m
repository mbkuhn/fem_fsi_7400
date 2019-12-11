function [fdat,ftxt,fmet] = output_init(namedat,nametxt,namemet,sconn,mesh)

% Open specified files
fdat = fopen(namedat,'w');
ftxt = fopen(nametxt,'w');
fmet = fopen(namemet,'w');

% Write size of domain (nnode, nelem, nselem)
fwrite(fdat,mesh.nnode,'int');
fwrite(fdat,mesh.nelem,'int');
fwrite(fdat,size(sconn,1),'int'); % number of elements in solid

% Write mesh and connectivity
fwrite(fdat,mesh.x,'double');
fwrite(fdat,mesh.y,'double');
fwrite(fdat,mesh.conn,'double');
fwrite(fdat,sconn,'double');

end