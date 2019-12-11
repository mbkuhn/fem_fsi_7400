function plots = plots_init(mesh,plots,ndim_out)
if (round((mesh.nx+1)/2) ~= (mesh.nx+1)/2)
    ix1 = round((mesh.nx+1)/2);
    ix2 = ix1-1;
else
    ix1 = (mesh.nx+1)/2;
    ix2 = ix1;
end
if (round((mesh.ny+1)/2) ~= (mesh.ny+1)/2)
    iy1 = round((mesh.ny+1)/2);
    iy2 = iy1-1;
else
    iy1 = (mesh.ny+1)/2;
    iy2 = iy1;
end
% u_2 as function of x, at ly/2
plots.off_iy1 = ndim_out*mesh.nx*(iy1-1); 
plots.off_iy2 = ndim_out*mesh.nx*(iy2-1);
plots.ux = zeros(mesh.nx,1);
plots.ytemp = mesh.y(1+((1:mesh.ny)-1)*mesh.nx);
plots.uy = zeros(mesh.nx,1);
plots.off_ix1 = ndim_out*ix1; 
plots.off_ix2 = ndim_out*ix2;

% Get connectivity of solid
plots.selem = (find(mesh.tag==1));
plots.sconn = mesh.conn(plots.selem,:);