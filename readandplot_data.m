clf;
read_flag = 1;

if (read_flag)
    clear; clc;

    fname = 'coarse_mesh_fsi';
    
    % Filenames
    fname_dat = [fname,'.dat'];
    fname_txt = [fname,'.txt'];
    fname_met = [fname,'.met'];
    
    % File IDs
    fdat = fopen(fname_dat,'r');
    ftxt = fopen(fname_txt,'r');
    fmet = fopen(fname_met,'r');
    
    % Read from text file
    A = fscanf(ftxt,'%d Time %e Iter %d %e Umax %e Pmax %e Ymax %e Ymin %e',[8 Inf]);
    tall = A(2,:); 
    ymax = A(7,:);
    ymin = A(8,:);
    
    % Read from metadata file
    B = fscanf(fmet,'%d %e',[2 Inf]);
    nout = max(B(1,:));
    tout = B(2,:);
    
    % Read size of domain (nnodes, nelem, nselem)
    nnode = fread(fdat,1,'int');
    nelem = fread(fdat,1,'int');
    nselem= fread(fdat,1,'int');
    
    % Read mesh and connectivity
    mx = fread(fdat,nnode,'double');
    my = fread(fdat,nnode,'double');
    conn = fread(fdat,[nelem  3],'double');
    sconn= fread(fdat,[nselem 3],'double');
    
    % Loop output timesteps and save data
    u = zeros(2*nnode,nout); y = u;
    p = zeros(nnode,nout);
    for ntime = 1:nout
        % Save temporal data
        u(:,ntime) = fread(fdat,2*nnode,'double');
        p(:,ntime) = fread(fdat,  nnode,'double');
        y(:,ntime) = fread(fdat,2*nnode,'double');
    end
    
    % Close files
    fclose(fdat);
    fclose(ftxt);
    fclose(fmet);
end

% Plot what I want to plot

np = 2;
figure(3);
triplot(conn,mx,my,'Color',[0.7 0.7 0.7]);
hold on
triplot(conn,mx+y(1:2:2*nnode-1,np),...
    my+y(2:2:2*nnode,np),'Color','black');
triplot(sconn,mx+y(1:2:2*nnode-1,np),...
    my+y(2:2:2*nnode,np),'Color','red');
quiver(mx+y(1:2:2*nnode-1,np),...
    my+y(2:2:2*nnode,np),u(1:2:2*nnode-1,np),...
    u(2:2:2*nnode,np),'Color','blue');
hold off
axis equal
% xlim([0,mesh.Lx]); ylim([0,mesh.Ly])
figure(4)
trimesh(conn,mx,my,p(:,np));

% Save to video