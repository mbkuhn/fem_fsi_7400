clf;
read_flag = 1;

if (read_flag)
    clear; clc;

    fname = 'mesh3_i20_dt005';
    
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
    
    xx = [mx(conn(:,[1 2 3 1])'); NaN(1,size(conn,1))];
    yy = [my(conn(:,[1 2 3 1])'); NaN(1,size(conn,1))];
    
    %npause = 200;
    % Loop output timesteps and save data
    v = VideoWriter([fname '.mp4']);
    open(v)
    for ntime = 1:nout
        
        % Save temporal data
        u = fread(fdat,2*nnode,'double');
        p = fread(fdat,  nnode,'double');
        y = fread(fdat,2*nnode,'double');
        
        %if (ntime==npause)
        plot(xx,yy)
        hold on
        % Plot pressure contour
        [CS,h]=tricontf(mx,my,conn,p);
        set(h,'edgecolor','none');
        [CS,h]=tricont(mx,my,conn,p,'-k');
        clabel(CS,h,'fontsize',14)
        % Plot velocity arrows
        quiver(mx+y(1:2:2*nnode-1),...
            my+y(2:2:2*nnode),u(1:2:2*nnode-1),...
            u(2:2:2*nnode),'Color','white');
        % Plot mesh of bar
        triplot(sconn,mx+y(1:2:2*nnode-1),...
            my+y(2:2:2*nnode),'Color','red');
        axis equal
        
        hold off
        ntime
            %pause()
        %end
        
        
        frame = getframe(gcf);
        writeVideo(v,frame)
    end
    
    % Close files
    %close(v);
    fclose(fdat);
    fclose(ftxt);
    fclose(fmet);
end