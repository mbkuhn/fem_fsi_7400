function [u,udot,p,yddot,ydot,y] = ...
    time_loop(gradG,resG,u,udot,p,yddot,ydot,y,solver)

% Unpack structures
time = solver.time;

% Open output files
[fdat,ftxt,fmet]=output_init(solver.fdat,solver.ftxt,solver.fmet,...
    solver.plots.sconn,solver.mesh);

% Get nodes for output
intf_nodes = [];
for iBC=1:size(solver.BC,1)
    if (solver.BC{iBC,1}(1)=='m')
        % Check element list for BC
        intf_nodes = [intf_nodes solver.BC{iBC,3}];
    end
end

time.nstep = 0;
time.ctime = 0.0;
while ((time.nstep<time.max_tstep)&&(time.ctime<time.endtime))
    % Calculate current value of force
    solver.coeff.f(1) = solver.Ff{1}(time.ctime);
    solver.coeff.f(2) = solver.Ff{2}(time.ctime);
    
    % Solve timestep
    time.nstep = time.nstep+1;
    [u,udot,p,yddot,ydot,y,niter,err] = ...
        time_adv(gradG,resG,u,udot,p,yddot,ydot,y,solver);
    time.ctime = time.ctime + time.dt;
    
    % Get stuff for printing
    umag = (u(1:2:length(u)-1).^2+u(2:2:length(u)).^2).^0.5;
    dnodes = [intf_nodes solver.plots.disp_nodes];
    y1intf= y(2*dnodes-1); y2intf=y(2*dnodes);
    if (isempty(y2intf))
        y2intf = 0;
    end
    
    % Print timestep info
    pstr = '%6i Time %9.2e Iter %2i %9.2e Umax %9.2e Pmax %9.2e Ymax %9.2e Ymin %9.2e \n';
    pnum = [time.nstep,time.ctime,niter,err,max(umag),max(p),max(y2intf),min(y2intf)];
    fprintf( pstr, pnum);
    fprintf(ftxt,pstr,pnum);
    % Plot instant
    if (mod(time.ctime+0.5*time.dt,time.outputint) < time.dt)
        output_data(fdat,fmet,pstr,pnum,...
            time.nstep,time.ctime,u,p,y);
        if (solver.plots.plotnow)
            plots_current(solver.mesh,solver.plots,solver.ndim_out,u,p,y);
        end
    end
end

% Close output file
output_close(fdat,ftxt,fmet);