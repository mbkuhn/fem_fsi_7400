function [u,udot,p,yddot,ydot,y,niter,err] = ...
    time_adv(gradG,resG,u,udot,p,yddot,ydot,y,solver)
% Save solution from last timestep
oldudot = udot;
oldu    = u;
oldyddot= yddot;
oldydot = ydot;
oldy    = y;

% Unpack variables from structures
gamma = solver.time.gamma;
dt = solver.time.dt;
alpha_f = solver.time.alpha_f;
alpha_m = solver.time.alpha_m;
beta = solver.time.beta;
tol = solver.iter.tol;
max_iter = solver.iter.max_iter;

% Predict solution at n+1
udot = (gamma-1)/gamma*oldudot;
yddot= (gamma-1)/gamma*oldyddot;
y    = oldy + dt*oldydot + 0.5*dt^2*oldyddot + beta*dt^2*(yddot-oldyddot);
% No alteration for predictor of ydot and u

% Initialize err and iter
err = tol+1; niter = 0;
while (err > tol && niter < max_iter)
    % Count iteration
    niter = niter + 1;
    
    % Update intermediate values of u and udot
    udot_am = oldudot + alpha_m*(udot-oldudot);
    u_af = oldu + alpha_f*(u-oldu);

    % Update intermediate values of y, ydot, yddot
    yddot_am = oldyddot + alpha_m*(yddot-oldyddot);
    ydot_af  = oldydot  + alpha_f*(ydot -oldydot );
    y_af     = oldy     + alpha_f*(y    -oldy    );
    
    % Solve timestep
    [du,dp,dy,err] = main_fem_step(gradG,resG,udot_am,u_af,p,...
        yddot_am,ydot_af,y_af,solver);
    
    % Converge based on relative error
    if (niter==1)
        err0 = max(err,1);
        % Don't use relative error if not much is changing
    end
    err = err/err0;
        
    % Update future values of u and udot
    udot = udot + du;
    u = u + gamma*dt*du;
    p = p + dp;
    % Update future values of y, ydot, yddot
    yddot = yddot + dy;
    ydot  = ydot  + gamma*dt*dy;
    y     = y     + beta*dt^2*dy;
    
end