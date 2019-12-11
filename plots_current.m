function plots_current(mesh,plots,ndim_out,u,p,y)
% % Collect data
% for i=1:mesh.nx
%     plots.ux(i) = (u(plots.off_iy1+ndim_out*i)+u(plots.off_iy2+ndim_out*i))/2;
% end
% % u_2 as function of y, at lx/2
% for i=1:mesh.ny
%     plots.uy(i) = (u(plots.off_ix1+(i-1)*ndim_out*mesh.nx) ...
%         + u(plots.off_ix2+(i-1)*ndim_out*mesh.nx))/2;
% end
% 
% % Plot results along x
% figure(1)
% plot(mesh.x(1:mesh.nx),plots.ux,'k')
% set(gca,'FontSize',14)
% xlabel('x')
% ylabel('magnitude')
% legend('u_2(x,Ly/2)')
% % Plot results along y
% figure(2)
% plot(plots.ytemp,plots.uy,'k')
% set(gca,'FontSize',14)
% xlabel('y')
% ylabel('magnitude')
% legend('u_2(Lx/2,y)')

% Plot results spatially as displacement in x
mf = 1e0; % magnification factor
figure(3)
if (ndim_out == 2)
    triplot(mesh.conn,mesh.x,mesh.y,'Color',[0.7 0.7 0.7]);
    hold on
    triplot(mesh.conn,mesh.x+mf*y(1:2:2*mesh.nnode-1),...
        mesh.y+mf*y(2:2:2*mesh.nnode),'Color','black');
    if (~isempty(plots.sconn))
        triplot(plots.sconn,mesh.x+y(1:2:2*mesh.nnode-1),...
            mesh.y+y(2:2:2*mesh.nnode),'Color','red');
    end
    quiver(mesh.x+y(1:2:2*mesh.nnode-1),...
        mesh.y+y(2:2:2*mesh.nnode),u(1:2:2*mesh.nnode-1),...
        u(2:2:2*mesh.nnode),'Color','blue');
    hold off
    xlim([0,mesh.Lx]); ylim([0,mesh.Ly])
    figure(4)
    trimesh(mesh.conn,mesh.x,mesh.y,p);
    xlim([0,mesh.Lx]); ylim([0,mesh.Ly])
else % (ndim_out == 1)
    trimesh(mesh.conn,mesh.x,mesh.y,u);
    xlim([0,mesh.Lx]); ylim([0,mesh.Ly])
end

drawnow