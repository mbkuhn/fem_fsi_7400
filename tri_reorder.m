function [econn,xl,yl] = tri_reorder(x,y,econn)
% Function to reorder the triangle points such that pt 3 is opposite the
% longest side of the triangle. Also outputs the coordinates of the nodes.

% Inputs:  econn = connectivity of a single element (1D array)
%          x = all x coordinates
%          y = all y coordinates
% Outputs: econn = reordered connectivity of a single element
%          xl = x coordinates of element nodes
%          yl = y coordinates of element nodes

xl= zeros(3,1); yl= zeros(3,1); dl= zeros(3,1);
for nn=1:3
    % Get coordinates
    xl(nn) = x(econn(nn));
    yl(nn) = y(econn(nn));
    % Calculate distance
    if (nn~=1)
        dl(nn-1) = sqrt((xl(nn)-xl(nn-1))^2+(yl(nn)-yl(nn-1))^2);
    end
end
dl(3) = sqrt((xl(3)-xl(1))^2+(yl(3)-yl(1))^2);
% Get index of maximum distance between nodes
[~,imax]=max(dl);
% Node opposite of max distance needs to be point 3
if (imax == 2) % pt 1 should be pt 3
    econn = [econn(2:3),econn(1)];
    xl = [xl(2:3);xl(1)];
    yl = [yl(2:3);yl(1)];
elseif (imax == 3) % pt 2 should be pt 3
    econn(2:3) = [econn(3),econn(2)];
    xl(2:3) = [xl(3); xl(2)];
    yl(2:3) = [yl(3); yl(2)];
end
% Make node ordering counterclockwise to get positive Jacobian
% (if Jacobian is negative, switch first two nodes)
if ((xl(1)-xl(3))*(yl(2)-yl(3)) < (yl(1)-yl(3))*(xl(2)-xl(3)))
    econn(1:2) = [econn(2),econn(1)];
    xl(1:2) = [xl(2),xl(1)];
    yl(1:2) = [yl(2),yl(1)];
end
end