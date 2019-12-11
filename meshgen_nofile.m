function mesh = meshgen_nofile(mesh)
% Total number of grid points in x and y direction
nx=13; %51; 13 21 31 41
ny=7; %51;  7 11 25 41

% Number of grid points for the bar
bnx=7; % 7 11 25 41 % This includes the obstacle
bny=3;

% Size of the domain
l=19.5;
h=12;

% Size of the obstacle
ol=1;
oh=ol;

% Size of the bar
bl=4;
bh=0.06;

% Size of the entrance
el=4.5;
enx=ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of nodes on the obstacle should be an integer
onx = round((bnx-1)*ol/(ol+bl)) + 1;
bnx = (onx-1)*(1+bl/ol) + 1;

ny = round((ny-bny)/2)*2 + bny;

% Total number of nodes
nNo = (bnx-2)*enx*2 + (2*enx+nx-onx)*ny - 2*(bnx-onx);
% Node coordinates
X=zeros(nNo,2);

% Boundary nodes list
wall=zeros(2*(ny+onx-2),1);
inlet=zeros(ny,1);
outlet=zeros(ny-2,1);
top=zeros(nx,1);
bottom=zeros(nx,1);
bar=zeros((bnx-onx)*bny,1);
bar_intf=zeros((bnx-onx-1)*2+bny,1);
bar_base=zeros(bny,1);

% Meshing the entrance
a=0;
b=0;
c=0;
for i=1:enx
   for j=2:ny-1
      a = a + 1;
      X(a,1) = (i-1)*el/(enx-1);
      X(a,2) = (j-1)/(ny-1)*(i-1)/(enx-1)*ol ...
             + (j-1)/(ny-1)*(enx-i)/(enx-1)*h...
             + (i-1)/(enx-1)*(h-ol)/2;
      if (i==1) 
        b = b + 1;
        inlet(b) = a;
      end
      if (i==enx)
        c = c+1;
        wall(c) = a;
      end
   end
end

% Meshing the region above and below the obstacle
% Entrance and wake length
wl = 2*el + ol + bl;
b=0;
for i=1:bnx
   for j=1:enx
      a = a + 1;
      X(a,1) = (i-1)/(bnx-1)*(j-1)/(enx-1)*(bl+ol) ...
             + (i-1)/(bnx-1)*(enx-j)/(enx-1)*wl ...
             + (j-1)/(enx-1)*el;
      X(a,2) = (j-1)/(enx-1)*(h-oh)/2;
      a = a + 1;
      X(a,1) = X(a-1,1);
      X(a,2) = h - X(a-1,2);
      if (j==1) 
        b = b + 1;
        top(b) = a;
        bottom(b) = a-1;
      end
      if (j==enx && i <= onx) 
        c = c + 2;
        wall(c-1:c) = [a-1 a];
      end
   end
end

% Meshing the region above and below the bar
for i=1:bnx-onx+1
   for j=2:(ny-bny)/2
      a = a + 1;
      X(a,1) = (i-1)/(bnx-onx)*bl + el + ol;
      X(a,2) = (j-1)/(ny-bny)*(ol-bh) + (h-ol)/2;
      a = a + 1;
      X(a,1) = X(a-1,1);
      X(a,2) = h - X(a-1,2);
      if (i == 1) 
        c = c + 2;
        wall(c-1:c) = [a-1 a];
      end
   end
end

% Meshing the bar
d=0; f=0; ff=0;
for i=1:bnx-onx+1
   for j=1:bny
      a = a + 1;
      X(a,1) = (i-1)/(bnx-onx)*bl + el + ol;
      X(a,2) = (j-1)/(bny - 1)*bh + (h-bh)/2;
      d = d + 1;
      bar(d) = a;
      if (i == 1) 
        c = c + 1;
        wall(c) = a;
      end
      if ( (i==1 && (j==1 || j==bny) ) || ( i > 1 && ( i==bnx-onx+1 || ...
              j==1 || j==bny ) ) )
          % Save edge of bar, no internal points
          f = f + 1;
          bar_intf(f) = a;
      end
      if ( i==1 )
          % Save base of bar, where it is fixed
          ff = ff + 1;
          bar_base(ff) = a;
      end
   end
end

% Meshing the top and bottom part of the wake
for i=2:enx
   for j=2:(ny-bny)/2 + 1
      a = a + 1;
      X(a,1) = (i-1)*el/(enx-1) + el + ol + bl;
      X(a,2) = (j-1)/(ny-bny)*(enx-i)/(enx-1)*(ol-bh) ...
             + (j-1)/(ny-1)*(i-1)/(enx-1)*h...
             + (enx-i)/(enx-1)*(h-oh)/2;
      a = a + 1;
      X(a,1) = X(a-1,1);
      X(a,2) = h - X(a-1,2);
   end
end

% Meshing the wake
for i=2:enx
   for j=2:bny - 1
      a = a + 1;
      X(a,1) = (i-1)*el/(enx-1) + el + ol + bl;
      X(a,2) = (j-1)/(bny-1)*(enx-i)/(enx-1)*bh ...
             + (j-1)/(ny-1)*(i-1)/(enx-1)*h ...
             + (enx-i)/(enx-1)*(h-bh)/2 ...
             + (i-1)/(enx-1)*h/2*(ny-bny)/(ny-1);
   end
end

% Meshing the exit and save outlet nodes
c=b;
f=0;
for i=2:nx-bnx+1
   for j=1:ny
      a = a + 1;
      X(a,1) = (i-1)/(nx-bnx)*(l-wl) + wl;
      X(a,2) = (j-1)/(ny-1)*h;
      if (j==1) 
        b = b + 1;
        bottom(b) = a;
      elseif(j==ny)
        c = c + 1;  
        top(b) = a;
      end
      if (i==nx-bnx+1) 
        f = f + 1;
        outlet(f) = a;
      end
   end
end

% Include first top and first bottom point in inlet
inlet(ny-1)=top(1); inlet(ny) = bottom(1);

% Constructing the connectivity array
IEN = delaunay(X(:,1),X(:,2));

% Removing elements inside the obstacle
tmp=zeros(nNo,1);
tmp(wall)=1;
rml=[];
for e=1:size(IEN,1)
  if (all(tmp(IEN(e,:))==1))
    rml = [rml; e];
  end
end
IEN(rml,:)=[];
nEl = size(IEN,1);

% Filling tag (flagging elements in the bar)
tag = zeros(nEl,1);
for e=1:nEl
    wNo=find(bar==IEN(e,1) | bar==IEN(e,2) | bar==IEN(e,3));
    if length(wNo)==3
        tag(e) = 1;
      else
        tag(e) = 0;
    end
end

% Get elements that are connected to bar interface
bar_intf_elem = [];
for e=1:nEl
    wNo=find(bar_intf==IEN(e,1) | bar_intf==IEN(e,2) | bar_intf==IEN(e,3));
    if length(wNo) > 0
        bar_intf_elem = [bar_intf_elem, e];
    end
end

% Put information into structure
mesh.inlet = inlet;
mesh.top = top;
mesh.bottom = bottom;
mesh.wall = wall;
mesh.bar = bar;
mesh.outlet = outlet;
mesh.bar_intf = bar_intf;
mesh.bar_base = bar_base;
mesh.bar_intf_elem = bar_intf_elem;
mesh.conn = IEN;
mesh.x = X(:,1);
mesh.y = X(:,2);
mesh.nelem = nEl;
mesh.npe = size(IEN,2);
mesh.Lx = l;
mesh.Ly = h;
mesh.nx = nx;
mesh.ny = ny;
mesh.tag = tag;
mesh.nnode = nNo;

% % Plotting the grid
% figure(1)
% plot(X(1:a,1),X(1:a,2),'.',X(bar,1),X(bar,2),'r.',...
%      X(wall,1),X(wall,2),'k.',X(inlet,1),X(inlet,2),'g.',...
%      X(top,1),X(top,2),'c.',X(bottom,1),X(bottom,2),'m.')
% axis equal
% 
% figure(2)
% triplot (IEN, X(:,1), X(:,2));
