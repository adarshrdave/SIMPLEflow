%% F18 CFD Project - incompressible flow around wing
%% Setup
%                     N
% (1,ny) ----------------------------- (nx,ny)
%      |                             |
%      |              2              |
%   W  |    1     |-------|   4      |   E
%      |          |-------|          |
%      |   (w.idx,w.idy)             |
%      |               3             |
%(1,1) ----------------------------- (nx,1)
%                     S

clear all;
% grid parameters
grid.Lx = 200;    %x-length of box
grid.Ly = 100;    %y-length of box
grid.dx = 1;   %cell dimensions
grid.dy = grid.dx;

%number of int. cells: nx*ny - ie. interior p nodes
grid.nx = grid.Lx/grid.dx;
grid.ny = grid.Ly/grid.dx;

%set up arrays for p,u,v - big arrays including ghost nodes
field.p = zeros(grid.nx+2,grid.ny+2);
field.u = zeros(grid.nx+1,grid.ny+2);
field.v = zeros(grid.nx+2,grid.ny+1);


% Inflow across all of W & S
param.alpha = pi/4;
param.vIN = 10; % magnitude of velocity of inflow

% Wing - 
    %desired size in 'm'
w.Ly = 10;
w.Lx = 40;
    %size in no cells
w.ldy = floor(w.Ly/grid.dy);
w.ldx = floor(w.Lx/grid.dx);
    %position in  grid - left-bottom corner-cell of wing
w.idy = ceil((grid.ny/2)-(w.ldy/2));
w.idx = ceil((grid.nx/2)-(w.ldx/2));

% Init BCs
field = setBC(field,grid,param,w);

% Parameters (arbitraryly chosen for now)

% Initialization, Parameters
param.Re = 100;
param.rho = 1.184;
param.mu = 1;
param.T = 5;
param.dt = .02;
param.tsteps = param.T/param.dt;
param.eps = .01;
param.Q=0.5;
param.s=10; %maximum speed in m/s

%% FSM

for t=1:param.tsteps
    
field = setBC(field,grid,param,w);

%figure(1)
%surface(field.u.')
%pause(.01)

%Intermediate velocity field
field = computeVF_sections(field,grid,param,w);

%field = setBC(field,grid,param,w);

field.p = solveP_GS(field,grid,param,w);

field = updateV(field,grid,param,w);

end
visualize(field,grid);
%% test space

%% Functions
%computes CD at v(i,j) node
function[dP_dy] =  dPdy(field,grid)
dP_dy = zeros(size(field.v));
for i=2:grid.nx+1
    for j=2:grid.ny
        dP_dy(i,j) = (field.p(i,j+1)-field.p(i,j))/grid.dy;
    end
end
end

%computes CD at u(i,j)
function[dP_dx] =  dPdx(field,grid)
dP_dx = zeros(size(field.u));
for i=2:grid.nx
    for j=2:grid.ny+1
        dP_dx(i,j) = (field.p(i+1,j)-field.p(i,j))/grid.dx;
    end
end
end

function[field] = setBC(field,grid,param,w)
%Corners
field.p(1,1) = 0;
field.p(grid.nx+2,1) = 0;
field.p(1,grid.ny+2) = 0;
field.p(grid.nx+2,grid.ny+2) = 0;

field.v(1,1) = 0;
field.v(grid.nx+2,1) = 0;
field.v(1,grid.ny+1) = 0;
field.v(grid.nx+2,grid.ny+1) = 0;

field.u(1,1) = 0;
field.u(grid.nx+1,1) = 0;
field.u(1,grid.ny+2) = 0;
field.u(grid.nx+1,grid.ny+2) = 0;

%West - Neumann
field.u(1,2:grid.ny+1) = param.vIN*sin(param.alpha);
field.v(1,2:grid.ny) = 2*param.vIN*cos(param.alpha) - field.v(2,2:grid.ny);
field.p(1,2:grid.ny+1) = field.p(2,2:grid.ny+1);

%South - Neumann
field.u(2:grid.nx,1) = 2*param.vIN*sin(param.alpha)-field.u(2:grid.nx,2);
field.v(2:grid.nx+1,1) = param.vIN*cos(param.alpha);
field.p(2:grid.nx+1,1) = field.p(2:grid.nx+1,2);

%North - far-field
%field.v(2:grid.nx+1,grid.ny+1) = 0;
field.v(2:grid.nx+1,grid.ny+1) = field.v(2:grid.nx+1,grid.ny);
field.u(2:grid.nx,grid.ny+2) = field.u(2:grid.nx,grid.ny+1);
field.p(2:grid.nx+1,grid.ny+2) = field.p(2:grid.nx+1,grid.ny+1);

%East - far-field
field.v(grid.nx+2,2:grid.ny) = field.v(grid.nx+1,2:grid.ny);
%field.u(grid.nx+1,2:grid.ny+1) = 0;
field.u(grid.nx+1,2:grid.ny+1) = field.u(grid.nx,2:grid.ny+1);
field.p(grid.nx+2,2:grid.ny+1) = field.p(grid.nx+1,2:grid.ny+1);

%BC for wing - 
%bot
field.v(w.idx:w.idx+w.ldx-1,w.idy-1)=0;
field.u(w.idx:w.idx+w.ldx-2,w.idy)= ...
    -field.u(w.idx:w.idx+w.ldx-2,w.idy-1);
field.p(w.idx:w.idx+w.ldx-1,w.idy)=field.p(w.idx:w.idx+w.ldx-1,w.idy-1);
%top
field.v(w.idx:w.idx+w.ldx-1,w.idy+w.ldy-1)=0;
field.u(w.idx:w.idx+w.ldx-2,w.idy+w.ldy-1)= ...
    -field.u(w.idx:w.idx+w.ldx-2,w.idy+w.ldy);
field.p(w.idx:w.idx+w.ldx-1,w.idy+w.ldy-1)=field.p(w.idx:w.idx+w.ldx-1,w.idy+w.ldy);

%left 
field.u(w.idx-1,w.idy:w.idy+w.ldy-1)=0;
field.v(w.idx,w.idy:w.idy+w.ldy-2)=...
    -field.v(w.idx-1,w.idy:w.idy+w.ldy-2);
field.p(w.idx,w.idy:w.idy+w.ldy-1) = field.p(w.idx-1,w.idy:w.idy+w.ldy-1);

%right
field.u(w.idx+w.ldx-1,w.idy:w.idy+w.ldy-1)=0;
field.v(w.idx+w.ldx-1,w.idy:w.idy+w.ldy-2)=...
    -field.v(w.idx+w.ldx,w.idy:w.idy+w.ldy-2);
field.p(w.idx+w.ldx-1,w.idy:w.idy+w.ldy-1) = ...
    field.p(w.idx+w.ldx,w.idy:w.idy+w.ldy-1);
end

function[field] = computeVF_sections(field,grid,param,w)
U=field.u;
V=field.v;

%U loops
%part 1
for i=2:w.idx-2
    for j=2:grid.ny+1    
        field.u(i,j)= computeAB_u(grid,param,U,V,i,j);
    end
end
%part 2 - top
for i=w.idx-1:w.idx+w.ldx-1
    for j=w.idy+w.ldy:grid.ny+1
        field.u(i,j) = computeAB_u(grid,param,U,V,i,j);
    end
end
%part 3 - bot
for i=w.idx-1:w.idx+w.ldx-1
    for j=2:w.idy-1
    field.u(i,j)= computeAB_u(grid,param,U,V,i,j);
    end
end
%part 4 - right
for i=w.idx+w.ldx:grid.nx
    for j=2:grid.ny+1
     field.u(i,j)= computeAB_u(grid,param,U,V,i,j);
    end
end

%V loops 
for i=2:w.idx-1
    for j=2:grid.ny
        field.v(i,j)= computeAB_v(grid,param,U,V,i,j);
    end
end
for i=w.idx:w.idx+w.ldx-1
    for j=w.idy+w.ldy:grid.ny
    field.v(i,j)= computeAB_v(grid,param,U,V,i,j);
    end
end
for i=w.idx:w.idx+w.ldx-1
    for j=2:w.idy-2
      field.v(i,j)= computeAB_v(grid,param,U,V,i,j);
    end
end
for i=w.idx+w.ldx:grid.nx+1
    for j=2:grid.ny
     field.v(i,j)= computeAB_v(grid,param,U,V,i,j);
    end
end  
end

function[u] = computeAB_u(grid,param,U,V,i,j)
 uuxCD=((U(i,j)+U(i+1,j))^2-(U(i,j)+U(i-1,j))^2)/(4*grid.dx);
        uvyCD=((V(i,j)+V(i+1,j))*(U(i,j)+U(i,j+1))-(V(i,j-1)+V(i+1,j-1))*(U(i,j-1)+U(i,j)))/(4*grid.dy);
        %Upwind
        uuxUW=(abs(U(i,j)+U(i+1,j))*(U(i,j)-U(i+1,j))-abs(U(i-1,j)+U(i,j))*(U(i-1,j)-U(i,j)))*param.Q/(4*grid.dx);
        uvyUW=(abs(V(i,j)+V(i+1,j))*(U(i,j)-U(i,j+1))-abs(V(i,j-1)+V(i+1,j-1))*(U(i,j-1)-U(i,j)))*param.Q/(4*grid.dx);     
        %Diffusion
        uxx=(U(i+1,j)-2*U(i,j)+U(i-1,j))/(grid.dx*grid.dx);
        uyy=(U(i,j+1)-2*U(i,j)+U(i,j-1))/(grid.dy*grid.dy);
        
        %Update
        u=U(i,j)+param.dt*((1/param.Re)*(uxx+uyy)-(uuxCD+uuxUW+uvyCD+uvyUW));
end
function[v] = computeAB_v(grid,param,U,V,i,j)
%Centered Difference
        uvxCD=((U(i,j+1)+U(i,j))*(V(i,j)+V(i+1,j))-(U(i-1,j+1)+U(i-1,j))*(V(i,j)+V(i-1,j)))/(4*grid.dx);
        vvyCD=((V(i,j)+V(i,j+1))^2-(V(i,j)+V(i,j-1))^2)/(4*grid.dy);
        %Upwind
        uvxUW=(abs(U(i,j)+U(i,j+1))*(V(i,j)-V(i+1,j))-abs(U(i-1,j)+U(i-1,j+1))*(V(i-1,j)-V(i,j)))*param.Q/(4*grid.dx);
        vvyUW=(abs(V(i,j)+V(i,j+1))*(V(i,j)-V(i,j+1))-abs(V(i,j-1)+V(i,j))*(V(i,j-1)-V(i,j)))*param.Q/(4*grid.dy);
        %Diffusion
        vxx=(V(i+1,j)-2*V(i,j)+V(i-1,j))/(grid.dx*grid.dx);
        vyy=(V(i,j+1)-2*V(i,j)+V(i,j-1))/(grid.dy*grid.dy);
        
        %update
        v=V(i,j)+param.dt*((1/param.Re)*(vxx+vyy)-(uvxCD+uvxUW+vvyCD+vvyUW));
end

function[p_new] = solveP_GS(field,grid,param,w)
p_new = field.p;

for i=2:w.idx-1
    for j=2:grid.ny+1
         dUdx = (field.u(i,j) - field.u(i-1,j))/grid.dx;
        dVdy = (field.v(i,j) - field.v(i,j-1))/grid.dy;
        p_new(i,j) = (1/4)*(p_new(i-1,j) + field.p(i+1,j) + p_new(i,j-1) + ... 
            field.p(i,j+1) - (grid.dx^2/param.dt)*(dUdx + dVdy));
    end
end
for i=w.idx:w.idx+w.ldx-1
    for j=w.idy+w.ldy:grid.ny+1
          dUdx = (field.u(i,j) - field.u(i-1,j))/grid.dx;
        dVdy = (field.v(i,j) - field.v(i,j-1))/grid.dy;
        p_new(i,j) = (1/4)*(p_new(i-1,j) + field.p(i+1,j) + p_new(i,j-1) + ... 
            field.p(i,j+1) - (grid.dx^2/param.dt)*(dUdx + dVdy));
    end
end
for i=w.idx:w.idx+w.ldx-1
  for j=2:w.idy-1
       dUdx = (field.u(i,j) - field.u(i-1,j))/grid.dx;
        dVdy = (field.v(i,j) - field.v(i,j-1))/grid.dy;
        p_new(i,j) = (1/4)*(p_new(i-1,j) + field.p(i+1,j) + p_new(i,j-1) + ... 
            field.p(i,j+1) - (grid.dx^2/param.dt)*(dUdx + dVdy));
  end
end
for i=w.idx+w.ldx:grid.nx+1
    for j=2:grid.ny+1
       dUdx = (field.u(i,j) - field.u(i-1,j))/grid.dx;
        dVdy = (field.v(i,j) - field.v(i,j-1))/grid.dy;
        p_new(i,j) = (1/4)*(p_new(i-1,j) + field.p(i+1,j) + p_new(i,j-1) + ... 
            field.p(i,j+1) - (grid.dx^2/param.dt)*(dUdx + dVdy));
    end
end

end

function[field] = updateV(field, grid, param,w)

dP_dy = dPdy(field,grid);
dP_dx = dPdx(field,grid);

%U loops- for everything computed at u(i,j)
%Part 1 - left
for i=2:w.idx-2
    for j=2:grid.ny+1   
        field.u(i,j) = field.u(i,j)-param.dt*dP_dx(i,j);
    end
end
%part 2 - top
for i=w.idx-1:w.idx+w.ldx-1
    for j=w.idy+w.ldy:grid.ny+1
        field.u(i,j) = field.u(i,j)-param.dt*dP_dx(i,j);
    end
end
%part 3 - bot
for i=w.idx-1:w.idx+w.ldx-1
    for j=2:w.idy-1
     field.u(i,j) = field.u(i,j)-param.dt*dP_dx(i,j);
    end
end
%part 4 - right
for i=w.idx+w.ldx:grid.nx
    for j=2:grid.ny+1
       field.u(i,j) = field.u(i,j)-param.dt*dP_dx(i,j);
    end
end

%V loops - for everything computed at V(i,j)
for i=2:w.idx-1
    for j=2:grid.ny
        field.v(i,j) = field.v(i,j)-param.dt*dP_dy(i,j);
    end
end
for i=w.idx:w.idx+w.ldx-1
    for j=w.idy+w.ldy:grid.ny
     field.v(i,j) = field.v(i,j)-param.dt*dP_dy(i,j);
    end
end
for i=w.idx:w.idx+w.ldx-1
    for j=2:w.idy-2
       field.v(i,j) = field.v(i,j)-param.dt*dP_dy(i,j);
    end
end
for i=w.idx+w.ldx:grid.nx+1
    for j=2:grid.ny
       field.v(i,j) = field.v(i,j)-param.dt*dP_dy(i,j);
    end
end

end

function[] = visualize(field,grid)
[x,y] = meshgrid(1:1:grid.nx,1:1:grid.ny);
u = zeros(grid.nx,grid.ny);
v = zeros(grid.nx,grid.ny);

for i = 2:grid.nx+1
    for j = 2:grid.ny+1
        u(i-1,j-1) = (field.u(i,j)+field.u(i-1,j))/2;
        v(i-1,j-1) = (field.v(i,j)+field.v(i,j-1))/2;    
    end
end
quiver(x.',y.',u,v)
end

% figure(2)
% surface(v.')

















