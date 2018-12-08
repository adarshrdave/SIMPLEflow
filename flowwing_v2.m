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

% grid parameters
grid.Lx = 10;    %x-length of box
grid.Ly = 5;    %y-length of box
grid.dx = .1;   %cell dimensions
grid.dy = grid.dx;

%number of int. cells: nx*ny - ie. interior p nodes
grid.nx = grid.Lx/grid.dx;
grid.ny = grid.Ly/grid.dx;

%set up arrays for p,u,v - big arrays including ghost nodes
field.p = zeros(grid.nx+2,grid.ny+2);
field.u = zeros(grid.nx+1,grid.ny+2);
field.v = zeros(grid.nx+2,grid.ny+1);

%information field: this is every value we need plus the corner values in
%each field, which are not used. If you take the fields without ghost
%nodes (but with values exactly on boundary), the dimensions are as in the lecture notes:
%P: ny*nx       U: ny*(nx+1)        V: (ny+1)*nx
%this will be the array we visualize

%if you look at only the interior nodes, you get:
%P: ny*nx       U: ny*nx-1         V: (ny-1)*nx
%this will be the arrays we iterate over after BCs have been set.
%when iterating over the complete interior (no sections yet), the loops
%over the big arrays would be (remember i is for y (rows) and j for x (cols):
%P: for i=2:ny+1
%       for j=2:nx+1
%           ...

%U: for i=2:ny+1
%        for j=2:nx
%           ...
%
%V: for i=2:ny
%       for j=2:ny+1
%           ...

% Inflow across all of W & S
param.alpha = pi/4;
param.vIN = 0.5; % magnitude of velocity of inflow

% Wing - 
    %desired size in 'm'
w.Ly = 1;
w.Lx = 2;
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
param.T = 1;
param.dt = .01;
param.tsteps = param.T/param.dt;
param.eps = .01;
param.Q=0.5;
param.s=10; %maximum speed in m/s

% etc

[A1,A2,A3,A4,ni] = init_pMatrices(grid,w);
M.A1 = A1;
M.A2 = A2;
M.A3 = A3;
M.A4 = A4;

%% FSM

%P: for i=2:nx+1
%       for j=2:ny+1
%           ...
for t=1:param.tsteps
    
field = setBC(field,grid,param,w);
figure(1)
surface(field.p.')
%pause(1.5)

%Intermediate velocity field
field = computeVF_sections(field,grid,param,w);

field = setBC(field,grid,param,w);

field = solve_P(field,grid,param,w,ni,M);

field.v = field.v-param.dt*dPdy(field,grid);
field.u = field.u-param.dt*dPdx(field,grid);


end

%% Loops with wing + test space
for dummy=1:1
%looping over interior nodes only
%for each timestep:

%field=setBC(field,grid,param);

%Intermediate velocity field
U=field.u;
V=field.v;

%U loops- for everything computed at u(i,j)
%Part 1 - left
for i=2:w.idx-2
    for j=2:grid.ny+1
        %calculations
        %field.u(i,j)=..
    end
end
%part 2 - top
for i=w.idx-1:w.idx+w.ldx-1
    for j=w.idy+w.ldy:grid.ny+1
        %calculations
        %field.u(i,j)=..
    end
end
%part 3 - bot
for i=w.idx-1:w.idx+w.ldx-1
    for j=2:w.idy-1
     %calculations
     %field.u(i,j)=..
    end
end
%part 4 - right
for i=w.idx+w.ldx:grid.nx
    for j=2:grid.ny+1
       %calculations
       %field.u(i,j)=..
    end
end

%V loops - for everything computed at V(i,j)
for i=2:w.idx-1
    for j=2:grid.ny
        %calculations
        %field.v(i,j)=..
    end
end
for i=w.idx:w.idx+w.ldx-1
    for j=w.idy+w.ldy:grid.ny
     %calculations
      %field.v(i,j)=..
    end
end
for i=w.idx:w.idx+w.ldx-1
    for j=2:w.idy-2
       %calculations
        %field.v(i,j)=..
    end
end
for i=w.idx+w.ldx:grid.nx+1
    for j=2:grid.ny
       %calculations
        %field.v(i,j)=..
    end
end

%P loops - for everything computed at P(i,j)
for i=2:w.idx-1
    for j=2:grid.ny+1
         %field.p(i,j)=..
    end
end
for i=w.idx:w.idx+w.ldx-1
    for j=w.idy+w.ldy:grid.ny+1
          %field.p(i,j)=..
    end
end
for i=w.idx:w.idx+w.ldx-1
  for j=2:w.idy-1
        %field.p(i,j)=..
  end
end
for i=w.idx+w.ldx:grid.nx+1
    for j=2:grid.ny+1
          %field.p(i,j)=..
    end
end
[b1,b2,b3,b4]  = compute_RHS(field,grid,param,w);
A=rand(100,100);
b=rand(100,1);
% x=pcg(A,b);
% x1 = A\b;
x=GS(M.A1,b1);
%x1=M.A1\b1;
norm(M.A1*x1-b1)
%field = setBC_nowing(field,grid,param);
%field = setBC(field,grid,param,w);
end
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

%at P(i,j)
function[du_dx] =  dudx(field,grid)
du_dx = zeros(grid.nx,grid.ny);
for i=2:grid.nx+1
    for j=2:grid.ny+1
        du_dx(i-1,j-1) = (field.u(i,j)-field.u(i-1,j))/grid.dx;
    end
end
end

%at P(i,j)
function[dv_dy] =  dvdy(field,grid)
dv_dy = zeros(grid.nx,grid.ny);
for i=2:grid.nx+1
    for j=2:grid.ny+1
        dv_dy(i-1,j-1) = (field.v(i,j)-field.v(i,j-1))/grid.dy;
    end
end
end

%computes matrices for poisson solvers for each section 
%only needs to be done once before starting timesteps
%needs b from compute_RHS()!
function[A1,A2,A3,A4,ni] = init_pMatrices(grid,w)
ni = zeros(4,2);
%--------------------PART 1-------------------------
nx1 = w.idx-2;    ny1 = grid.ny;
ni(1,1) = nx1;  ni(1,2) = ny1;
A1 = zeros(nx1*ny1,nx1*ny1);
for i=2:nx1-1
    %south
    A1((i-1)*ny1+1,(i-1)*ny1+1)=-3;
    A1((i-1)*ny1+1,(i-1)*ny1+2)=1;
    A1((i-1)*ny1+1,(i-1)*ny1+1+ny1)=1;
    A1((i-1)*ny1+1,(i-1)*ny1+1-ny1)=1;
    
    %north
    A1((i)*ny1,(i)*ny1)=-3;
    A1((i)*ny1,(i)*ny1-1)=1;
    A1((i)*ny1,(i)*ny1+ny1)=1;
    A1((i)*ny1,(i)*ny1-ny1)=1;
     for j=2:ny1-1
        %interior
        k = (i-1)*ny1+j;
        A1(k,k) = -4;
        A1(k,k-1) = 1;
        A1(k,k+1) = 1;
        A1(k,k+ny1) = 1;
        A1(k,k-ny1) = 1;
    end   
    
end 

for j=2:ny1-1
    %west
    A1(j,j) = -3;
    A1(j,j+1) = 1;
    A1(j,j-1) = 1;
    A1(j,j+ny1) = 1;
    
    %east
    if((j<w.idy-1)||(j>=w.idy+w.ldy-1))
    A1((nx1-1)*ny1+j,(nx1-1)*ny1+j) = -4;
    A1((nx1-1)*ny1+j,(nx1-1)*ny1+j-1) = 1;
    A1((nx1-1)*ny1+j,(nx1-1)*ny1+j+1) = 1; 
    A1((nx1-1)*ny1+j,(nx1-2)*ny1+j) = 1;
    else
    A1((nx1-1)*ny1+j,(nx1-1)*ny1+j) = -3;
    A1((nx1-1)*ny1+j,(nx1-1)*ny1+j-1) = 1;
    A1((nx1-1)*ny1+j,(nx1-1)*ny1+j+1) = 1; 
    A1((nx1-1)*ny1+j,(nx1-2)*ny1+j) = 1;
    end
end
%corners
%SW
A1(1,1) = -2;
A1(1,2) = 1;
A1(1,1+ny1) = 1;
%NW
A1(ny1,ny1) = -2;
A1(ny1,ny1-1) = 1;
A1(ny1,ny1*2) = 1;
%SE
A1((nx1-1)*ny1+1,(nx1-1)*ny1+1) = -3;
A1((nx1-1)*ny1+1,(nx1-1)*ny1+2) = 1;
A1((nx1-1)*ny1+1,(nx1-2)*ny1+1) = 1;
%NE
A1(nx1*ny1,nx1*ny1) = -3;
A1(nx1*ny1,nx1*ny1-1) = 1;
A1(nx1*ny1,nx1*ny1-ny1) = 1;

A1=(1/(grid.dx^2))*A1;
%-------------------PART 2---------------------------
nx2 = w.ldx;    ny2 = grid.ny+1-(w.idy+w.ldy)+1;
ni(2,1) = nx2;  ni(2,2) = ny2;
A2 = zeros(nx2*ny2,nx2*ny2);
for i=2:nx2-1
    %south
    A2((i-1)*ny2+1,(i-1)*ny2+1)=-3;
    A2((i-1)*ny2+1,(i-1)*ny2+2)=1;
    A2((i-1)*ny2+1,(i-1)*ny2+1+ny2)=1;
    A2((i-1)*ny2+1,(i-1)*ny2+1-ny2)=1;
    
    %north
    A2((i)*ny2,(i)*ny2)=-3;
    A2((i)*ny2,(i)*ny2-1)=1;
    A2((i)*ny2,(i)*ny2+ny2)=1;
    A2((i)*ny2,(i)*ny2-ny2)=1;
     for j=2:ny2-1
        %interior
        k = (i-1)*ny2+j;
        A2(k,k) = -4;
        A2(k,k-1) = 1;
        A2(k,k+1) = 1;
        A2(k,k+ny2) = 1;
        A2(k,k-ny2) = 1;
    end   
    
end 

for j=2:ny2-1
    %west
    A2(j,j) = -4;
    A2(j,j+1) = 1;
    A2(j,j-1) = 1;
    A2(j,j+ny2) = 1;
    
    %east
    A2((nx2-1)*ny2+j,(nx2-1)*ny2+j) = -4;
    A2((nx2-1)*ny2+j,(nx2-1)*ny2+j-1) = 1;
    A2((nx2-1)*ny2+j,(nx2-1)*ny2+j+1) = 1; 
    A2((nx2-1)*ny2+j,(nx2-2)*ny2+j) = 1;
end

%corners
%SW
A2(1,1) = -3;
A2(1,2) = 1;
A2(1,1+ny2) = 1;
%NW
A2(ny2,ny2) = -3;
A2(ny2,ny2-1) = 1;
A2(ny2,ny2*2) = 1;
%SE
A2((nx2-1)*ny2+1,(nx2-1)*ny2+1) = -3;
A2((nx2-1)*ny2+1,(nx2-1)*ny2+2) = 1;
A2((nx2-1)*ny2+1,(nx2-2)*ny2+1) = 1;
%NE
A2(nx2*ny2,nx2*ny2) = -3;
A2(nx2*ny2,nx2*ny2-1) = 1;
A2(nx2*ny2,nx2*ny2-ny2) = 1;

A2=(1/(grid.dx^2))*A2;

%--------------------PART 3----------------------
nx3 = w.ldx;    ny3 = w.idy-2;
ni(3,1) = nx3;  ni(3,2) = ny3;
A3 = zeros(nx3*ny3,nx3*ny3);
for i=2:nx3-1
    %south
    A3((i-1)*ny3+1,(i-1)*ny3+1)=-3;
    A3((i-1)*ny3+1,(i-1)*ny3+2)=1;
    A3((i-1)*ny3+1,(i-1)*ny3+1+ny3)=1;
    A3((i-1)*ny3+1,(i-1)*ny3+1-ny3)=1;
    
    %north
    A3((i)*ny3,(i)*ny3)=-3;
    A3((i)*ny3,(i)*ny3-1)=1;
    A3((i)*ny3,(i)*ny3+ny3)=1;
    A3((i)*ny3,(i)*ny3-ny3)=1;
     for j=2:ny3-1
        %interior
        k = (i-1)*ny3+j;
        A3(k,k) = -4;
        A3(k,k-1) = 1;
        A3(k,k+1) = 1;
        A3(k,k+ny3) = 1;
        A3(k,k-ny3) = 1;
    end   
    
end 

for j=2:ny3-1
    %west
    A3(j,j) = -4;
    A3(j,j+1) = 1;
    A3(j,j-1) = 1;
    A3(j,j+ny3) = 1;
    
    %east
    A3((nx3-1)*ny3+j,(nx3-1)*ny3+j) = -4;
    A3((nx3-1)*ny3+j,(nx3-1)*ny3+j-1) = 1;
    A3((nx3-1)*ny3+j,(nx3-1)*ny3+j+1) = 1; 
    A3((nx3-1)*ny3+j,(nx3-2)*ny3+j) = 1;
end

%corners
%SW
A3(1,1) = -3;
A3(1,2) = 1;
A3(1,1+ny3) = 1;
%NW
A3(ny3,ny3) = -3;
A3(ny3,ny3-1) = 1;
A3(ny3,ny3*2) = 1;
%SE
A3((nx3-1)*ny3+1,(nx3-1)*ny3+1) = -3;
A3((nx3-1)*ny3+1,(nx3-1)*ny3+2) = 1;
A3((nx3-1)*ny3+1,(nx3-2)*ny3+1) = 1;
%NE
A3(nx3*ny3,nx3*ny3) = -3;
A3(nx3*ny3,nx3*ny3-1) = 1;
A3(nx3*ny3,nx3*ny3-ny3) = 1;

A3=(1/(grid.dx^2))*A3;

%-----------------PART 4-----------------
nx4 = grid.nx+2-(w.idx+w.ldx);   ny4 = grid.ny; 
ni(4,1) = nx4;  ni(4,2) = ny4;
A4 = zeros(nx4*ny4,nx4*ny4);

for i=2:nx4-1
    %south
    A4((i-1)*ny4+1,(i-1)*ny4+1)=-3;
    A4((i-1)*ny4+1,(i-1)*ny4+2)=1;
    A4((i-1)*ny4+1,(i-1)*ny4+1+ny4)=1;
    A4((i-1)*ny4+1,(i-1)*ny4+1-ny4)=1;
    
    %north
    A4((i)*ny4,(i)*ny4)=-3;
    A4((i)*ny4,(i)*ny4-1)=1;
    A4((i)*ny4,(i)*ny4+ny4)=1;
    A4((i)*ny4,(i)*ny4-ny4)=1;
    
     for j=2:ny4-1
        %interior
        k = (i-1)*ny4+j;
        A4(k,k) = -4;
        A4(k,k-1) = 1;
        A4(k,k+1) = 1;
        A4(k,k+ny4) = 1;
        A4(k,k-ny4) = 1;
    end   
    
end 

for j=2:ny4-1
    %west
    if((j<w.idy-1)||(j>=w.idy+w.ldy-1))
    A4(j,j) = -4;
    A4(j,j+1) = 1;
    A4(j,j-1) = 1;
    A4(j,j+ny4) = 1;  
    else
    A4(j,j) = -3;
    A4(j,j+1) = 1;
    A4(j,j-1) = 1;
    A4(j,j+ny4) = 1;
    end
    %east
    A4((nx4-1)*ny4+j,(nx4-1)*ny4+j) = -3;
    A4((nx4-1)*ny4+j,(nx4-1)*ny4+j-1) = 1;
    A4((nx4-1)*ny4+j,(nx4-1)*ny4+j+1) = 1; 
    A4((nx4-1)*ny4+j,(nx4-2)*ny4+j) = 1;
    
end
%corners
%SW
A4(1,1) = -3;
A4(1,2) = 1;
A4(1,1+ny4) = 1;
%NW
A4(ny4,ny4) = -3;
A4(ny4,ny4-1) = 1;
A4(ny4,ny4*2) = 1;
%SE
A4((nx4-1)*ny4+1,(nx4-1)*ny4+1) = 2;
A4((nx4-1)*ny4+1,(nx4-1)*ny4+2) = 1;
A4((nx4-1)*ny4+1,(nx4-2)*ny4+1) = 1;
%NE
A4(nx4*ny4,nx4*ny4) = -2;
A4(nx4*ny4,nx4*ny4-1) = 1;
A4(nx4*ny4,nx4*ny4-ny4) = 1;

A4=(1/(grid.dx^2))*A4;


end

function[field] = setBC(field,grid,param,w)
field.p(2,2)=0;
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
field.v(2:grid.nx+1,grid.ny+1) = 0;
field.u(2:grid.nx,grid.ny+2) = field.u(2:grid.nx,grid.ny+1);
field.p(2:grid.nx+1,grid.ny+2) = field.p(2:grid.nx+1,grid.ny+1);

%East - far-field
field.v(grid.nx+2,2:grid.ny) = field.v(grid.nx+1,2:grid.ny);
field.u(grid.nx+1,2:grid.ny+1) = 0;
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
%computes b= (1/dt)*(du/dx+dy/dy)+g for each section
%g being the BCs as within section boundaries (P from last ts)
%this is called every iteration in solver
function[b1,b2,b3,b4] = compute_RHS(field,grid,param,w)

%-----------PART 1---------------
nx1 = w.idx-2;    ny1 = grid.ny; %can be set to ni(i,j) later 
b1 = zeros(nx1,ny1);
for i=2:w.idx-1
    for j=2:grid.ny+1
         b1(i-1,j-1) = (1/param.dt)*((field.u(i,j)-field.u(i-1,j))/grid.dx +...
             (field.v(i,j)-field.v(i,j-1))/grid.dy);    
         if (i==w.idx-1&&(j<w.idy||j>=w.idy+w.ldy)) %open boundaries
         b1(i-1,j-1) = b1(i-1,j-1)-(1/(grid.dx^2))*field.p(i+1,j);
         end
    end
end
b1 = reshape(b1.',[(nx1)*ny1 1]); 

%-------------Part 2-------------
nx2 = w.ldx;    ny2 = grid.ny+1-(w.idy+w.ldy)+1;
b2 = zeros(nx2,ny2);
for i=w.idx:w.idx+w.ldx-1
    for j=w.idy+w.ldy:grid.ny+1
         b2(i-w.idx+1,j-(w.idy+w.ldy)+1) = (1/param.dt)*((field.u(i,j)-field.u(i-1,j))/grid.dx +...
             (field.v(i,j)-field.v(i,j-1))/grid.dy);    
         if (i==w.idx) %open boundaries
         b2(i-w.idx+1,j-(w.idy+w.ldy)+1) = ...
             b2(i-w.idx+1,j-(w.idy+w.ldy)+1)-(1/(grid.dx^2))*field.p(i-1,j);
         end
         if (i==w.idx+w.ldx-1)
            b2(i-w.idx+1,j-(w.idy+w.ldy)+1) = ...
             b2(i-w.idx+1,j-(w.idy+w.ldy)+1)-(1/(grid.dx^2))*field.p(i+1,j);  
         end
    end
end
b2 = reshape(b2.',[(nx2)*ny2 1]); 

%-------------Part 3-----------------
nx3 = w.ldx;    ny3 = w.idy-2;
b3 = zeros(nx3,ny3);
for i=w.idx:w.idx+w.ldx-1
    for j=2:w.idy-1
         b3(i-w.idx+1,j-1) = (1/param.dt)*((field.u(i,j)-field.u(i-1,j))/grid.dx +...
             (field.v(i,j)-field.v(i,j-1))/grid.dy);    
         if (i==w.idx) %open boundaries
         b3(i-w.idx+1,j-1) = ...
             b3(i-w.idx+1,j-1)-(1/(grid.dx^2))*field.p(i-1,j);
         end
         if (i==w.idx+w.ldx-1)
            b3(i-w.idx+1,j-1) = ...
             b3(i-w.idx+1,j-1)-(1/(grid.dx^2))*field.p(i+1,j);  
         end
    end
end
b3 = reshape(b3.',[nx3*ny3 1]); 

%-----------------PART 4-------------------
nx4 = grid.nx+2-(w.idx+w.ldx);   ny4 = grid.ny;  %- dimension of part 1
b4 = zeros(nx4,ny4);

for i=w.idx+w.ldx:grid.nx+1
    for j=2:grid.ny+1
         b4(i-(w.idx+w.ldx)+1,j-1) = (1/param.dt)*((field.u(i,j)-field.u(i-1,j))/grid.dx +...
             (field.v(i,j)-field.v(i,j-1))/grid.dy);    
         if ((i==w.idx+w.ldx)&&(j<w.idy||j>=w.idy+w.ldy)) %open boundaries
         b4(i-(w.idx+w.ldx)+1,j-1) = ...
             b4(i-(w.idx+w.ldx)+1,j-1)-(1/(grid.dx^2))*field.p(i+1,j);
         end
    end
end
b4 = reshape(b4.',[nx4*ny4 1]); 


end

function[field] = solve_P(field,grid,param,w,ni,M)

[b1,b2,b3,b4]  = compute_RHS(field,grid,param,w);
p1 = M.A1\b1;%GS(A1,b1);
p2 = M.A2\b2;%GS(A2,b2);
p3 = M.A3\b3;%GS(A3,b3);
p4 = M.A4\b4;%GS(A4,b4);

p1 = reshape(p1,[ni(1,2) ni(1,1)]).';
p2 = reshape(p2,[ni(2,2) ni(2,1)]).';
p3 = reshape(p3,[ni(3,2) ni(3,1)]).';
p4 = reshape(p4,[ni(4,2) ni(4,1)]).';

%put solutions for sections back into field

field.p(2:w.idx-1,2:grid.ny+1) = p1;
field.p(w.idx:w.idx+w.ldx-1,w.idy+w.ldy:grid.ny+1) = p2;
field.p(w.idx:w.idx+w.ldx-1,2:w.idy-1) = p3;
field.p(w.idx+w.ldx:grid.nx+1,2:grid.ny+1) =p4;

end

function[x] = GS(A,b)
n=size(b,1);
normVal=Inf; 
tol=1e-1; itr=0;
x= zeros(size(b));
while normVal>tol
    x_old=x;
    for i=1:n
        sigma=0;
        for j=1:i-1
                sigma=sigma+A(i,j)*x(j);
        end
        for j=i+1:n
                sigma=sigma+A(i,j)*x_old(j);
        end
        x(i)=(1/A(i,i))*(b(i)-sigma);
    end
    itr=itr+1;
    normVal=norm(x_old-x);
end
end

function[field] = computeVF(field,grid,param,w)
U=field.u;
V=field.v;
for i=2:grid.nx
    for j=2:grid.ny+1
        %Centered-Difference
        uuxCD=((U(i,j)+U(i+1,j))^2-(U(i,j)+U(i-1,j))^2)/(4*grid.dx);
        uvyCD=((V(i,j)+V(i+1,j))*(U(i,j)+U(i,j+1))-(V(i,j-1)+V(i+1,j-1))*(U(i,j-1)+U(i,j)))/(4*grid.dy);
        %Upwind
        uuxUW=(abs(U(i,j)+U(i+1,j))*(U(i,j)-U(i+1,j))-abs(U(i-1,j)+U(i,j))*(U(i-1,j)-U(i,j)))*param.Q/(4*grid.dx);
        uvyUW=(abs(V(i,j)+V(i+1,j))*(U(i,j)-U(i,j+1))-abs(V(i,j-1)+V(i+1,j-1))*(U(i,j-1)-U(i,j)))*param.Q/(4*grid.dx);     
        %Diffusion
        uxx=(U(i+1,j)-2*U(i,j)+U(i-1,j))/(grid.dx*grid.dx);
        uyy=(U(i,j+1)-2*U(i,j)+U(i,j-1))/(grid.dy*grid.dy);
        
        %Update
        field.u(i,j)=U(i,j)+param.dt*((1/param.Re)*(uxx+uyy)-(uuxCD+uuxUW+uvyCD+uvyUW));
    end
end

%V: for i=2:nx+1
%       for j=2:ny
%           ...
for i=2:grid.nx+1
    for j=2:grid.ny
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
        field.v(i,j)=V(i,j)+param.dt*((1/param.Re)*(vxx+vyy)-(uvxCD+uvxUW+vvyCD+vvyUW));
    end
end

end
function[field] = setBC_nowing(field,grid,param)
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
field.v(2:grid.nx+1,grid.ny+1) = 0;
field.u(2:grid.nx,grid.ny+2) = field.u(2:grid.nx,grid.ny+1);
field.p(2:grid.nx+1,grid.ny+2) = field.p(2:grid.nx+1,grid.ny+1);

%East - far-field
field.v(grid.nx+2,2:grid.ny) = field.v(grid.nx+1,2:grid.ny);
field.u(grid.nx+1,2:grid.ny+1) = 0;
field.p(grid.nx+2,2:grid.ny+1) = field.p(grid.nx+1,2:grid.ny+1);

end

%solver without wing 
function[p_new] = pressure_poisson(field,param,grid,w)
%first right hand side b: du/dx + dv/dy at P(i,j)
b = zeros(grid.nx,grid.ny);
b = dvdy(field,grid) +dudx(field,grid);
b = b/param.dt;
%transform to vector for solving
b = reshape(b.',[grid.nx*grid.ny 1]);

%set up matrix for solving A*p = b - this should be removed from function
%later as it does not change - for partitions there will be four different
%ones
Ap = zeros(grid.nx*grid.ny,grid.nx*grid.ny);

for i=2:grid.nx-1
    %south
    Ap((i-1)*grid.ny+1,(i-1)*grid.ny+1)=-3;
    Ap((i-1)*grid.ny+1,(i-1)*grid.ny+2)=1;
    Ap((i-1)*grid.ny+1,(i-1)*grid.ny+1+grid.ny)=1;
    Ap((i-1)*grid.ny+1,(i-1)*grid.ny+1-grid.ny)=1;
    
    %north
    Ap((i)*grid.ny,(i)*grid.ny)=-3;
    Ap((i)*grid.ny,(i)*grid.ny-1)=1;
    Ap((i)*grid.ny,(i)*grid.ny+grid.ny)=1;
    Ap((i)*grid.ny,(i)*grid.ny-grid.ny)=1;
    
    for j=2:grid.ny-1
        %interior
        k = (i-1)*grid.ny+j;
        Ap(k,k) = -4;
        Ap(k,k-1) = 1;
        Ap(k,k+1) = 1;
        Ap(k,k+grid.ny) = 1;
        Ap(k,k-grid.ny) = 1;
    end   
end

for j=2:grid.ny-1
    %west
    Ap(j,j) = -3;
    Ap(j,j+1) = 1;
    Ap(j,j-1) = 1;
    Ap(j,j+grid.ny) = 1;
    
    %east
    Ap((grid.nx-1)*grid.ny+j,(grid.nx-1)*grid.ny+j) = 3;
    Ap((grid.nx-1)*grid.ny+j,(grid.nx-1)*grid.ny+j-1) = 1;
    Ap((grid.nx-1)*grid.ny+j,(grid.nx-1)*grid.ny+j+1) = 1; 
    Ap((grid.nx-1)*grid.ny+j,(grid.nx-2)*grid.ny+j) = 1;
end

%corners
%SW
Ap(1,1) = -2;
Ap(1,2) = 1;
Ap(1,1+grid.ny) = 1;
%NW
Ap(grid.ny,grid.ny) = -2;
Ap(grid.ny,grid.ny-1) = 1;
Ap(grid.ny,grid.ny*2) = 1;
%SE
Ap((grid.nx-1)*grid.ny+1,(grid.nx-1)*grid.ny+1) = -2;
Ap((grid.nx-1)*grid.ny+1,(grid.nx-1)*grid.ny+2) = 1;
Ap((grid.nx-1)*grid.ny+1,(grid.nx-2)*grid.ny+1) = 1;
%NE
Ap(grid.nx*grid.ny,grid.nx*grid.ny) = -2;
Ap(grid.nx*grid.ny,grid.nx*grid.ny-1) = 1;
Ap(grid.nx*grid.ny,grid.nx*grid.ny-grid.ny) = 1;

Ap =(1/grid.dx^2)*Ap;
p_new = Ap/b.';

p_new = reshape(p_new,grid.ny, grid.nx).';
end






















