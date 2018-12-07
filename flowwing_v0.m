%% F18 CFD Project - incompressible flow around wing
%% Setup
%                     N
% (1,ny) ----------------------------- (nx,ny)
%      |                             |
%      |                             |
%   W  |          |-------|          |   E
%      |          |-------|          |
%      |   (w.idx,w.idy)             |
%      |                             |
%(1,1) ----------------------------- (nx,1)
%                     S

% grid parameters
grid.Lx = 20;    %x-length of box
grid.Ly = 10;    %y-length of box
grid.dx = .5;   %cell dimensions
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
param.alpha = pi/3;
param.vIN = 5; % magnitude of velocity of inflow

% BCs - far-field for N & E, Neumann for W & S,
%       Dirchilet at wing
%I am now setting the ghost nodes and boundary values directly in the big
%array, I think this is better than having separate vectors. 
%all this code will be put in a setBC func, together with the wing BC,
%which is called at every timestep before calculating

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
field.p(2:grid.nx+1,grid.ny+2) = field.p(2:grid.nx+1,grid.ny+2);

%East - far-field
field.v(grid.nx+2,2:grid.ny) = field.v(grid.nx+1,2:grid.ny);
field.u(grid.nx+1,2:grid.ny+1) = 0;
field.p(grid.nx+2,2:grid.ny+1) = field.p(grid.nx+1,2:grid.ny+1);


%% Wing - 
    %desired size in 'm'
w.Ly = .5;
w.Lx = 1;
    %size in no cells
w.ldy = floor(w.Ly/grid.dy);
w.ldx = floor(w.Lx/grid.dx);
    %position in  grid - left-bottom corner-cell of wing
w.idy = ceil((grid.ny/2)-(w.ldy/2));
w.idx = ceil((grid.nx/2)-(w.ldx/2));

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

%% Parameters (arbitraryly chosen for now)

% Initialization, Parameters
param.Re = 100;
param.rho = 1.184;
param.mu = 1;
param.T = 5;
param.dt = .01;
param.tsteps = param.T/param.dt;
param.eps = .01;
param.Q=0.5;
param.s=10; %maximum speed in m/s

% etc


%% FSM
%P: for i=2:nx+1
%       for j=2:ny+1
%           ...

for t=1:param.tsteps
    
field=setBC_nowing(field,grid,param);

%Intermediate velocity field
U=field.u;
V=field.v;

%Update U
%U: for i=2:nx
%        for j=2:ny+1
%           ...
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
        field.u(i,j)=U(i,j)+param.dt*(1/param.Re*(uxx+uyy)-(uuxCD+uuxUW+uvyCD+uvyUW));
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
        field.v(i,j)=V(i,j)+param.dt*(1/param.Re*(vxx+vyy)-(uvxCD+uvxUW+vvyCD+vvyUW));
    end
end

p_new = pressure_poisson(field,param,grid,w);
field.p(2:grid.nx+1,2:grid.ny+1) = p_new;

field.v = field.v-param.dt*dPdy(field,grid);
field.u = field.u-param.dt*dPdx(field,grid);

figure(1)
surface(field.u.')
end
% for i=1:param.tsteps
%     
%     %set boundary conditions BC+ Wing loop over sections of interior nodes
%     %VF = Vn + param.dt*(A+B) - for u, v separately
%     
%     %compute P_(n+1) from grad(grad(Pp))=(1/param.dt)*(grad(VF)) - at P
%     %nodes
%    
%     %V = VF-param.dt*grad(P_(n+1));
% end
%% Loops with wing
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
du_dx = zeros(grid.ny,grid.ny);
for i=2:grid.nx+1
    for j=2:grid.ny+1
        du_dx(i-1,j-1) = (field.u(i,j)-field.u(i-1,j))/grid.dx;
    end
end
end

%at P(i,j)
function[dv_dy] =  dvdy(field,grid)
dv_dy = zeros(grid.ny,grid.ny);
for i=2:grid.nx+1
    for j=2:grid.ny+1
        dv_dy(i-1,j-1) = (field.v(i,j)-field.v(i,j-1))/grid.dy;
    end
end
end


%new solver with sections - just starting this here can merge later
function[p_new] = pressure_solver(field,param,w)

%---------part 1: left-------------
nx1 = w.idx-2;    ny1 = grid.ny; %- dimension of part 1
b1 = zeros(nx1,ny1);
A1 = zeros((nx1)*ny1,(nx1)*ny1);

%first right hand side b: du/dx + dv/dy + P- at open boundaries
for i=2:w.idx-1
    for j=2:grid.ny+1
         b1(i-1,j-1) = (field.u(i,j)-field.u(i-1,j))/grid.dx +...
             (field.v(i,j)-field.v(i,j-1))/grid.dy;    
         if (i==w.idx-1&&(j<w.idy||j>=w.idy+w.ldy)) %open boundaries
         b1(i-1,j-1) = b1(i-1,j-1)-(1/(grid.dx^2))*field.p(i+1,j);
         end
    end
end
b1=b1/param.dt;
b1 = reshape(b1.',[(nx1)*ny1 1]); 

%now matrix for solver

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
        A1(k,k+grid.ny) = 1;
        A1(k,k-grid.ny) = 1;
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
    A1((nx1-1)*ny1+j,(nx1-1)*ny1+j) = 3;
    A1((nx1-1)*ny1+j,(nx1-1)*ny1+j-1) = 1;
    A1((nx1-1)*ny1+j,(nx1-1)*ny1+j+1) = 1; 
    A1((nx1-1)*ny1+j,(nx1-2)*ny1+j) = 1;
    end
end

%p1 = GaussSeidel(A1,b1);

%------------part 2: top----------------



end


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
























