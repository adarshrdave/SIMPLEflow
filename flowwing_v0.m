%% F18 CFD Project - incompressible flow around wing
%% Setup
%                     N
% (1,ny) ----------------------------- (nx,ny)
%      |                             |
%      |    (w.idx,w.idy)            |
%   W  |          |-------|          |   E
%      |          |-------|          |
%      |                             |
%      |                             |
%(1,1) ----------------------------- (nx,1)
%                     S

% grid parameters
grid.Lx = 5;    %x-length of box
grid.Ly = 2;    %y-length of box
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
%when using values from other fields in the formulas, each cell
%corresponds to:  P(i,j)<->U(i,j)<->V(i-1,j)  - this is what's slightly
%annoying about using the big arrays for iterating. Could be fixed by
%adding one useless row to V


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
%set all corners (P,U,V) to 0, am too lazy to type this out atm

%West - Neumann
field.u(1,2:grid.ny+1) = param.vIN*sin(param.alpha);
field.v(1,2:grid.ny) = 2*param.vIN*cos(param.alpha) - field.v(2,2:grid.ny);
field.p(1,2:grid.ny+1) = field.p(2,2:grid.ny+1);

%South - Neumann
field.u(2:grid.nx,grid.ny+2) = 2*param.vIN*cos(param.alpha)-field.u(2:grid.nx,grid.ny+1);
field.v(2:grid.nx+1,grid.ny+1) = param.vIN*sin(param.alpha);
field.p(2:grid.nx+1,grid.ny+2) = field.p(2:grid.nx+1,grid.ny+1);

%North - far-field
field.v(2:grid.nx+1,1) = 0;
field.u(2:grid.nx,1) = field.u(2:grid.nx,2);
field.p(2:grid.nx+1,1) = field.p(2:grid.nx+1,2);

%East - far-field
field.v(grid.nx+2,2:grid.ny) = field.v(grid.nx+1,2:grid.ny);
field.u(grid.nx+1,2:grid.ny+1) = 0;
field.p(grid.nx+2,2:grid.ny+1) = field.p(grid.nx+1,2:grid.ny+1);


%% Wing - 
%this is if v(i,j) corresponds to P(i,j) - change v(i,j) to v(i-1,j) 
    %desired size in 'm'
w.Ly = .5;
w.Lx = 1;
    %size in no cells
w.ldy = floor(w.Ly/grid.dy);
w.ldx = floor(w.Lx/grid.dx);
    %position in  grid
w.idy = ceil((grid.ny/2)-(w.ldy/2));
w.idx = ceil((grid.nx/2)-(w.ldx/2));

%BC for wing - 
%top
field.v(w.idx:w.idx+w.ldx-1,w.idy)=0;
field.u(w.idx:w.idx+w.ldx-2,w.idy)= ...
    -field.u(w.idx:w.idx+w.ldx-2,w.idy-1);
field.p(w.idx:w.idx+w.ldx-1,w.idy)=field.p(w.idx:w.idx+w.ldx-1,w.idy-1);
%bot
field.v(w.idx:w.idx+w.ldx-1,w.idy+w.ldy)=0;
field.u(w.idx:w.idx+w.ldx-2,w.idy+w.ldy-1)= ...
    -field.u(w.idx:w.idx+w.ldx-2,w.idy+w.ldy-1);
field.p(w.idx:w.idx+w.ldx-1,w.idy+w.ldy-1)=field.p(w.idx:w.idx+w.ldx-1,w.idy+w.ldy);

%left - only v,u for now
field.u(w.idx-1,w.idy:w.idy+w.ldy-1)=0;
field.v(w.idx,w.idy+1:w.idy+w.ldy-1)=...
    -field.v(w.idx-1,w.idy+1:w.idy+w.ldy-1);
field.p(w.idx,w.idy:w.idy+w.ldy-1) = field.p(w.idx-1,w.idy:w.idy+w.ldy-1);

%right, only v,u for now
field.u(w.idx+w.ldx-1,w.idy:w.idy+w.ldy-1)=0;
field.v(w.idx+w.ldx-1,w.idy+1:w.idy+w.ldy-1)=...
    -field.v(w.idx+w.ldx,w.idy+1:w.idy+w.ldy-1);
field.p(w.idx+w.ldx-1,w.idy:w.idy+w.ldy-1) = ...
    field.p(w.idx+w.ldx,w.idy:w.idy+w.ldy-1);

%% Parameters (arbitraryly chosen for now)

% Initialization, Parameters
param.Re = 100;
param.rho = 1.184;
param.mu = 1;
param.T = 5;
param.dt = .1;
param.tsteps = param.T/param.dt;
param.eps = .01;
param.Q=0.5;
param.s=10; %maximum speed in m/s

% etc

%% FSM
%P: for i=2:nx+1
%       for j=2:ny+1
%           ...

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
































