%% F18 CFD Project - incompressible flow around wing
%% Setup
%                     N
% (1,1) ----------------------------- (1,nx)
%      |                             |
%      |    (w.idy,w.idx)            |
%   W  |          |-------|          |   E
%      |          |-------|          |
%      |                             |
%      |                             |
%(ny,1) ----------------------------- (ny,nx)
%                     S

% grid parameters
grid.Lx = 2;    %x-length of box
grid.Ly = 1;    %y-length of box
grid.dx = .1;   %cell dimensions
grid.dy = grid.dy;

%number of int. cells: nx*ny - ie. interior p nodes
grid.nx = grid.Lx/grid.dx;
grid.ny = grid.Ly/grid.dx;

%set up arrays for p,u,v - big arrays including ghost nodes
field.p = zeros(grid.ny+2,grid.nx+2);
field.u = zeros(grid.ny+2,grid.nx+1);
field.v = zeros(grid.ny+1,grid.nx+2);

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
field.u(2:grid.ny+1,1) = param.vIN*sin(param.alpha);
field.v(2:grid.ny,1) = 2*param.vIN*cos(param.alpha) - field.v(2:grid.ny-1,2);
field.p(2:grid.ny+1,1) = field.p(2:grid.ny+1,2);

%South - Neumann
field.u(grid.ny+2,2:grid.nx) = 2*param.vIN*cos(param.alpha)-field.u(grid.ny+1,2:grid.nx);
field.v(grid.ny+1,2:grid.nx+1) = param.vIN*sin(param.alpha);
field.p(grid.ny+2,2:grid.nx+1) = field.p(grid.ny+1,2:grid.nx+1);

%North - far-field
field.v(1,2:grid.nx+1) = 0;
field.u(1,2:grid.nx) = field.u(2,2:grid.nx);
field.p(1,2:grid.nx+1) = field.p(2,2:grid.nx+1);

%East - far-field
field.v(2:grid.ny,grid.nx+2) = field.v(2:grid.ny,grid.nx+1);
field.u(2:grid.ny+1,gird.nx+1) = 0;
field.p(2:grid.ny+1,grid.nx+2) = field.p(2:grid.ny+1,grid.nx+1);


%% Wing - we should maybe test without first(?) - NOT CHANGED YET
    %desired size in 'm'
w.Ly = .25;
w.Lx = .5;
    %size in no cells
w.ldy = floor(w.Ly/grid.dy);
w.ldx = floor(w.Lx/grid.dx);
    %position in  grid
w.idy = ceil((grid.ny/2)-(w.ldy/2));
w.idx = ceil((grid.nx/2)-(w.ldx/2));

%BC for wing - there is probably a more elegenat way to code this
%Note: I am not sure how to properly set values for the corners of 
%the wing, it seems like there will be conflicting p values
%top
field.v(w.idy-1,w.idx:w.idx+w.ldx-1)=0;
field.u(w.idy,w.idx:w.idx+w.ldx-2)= ...
    -field.u(w.idy-1,w.idx:w.idx+w.ldx-2);
field.p(w.idy,w.idx:w.idx+w.ldx-1)=field.p(w.idy-1,w.idx:w.idx+w.ldx-1);
%bot
field.v(w.idy+w.ldy-1,w.idx:w.idx+w.ldx-1)=0;
field.u(w.idy+w.ldy-1,w.idx:w.idx+w.ldx-2)= ...
    -field.u(w.idy+w.ldy-1,w.idx:w.idx+w.ldx-2);
field.p(w.idy+w.ldy-1,w.idx:w.idx+w.ldx-1)=field.p(w.idy+w.ldy,w.idx:w.idx+w.ldx-1);

%left - only v,u for now
field.u(w.idy:w.idy+w.ldy-1,w.idx-1)=0;
field.v(w.idy:w.idy+w.ldy-1,w.idx)=...
    -field.v(w.idy:w.idy+w.ldy-1,w.idx-1);

%right, only v,u for now
field.u(w.idy:w.idy+w.ldy-1,w.idx+w.ldx-1)=0;
field.v(w.idy:w.idy+w.ldy-1,w.idx)=...
    -field.v(w.idy:w.idy+w.ldy-1,w.idx+w.ldx);

%% Parameters (arbitraryly chosen for now)

% Initialization, Parameters
param.Re = 100;
param.rho = 1;
param.mu = 1;
param.T = 5;
param.dt = .1;
param.tsteps = param.T/param.dt;
param.eps = .01;
% etc

%% SIMPLE
%fields are initialized

% TODO (i):
%- agree on matrix structure (sparse/dense, indexing etc.)
%- agree on handling of wing (I am unsure about this)
%- combine fields.[] and BC.[] in complete data object (Z for now)
% to be passed to functions

% (ii) write Z = setBC(Z,param) to set BC at each step (this will be
% similar to what I did above) - BCs also need to be incorporated in
%the solvers, I suspect

%Z = setBC(Z,param);
%Vp = inf;  this is the velocity field (u,v) after each update.
%V(n+1) has to converge to Vp
%Vn = Z.V;
%Pn = Z.P;

%TODO (iii) write iterative solver (eg GS) x = GS(A,b)
% this can be taken from somewhere


%TODO (iv): write function to compute A = -grad(V*V)
    % ie.(for conservative form): Au = -d(u*u)/dx - d(u*v)/dy
    % and Av = -d(u*v)/dx -d(v*v)/dy , both computed at cell centers
    % This involves the blended schemes from lecture 
    %----> A = advection(Z,param);
    
%TODO (v): write function to compute B = (mu/rho)*grad(grad(V))
    % this is just the laplace equation from HWs (+ BCs)
    % (ie.CD twice) 
    
%TODO (vi): write function to compute grad(P). I think CD is fine
    % or even forward/backward
    
%for u, v seperately (I'm gonna use V for both for the comments):
% once here for first step
%A = advection(Z.V,param);
%B = diffusion(Z.V,param);
%dP = grad(Z.P,param);
%Vo=Vn    V from previous step
%Vn = =0; % V at timestep n+1



for i=1:param.tsteps
   
    
   while(abs(Vp-Vn)>param.eps) 
    % calculate V_prime from previous timestep
    %if AB, save additional timestep after while loop
    %Vp = Vo + param.dt*(A+B-grad(Pn)) 
    
    %compute P_prime from grad(grad(Pp))=(1/param.dt)*(grad(vp))
    %I think functions from (iii),(v) &(vi) can be 'overloaded' for this
    %Pp = diffusion((1/param.dt)*(grad(vp))) ~ sth like this
    
    %Pn = Pn + Pp
    %Vn = Vp -param.dt*grad(Pp)  
   end
   
  %solution at timestep i found -> update Z &compute new A,B, gradP
%Z.P = Pn;
%Z.V = Vn;
%A = advection(Z.V,param);
%B = diffusion(Z.V,param);
%dP = grad(Z.P,param);

%output for logging, visualization etc
%TODO (vii): Visualization procedure/ calculation/logging of results

end
































