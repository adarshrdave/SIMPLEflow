%% F18 CFD Project - incompressible flow around wing
%TLDR: We should agree on/discuss indexing,BC/grid setup. I've started
%with the grid and BC equations, but it's already kinda messy, and 
%I have some questions about the BCs.
%Then we can split up on functions, which should work well here,
%I think. see TODOs in SIMPLE section
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
grid.Lx = 2;
grid.Ly = 1;
grid.dx = .1;
grid.dy = .1;

%number of int. cells: nx*ny - ie. interior p nodes
grid.nx = grid.Lx/grid.dx;
grid.ny = grid.Ly/grid.dx;

%set up arrays for p,u,v - only interior
field.p = zeros(grid.ny,grid.nx);
field.u = zeros(grid.ny,grid.nx-1);
field.v = zeros(grid.ny-1,grid.nx);

% BCs - far-field for N & E, Neumann for W & S,
%       Dirchilet at wing
% Note: I don't know how these BCs will behave in the NW & SE corners

% Inflow across all of W & S
param.alpha = pi/3;
param.vIN = 5; % magnitude of velocity of inflow

%W vectors for ghost nodes at W, 
%starting at cell P(1,0) until cell P(ny,0)(v not used)
BC.uW = zeros(grid.ny,1);
BC.uW(:) = param.vIN*sin(param.alpha);
BC.vW = zeros(grid.ny,1);
BC.vW(:) = -[field.v(:,1); zeros(1,1)] + 2*param.vIN*cos(param.alpha); 
BC.pW = zeros(grid.ny,1);
BC.pW(:) = field.p(:,1);

%S vectors for ghost nodes at S
%starting at cell P(ny+1,1)(not used) until cell P(ny+1,nx)(u not used)
BC.uS = zeros(1,grid.nx);
BC.uS(:) = -[field.u(grid.ny,:) zeros(1,1)] + 2*param.vIN*cos(param.alpha); 
BC.vS = zeros(1,grid.nx);
BC.vS(:) = param.vIN*cos(param.alpha);
BC.pS = zeros(1,grid.nx);
BC.pS(:) = field.p(grid.ny,:);


%far-field at N&E 

%N vectors for ghost nodes at N
%starting cell P(0,1)until P(0,nx)(u not used)
BC.uN = zeros(1,grid.nx);
BC.uN(:) = [field.u(1,:) zeros(1,1)];
BC.vN = zeros(1,grid.nx);
BC.vN(:) = 0;
BC.pN = zeros(1,grid.nx);
BC.pN(:) = field.p(1,:);

%E vectors for ghost nodes at E
%starting cell P(nx+1,1)until P(nx+1,ny)(v not used)
BC.uE = zeros(grid.ny,1);
BC.uE(:) = 0;
BC.vE = zeros(grid.ny,1);
BC.vE(:) =  [field.v(:,grid.nx); zeros(1,1)];
BC.pE = zeros(grid.ny,1);
BC.pE(:) = field.p(:,grid.nx);

%Note: at the moment, ghost cell vectors are one element too short
%to set up the field with a complete ghost shell.
%If that is needed for the algorithm (I'm not sure atm), additional
%dummy elements need to be inserted to set up the full (ny+2)*(nx+2) 
%matrix

%% Wing - we should maybe test without first(?)
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
































