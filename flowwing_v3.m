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
% GRID
grid.Lx = 16;    %x-length of box
grid.Ly = 8;    %y-length of box
grid.dx = .1;   %cell dimensions
grid.dy = grid.dx;

% WING
w.Ly = 0.4;
w.Lx = 4;

%PARAMETERS
param.alpha = (55/180)*(pi);
param.vIN = 1; 
param.Re = 100;
param.rho = 1.184;
param.mu = 1;
param.T = 10;
param.dt = .01;
param.tsteps = param.T/param.dt;
param.eps = .01;
param.Q=0.5;
param.s=10; %maximum speed in m/s

%------------------INITIALIZATION--------------------------------
%number of int. cells: nx*ny - ie. interior p nodes
grid.nx = ceil(grid.Lx/grid.dx);
grid.ny = ceil(grid.Ly/grid.dy);

%set up arrays for p,u,v - big arrays including ghost nodes
field.p = zeros(grid.nx+2,grid.ny+2);
field.u = zeros(grid.nx+1,grid.ny+2);
field.v = zeros(grid.nx+2,grid.ny+1);

%size of wing in no cells
w.ldy = floor(w.Ly/grid.dy);
w.ldx = floor(w.Lx/grid.dx);
    %position in  grid - left-bottom corner-cell of wing
w.idy = ceil((grid.ny/2)-(w.ldy/2));
w.idx = ceil((grid.nx/2)-(w.ldx/2));

% Init BCs
field = setBC(field,grid,param,w);
%-----------------------------------------------------------

%% Steady State run
% determine T s. t. field relatively constant ~ check diff

for t=1:param.tsteps

field = setBC(field,grid,param,w);
p_old=field.p;

field = computeVF_sections(field,grid,param,w);

field.p = solveP_GS(field,grid,param,w);

field = updateV(field,grid,param,w);

%some outputs - comment/uncomment, set frequency
if (mod(t,10)==0)
%       figure(1)
%      axis tight manual
%      visualize(field,grid,param,w)

%rel. change in p field if wanted
% diff = norm(p_old-field.p)/norm(p_old);
% disp(diff);
% timestep
disp(t*param.dt)
end
end
visualize(field,grid,param,w);
[drag,cd, lift,cl] = solve_DragLift2(field,grid,param,w);

%save field for reinitilization
field_ss=field;
%% Varying angle - not really needed anymore
for dummy=1:1
% bring to steady state at param.alpha=alpha_start first.
field = field_ss;
alpha_start = pi/2;
alpha_end = pi/4;
d = 2;
param.T = .5; %for each angle
param.dt = .01;
param.tsteps = 500;
ax = gca;
ax.NextPlot = 'replaceChildren';
F2(param.tsteps) = struct('cdata',[],'colormap',[]);
x =1;
for t=1:param.tsteps
field = setBC(field,grid,param,w);
field = computeVF_sections(field,grid,param,w);
field.p = solveP_GS(field,grid,param,w);
field = updateV(field,grid,param,w);
if(mod(t,10)==0)
     visualize(field,grid,param,w);
     drawnow
     F2((di-1)*2 +x) = getframe; 
     x=x+1;
end

end
end
%% test space - some good plots in here, don't change
for dummy = 1:1
% figure(1)
% plot(lifts2)
% figure(2)
% plot(drags2)
% [drag,cd, lift,cl] = solve_DragLift2(field,grid,param,w);
% figure(1)
% axis tight manual
% visualize(field,grid,param,w);
% figure(1)
% plot(lifts)
% figure(2)
% plot(drags)
% v = VideoWriter('anmtn3_fp8.avi','Motion JPEG AVI');
% v.FrameRate = 8;
% open(v)
% writeVideo(v,F2)
% close(v)
% figure(1)
% l = [1 2 4 6 8];
% lft = [0.3 .5 5.1 9.4 13];
% drg = [ 0.5 0.8 3.5 5.9 8.5];
% plot(l,lft,l,drg)
%  legend('lift','drag','Location','southeast')
%  xlabel('length of wing in m')
%  ylabel('lift/drag force in N')
% figure(1);hold on;
% alph1 = [0 10 20 30 35 45 50 55 60 75];
% drg1 = [0 0.488 1.54 2.8 3.4 4.4 4.64 5.0 5.6 7.6];
% lft1 = [0 2.29 4.15 4.9 5.0 4.8 4.48 4.1 3.7 2.15];
% plot(alph1,lft1,'-r+')
% plot(alph1,drg1,'-bo')
% legend('lift','drag','Location','northwest')
% title('freestream velocity: 1 m/s')
% % hold off
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% set(findall(gca, 'Type', 'Line'),'MarkerSize',8);
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
% figure(2);hold on;
% alph2 = [0 5 10 12 15 20 30];
% drg2 = [0 0.4 1.5 2.0 2.5 2.9 4.4];
% lft2 = [0 3.5 6.4 6.8 6.14 3.9 3.5];
% plot(alph2,lft2,'-r+')
% plot(alph2,drg2,'-b*')
% legend('lift','drag','Location','southeast')
% title('freestream velocity: 2 m/s')
% xlabel('angle of attack in deg.')
% ylabel('lift/drag force in N')
% hold off
% %axis([0 75 0 8])
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% set(findall(gca, 'Type', 'Line'),'MarkerSize',8);
% set(findall(gcf,'-property','FontSize'),'FontSize',12)
end
%% Functions
function[dP_dy] =  dPdy(field,grid)
dP_dy = zeros(size(field.v));
for i=2:grid.nx+1
    for j=2:grid.ny
        dP_dy(i,j) = (field.p(i,j+1)-field.p(i,j))/grid.dy;
    end
end
end

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

function[] = visualize(field,grid,param,w)
[x,y] = meshgrid(1:1:grid.nx,1:1:grid.ny);
u = zeros(grid.nx,grid.ny);
v = zeros(grid.nx,grid.ny);
p = zeros(grid.nx,grid.ny);
for i = 2:grid.nx+1
    for j = 2:grid.ny+1
        u(i-1,j-1) = (field.u(i,j)+field.u(i-1,j))/2;
        v(i-1,j-1) = (field.v(i,j)+field.v(i,j-1))/2;  
        p(i-1,j-1) = field.p(i,j);
    end
end
u(w.idx-1:w.idx+w.ldx-2,w.idy-1:w.idy+w.ldy-2)=0;
v(w.idx-1:w.idx+w.ldx-2,w.idy-1:w.idy+w.ldy-2)=0;

hold on;
clf
surface(p','EdgeColor','none','LineStyle','none','FaceLighting','phong')
xlabel(['x in cells (dx= ',num2str(grid.dx),' m)'])
ylabel(['y in cells (dy= ',num2str(grid.dy),' m)'])
c = colorbar;
c.Label.String = 'pressure (in kPa)';
% contour(p',)
set(c, 'ylim', [1 5])
h=streamslice(x,y,u.',v.',5);
set(h,'color','w')
for i=1:length(h)
    zi = interp2(p',get(h(i),'xdata'),get(h(i),'ydata'));
    set(h(i),'zdata',zi);
end
adeg = 90-ceil(((param.alpha/pi) *180));
title(['angle of attack: ',num2str(adeg),' degrees'])
axis([0 grid.nx 0 grid.ny])
hold off
end

%TODO: add shear drag
function[drag,Cd,lift,Cl] = solve_DragLift2(field,grid,param,w)
drag = 0;
lift = 0;
% bottom
for i = w.idx:w.idx+w.ldx-1 % w.idx:w.idx+w.ldx-1,w.idy
    drag = drag + field.p(i,w.idy-1)*cos(param.alpha)*grid.dx;
    lift = lift + field.p(i,w.idy-1)*sin(param.alpha)*grid.dx;
end

% top
for i = w.idx:w.idx+w.ldx-1 % w.idx:w.idx+w.ldx-1,w.idy+w.ldy-1
    drag = drag - field.p(i,w.idy+w.ldy)*cos(param.alpha)*grid.dx;
    lift = lift - field.p(i,w.idy+w.ldy)*sin(param.alpha)*grid.dx;
end

% left
for j = w.idy:w.idy+w.ldy-1 % w.idx,w.idy:w.idy+w.ldy-1
    drag = drag + field.p(w.idx-1,j)*cos(param.alpha)*grid.dy;
    lift = lift - field.p(w.idx-1,j)*sin(param.alpha)*grid.dy;
end

% right
for j = w.idy:w.idy+w.ldy-1 % w.idx+w.ldx-1,w.idy:w.idy+w.ldy-1
    drag = drag - field.p(w.idx+w.ldx,j)*cos(param.alpha)*grid.dy;
    lift = lift + field.p(w.idx+w.ldx,j)*sin(param.alpha)*grid.dy;
end
Cd = (2*drag)/(param.rho*(param.vIN^2)*w.Lx);
Cl = (2*lift)/(param.rho*(param.vIN^2)*w.Lx);

end

%% deprecated
function[] = makeGIF(field,grid,w,param,a,filename2)
% a=param.alpha;
fig1 = figure(1);
% adeg = ceil(((a/pi) *180));
% title(['Inflow angle: ',num2str(adeg),' degree'])
visualize(field,grid,param,w);
drawnow  
frame2 = getframe(fig1); 
im2 = frame2im(frame2); 
[imind2,cm2] = rgb2ind(im2,256); 
% Write to the GIF File 
if a == 1
  imwrite(imind2,cm2,filename2,'gif', 'Loopcount',inf); 
else 
imwrite(imind2,cm2,filename2,'gif','WriteMode','append'); 
end 
end

