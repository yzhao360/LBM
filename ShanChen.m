clear all;clc;close all;
% define numerical parameters
N=256;
nx=1*N; ny=N;
tau1=1; tau2=1; % relexation time
G=-1.5;
% define weight coefficient(D2Q9)
w0=4/9;
w1=1/9; w2=1/9; w3=1/9; w4=1/9;
w5=1/36; w6=1/36; w7=1/36; w8=1/36;
% initialize variable values in the field
c=1; %lattice speed
dt=1;% delta t
% initialize distribution functions for two components
delta_rho=0.001*(1-2*rand(ny+1,nx+1));
rho1=1+delta_rho;
rho2=1-delta_rho;
% distribution function for component 1
f(:,:,1)=w1*rho1;
f(:,:,2)=w2*rho1;
f(:,:,3)=w3*rho1;
f(:,:,4)=w4*rho1;
f(:,:,5)=w5*rho1;
f(:,:,6)=w6*rho1;
f(:,:,7)=w7*rho1;
f(:,:,8)=w8*rho1;
f(:,:,9)=w0*rho1;
% distribution function for component 2
g(:,:,1)=w1*rho2;
g(:,:,2)=w2*rho2;
g(:,:,3)=w3*rho2;
g(:,:,4)=w4*rho2;
g(:,:,5)=w5*rho2;
g(:,:,6)=w6*rho2;
g(:,:,7)=w7*rho2;
g(:,:,8)=w8*rho2;
g(:,:,9)=w0*rho2;

for it=1:10000
    % macropic properties
    % calculate interaction body forces
    rho1=sum(f,3); %density of fluid 1
    rho2=sum(g,3); %density of fluid 2
    rho_tot=rho1+rho2;% total local density
    % body forces
    F221_x=-rho1.*G.*(w1*circshift(rho2,[0,1])-w3*circshift(rho2,[0,-1])+w5*circshift(rho2,[-1,1])-w6*circshift(rho2,[-1,-1])...
        -w7*circshift(rho2,[1,-1])+w8*circshift(rho2,[1,1]));
    F221_y=-rho1.*G.*(w2*circshift(rho2,[-1 0])-w4*circshift(rho2,[1,0])+w5*circshift(rho2,[-1,1])+w6*circshift(rho2,[-1,-1])...
        -w7*circshift(rho2,[1,-1])-w8*circshift(rho2,[1,1]));
    F122_x=-rho2.*G.*(w1*circshift(rho1,[0,1])-w3*circshift(rho1,[0,-1])+w5*circshift(rho1,[-1,1])-w6*circshift(rho1,[-1,-1])...
        -w7*circshift(rho1,[1,-1])+w8*circshift(rho1,[1,1]));
    F122_y=-rho2.*G.*(w2*circshift(rho1,[-1 0])-w4*circshift(rho1,[1,0])+w5*circshift(rho1,[-1,1])+w6*circshift(rho1,[-1,-1])...
        -w7*circshift(rho1,[1,-1])-w8*circshift(rho1,[1,1]));
    % velocity field
    u=(sum(f(:,:,[1 5 8]),3)-sum(f(:,:,[3 6 7]),3)+sum(g(:,:,[1 5 8]),3)-sum(g(:,:,[3 6 7]),3))./rho_tot+(F221_x+F122_x)./rho_tot./2;
    v=(sum(f(:,:,[2 5 6]),3)-sum(f(:,:,[4 7 8]),3)+sum(g(:,:,[2 5 6]),3)-sum(g(:,:,[4 7 8]),3))./rho_tot+(F221_y+F122_y)./rho_tot./2;
    
    % collision step
    % calculate equilibrium distribution for fluid 1
    feq(:,:,1)=w1*rho1.*(1+3*u/c+9/2*u.^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,2)=w2*rho1.*(1+3*v/c+9/2*v.^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,3)=w3*rho1.*(1+3*-u/c+9/2*u.^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,4)=w4*rho1.*(1+3*-v/c+9/2*v.^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,5)=w5*rho1.*(1+3*(u+v)/c+9/2*(u+v).^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,6)=w6*rho1.*(1+3*(-u+v)/c+9/2*(-u+v).^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,7)=w7*rho1.*(1+3*(-u-v)/c+9/2*(-u-v).^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,8)=w8*rho1.*(1+3*(u-v)/c+9/2*(u-v).^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,9)=w0*rho1.*(1-3/2*(u.^2+v.^2)/c^2);
    % for fluid 2
    geq(:,:,1)=w1*rho2.*(1+3*u/c+9/2*u.^2/c^2-3/2*(u.^2+v.^2)/c^2);
    geq(:,:,2)=w2*rho2.*(1+3*v/c+9/2*v.^2/c^2-3/2*(u.^2+v.^2)/c^2);
    geq(:,:,3)=w3*rho2.*(1+3*-u/c+9/2*u.^2/c^2-3/2*(u.^2+v.^2)/c^2);
    geq(:,:,4)=w4*rho2.*(1+3*-v/c+9/2*v.^2/c^2-3/2*(u.^2+v.^2)/c^2);
    geq(:,:,5)=w5*rho2.*(1+3*(u+v)/c+9/2*(u+v).^2/c^2-3/2*(u.^2+v.^2)/c^2);
    geq(:,:,6)=w6*rho2.*(1+3*(-u+v)/c+9/2*(-u+v).^2/c^2-3/2*(u.^2+v.^2)/c^2);
    geq(:,:,7)=w7*rho2.*(1+3*(-u-v)/c+9/2*(-u-v).^2/c^2-3/2*(u.^2+v.^2)/c^2);
    geq(:,:,8)=w8*rho2.*(1+3*(u-v)/c+9/2*(u-v).^2/c^2-3/2*(u.^2+v.^2)/c^2);
    geq(:,:,9)=w0*rho2.*(1-3/2*(u.^2+v.^2)/c^2);  
    %calculate body force terms in collision step
    % for fluid 1
    F1(:,:,1)=w1*(1-1/2/tau1)*((6*u+3).*F221_x+-3*v.*F221_y);
    F1(:,:,2)=w2*(1-1/2/tau1)*(-3*u.*F221_x+(3+6*v).*F221_y);
    F1(:,:,3)=w3*(1-1/2/tau1)*((6*u-3).*F221_x+-3*v.*F221_y);
    F1(:,:,4)=w4*(1-1/2/tau1)*(-3*u.*F221_x+(-3+6*v).*F221_y);
    F1(:,:,5)=w5*(1-1/2/tau1)*((3+6*u+9*v).*F221_x+(3+9*u+6*v).*F221_y);
    F1(:,:,6)=w6*(1-1/2/tau1)*((-3+6*u-9*v).*F221_x+(3-9*u+6*v).*F221_y);
    F1(:,:,7)=w7*(1-1/2/tau1)*((-3+6*u+9*v).*F221_x+(-3+9*u+6*v).*F221_y);
    F1(:,:,8)=w8*(1-1/2/tau1)*((3+6*u-9*v).*F221_x+(-3-9*u+6*v).*F221_y);
    F1(:,:,9)=w0*(1-1/2/tau1)*(-3*u.*F221_x+-3*v.*F221_y);
    % for fluid 2
    F2(:,:,1)=w1*(1-1/2/tau2)*((6*u+3).*F122_x+-3*v.*F122_y);
    F2(:,:,2)=w2*(1-1/2/tau2)*(-3*u.*F122_x+(3+6*v).*F122_y);
    F2(:,:,3)=w3*(1-1/2/tau2)*((6*u-3).*F122_x+-3*v.*F122_y);
    F2(:,:,4)=w4*(1-1/2/tau2)*(-3*u.*F122_x+(-3+6*v).*F122_y);
    F2(:,:,5)=w5*(1-1/2/tau2)*((3+6*u+9*v).*F122_x+(3+9*u+6*v).*F122_y);
    F2(:,:,6)=w6*(1-1/2/tau2)*((-3+6*u-9*v).*F122_x+(3-9*u+6*v).*F122_y);
    F2(:,:,7)=w7*(1-1/2/tau2)*((-3+6*u+9*v).*F122_x+(-3+9*u+6*v).*F122_y);
    F2(:,:,8)=w8*(1-1/2/tau2)*((3+6*u-9*v).*F122_x+(-3-9*u+6*v).*F122_y);
    F2(:,:,9)=w0*(1-1/2/tau2)*(-3*u.*F122_x+-3*v.*F122_y);
    % collision
    f=f-1/tau1.*(f-feq)+F1;
    g=g-1/tau2.*(g-geq)+F2;
    
    %streaming
    f(:,:,1)=circshift(f(:,:,1),[0,1]); 
    f(:,:,2)=circshift(f(:,:,2),[-1,0]); 
    f(:,:,3)=circshift(f(:,:,3),[0,-1]);
    f(:,:,4)=circshift(f(:,:,4),[1,0]);
    f(:,:,5)=circshift(f(:,:,5),[-1,1]);
    f(:,:,6)=circshift(f(:,:,6),[-1,-1]);
    f(:,:,7)=circshift(f(:,:,7),[1,-1]);
    f(:,:,8)=circshift(f(:,:,8),[1,1]);
    
    g(:,:,1)=circshift(g(:,:,1),[0,1]); 
    g(:,:,2)=circshift(g(:,:,2),[-1,0]); 
    g(:,:,3)=circshift(g(:,:,3),[0,-1]);
    g(:,:,4)=circshift(g(:,:,4),[1,0]);
    g(:,:,5)=circshift(g(:,:,5),[-1,1]);
    g(:,:,6)=circshift(g(:,:,6),[-1,-1]);
    g(:,:,7)=circshift(g(:,:,7),[1,-1]);
    g(:,:,8)=circshift(g(:,:,8),[1,1]);
    % add BC
    %L
    f(:,1,1)=f(:,nx+1,1);
    f(:,1,5)=f(:,nx+1,5);
    f(:,1,8)=f(:,nx+1,8);
    g(:,1,1)=g(:,nx+1,1);
    g(:,1,5)=g(:,nx+1,5);
    g(:,1,8)=g(:,nx+1,8);   
    %R
    f(:,nx+1,3)=f(:,1,3);
    f(:,nx+1,6)=f(:,1,6);
    f(:,nx+1,7)=f(:,1,7);
    g(:,nx+1,3)=g(:,1,3);
    g(:,nx+1,6)=g(:,1,6);
    g(:,nx+1,7)=g(:,1,7);
    %T
    f(1,:,4)=f(ny+1,:,4);
    f(1,:,7)=f(ny+1,:,7);
    f(1,:,8)=f(ny+1,:,8);
    g(1,:,4)=g(ny+1,:,4);
    g(1,:,7)=g(ny+1,:,7);
    g(1,:,8)=g(ny+1,:,8);
    %B
    f(ny+1,:,2)=f(1,:,2);
    f(ny+1,:,5)=f(1,:,5);
    f(ny+1,:,6)=f(1,:,6);
    g(ny+1,:,2)=g(1,:,2);
    g(ny+1,:,5)=g(1,:,5);
    g(ny+1,:,6)=g(1,:,6);    
    
    if rem(it,5)==0       
        imagesc(rho1);
        fff=getframe;
        axis equal off;
    end
end





