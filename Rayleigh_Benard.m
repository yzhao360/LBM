clc;clear all;close all;
% define numerical size of domain
N=128; % grid sizes
r=4; % ratio of Long-side/Short-side
nx=r*N; W=nx;% long side of the domain-horizontal
ny=N; H=ny;% short side of the domain-vertical

% temperature features
T1=0; % temperature on the upper boundary
T2=1; % temperature on the lower boundary
dT=T2-T1; 


% define weight coefficient(D2Q9)
w0=4/9;
w1=1/9; w2=1/9; w3=1/9; w4=1/9;
w5=1/36; w6=1/36; w7=1/36; w8=1/36;

% initialize variable values in the field
c=1; %lattice speed
dt=1;% delta t
% for D2Q9 model-fluids part
f=1/9*ones(ny+1,nx+1,9);
feq=f;
% for D2Q5 model-temperature part
g=zeros(ny+1,nx+1,5); % temperature distribution function
geq=g;
% upper bound
g(1,:,1)=T1*ones(1,nx+1)/6;
g(1,:,2)=T1*ones(1,nx+1)/6;
g(1,:,3)=T1*ones(1,nx+1)/6;
g(1,:,4)=T1*ones(1,nx+1)/6;
g(1,:,5)=T1*ones(1,nx+1)/3;
% lower bound
g(ny+1,:,1)=T2*ones(1,nx+1)/6;
g(ny+1,:,2)=T2*ones(1,nx+1)/6;
g(ny+1,:,3)=T2*ones(1,nx+1)/6;
g(ny+1,:,4)=T2*ones(1,nx+1)/6;
g(ny+1,:,5)=T2*ones(1,nx+1)/3;
% perturbation needed to trigger the convection
g(ny,nx/2+1,1)=1.1*T2/6;
g(ny,nx/2+1,2)=1.1*T2/6;
g(ny,nx/2+1,3)=1.1*T2/6;
g(ny,nx/2+1,4)=1.1*T2/6;
g(ny,nx/2+1,5)=1.1*T2/3;
% important dimensionless number
Ra=5e4; % Raleigh number
Pr=0.71;
gr=1e-3; % gravitational acceleration( acting like a stability coeff)
beta=1;
nu=sqrt((Pr/Ra)*gr*beta*dT*H^3); %kinematic viscosity of the fluids
tau=3*nu+0.5; % single relaxation time for fluids
kappa=beta*gr*dT*H^3/Ra/nu;% heat diffusivity coefficient
tau2=2*kappa+0.5; % single relaxation time for heat transfer

domy=2:ny;

% myobj=VideoWriter('convection.avi');
% writerObj.FrameRate=30;
% open(myobj);

for it=1:50000

    %macroscopic properties
    rho=sum(f,3);
    T=sum(g,3);
    u=(sum(f(:,:,[1 5 8]),3)-sum(f(:,:,[3 6 7]),3))./rho;
    v=(sum(f(:,:,[2 5 6]),3)-sum(f(:,:,[4 7 8]),3))./rho;
    %collision
    % calculate temperature equilibrium distri-function
    geq(:,:,1)=T.*(1+3*u/c)./6;
    geq(:,:,2)=T.*(1+3*v/c)./6;
    geq(:,:,3)=T.*(1-3*u/c)./6;
    geq(:,:,4)=T.*(1-3*v/c)./6; 
    geq(:,:,5)=T./3.*(1-0);
    % calculate equilibrium distribution
    feq(:,:,1)=w1*rho.*(1+3*u/c+9/2*u.^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,2)=w2*rho.*(1+3*v/c+9/2*v.^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,3)=w3*rho.*(1+3*-u/c+9/2*u.^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,4)=w4*rho.*(1+3*-v/c+9/2*v.^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,5)=w5*rho.*(1+3*(u+v)/c+9/2*(u+v).^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,6)=w6*rho.*(1+3*(-u+v)/c+9/2*(-u+v).^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,7)=w7*rho.*(1+3*(-u-v)/c+9/2*(-u-v).^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,8)=w8*rho.*(1+3*(u-v)/c+9/2*(u-v).^2/c^2-3/2*(u.^2+v.^2)/c^2);
    feq(:,:,9)=w0*rho.*(1-3/2*(u.^2+v.^2)/c^2); 
    %calculate body forces
    T0=mean(T(:));
    G=gr*beta*rho.*(T-T0);
    F(:,:,1)=zeros(ny+1,nx+1);
    F(:,:,2)=3*w2.*G;
    F(:,:,3)=zeros(ny+1,nx+1);
    F(:,:,4)=3*w4.*-G;
    F(:,:,5)=3*w5.*G;
    F(:,:,6)=3*w6.*G;
    F(:,:,7)=3*w7.*-G;
    F(:,:,8)=3*w8.*-G;
    F(:,:,9)=zeros(ny+1,nx+1);
    % collision process
    g=g-1/tau2.*(g-geq);
    f=f-1/tau.*(f-feq)+F;        

    % streaming process for temperature
    g(:,:,1)=circshift(g(:,:,1),[0,1]);
    g(:,:,2)=circshift(g(:,:,2),[-1,0]); 
    g(:,:,3)=circshift(g(:,:,3),[0,-1]);
    g(:,:,4)=circshift(g(:,:,4),[1,0]);
    % streaming process for fluid particles
    f(:,:,1)=circshift(f(:,:,1),[0,1]); 
    f(:,:,2)=circshift(f(:,:,2),[-1,0]); 
    f(:,:,3)=circshift(f(:,:,3),[0,-1]);
    f(:,:,4)=circshift(f(:,:,4),[1,0]);
    f(:,:,5)=circshift(f(:,:,5),[-1,1]);
    f(:,:,6)=circshift(f(:,:,6),[-1,-1]);
    f(:,:,7)=circshift(f(:,:,7),[1,-1]);
    f(:,:,8)=circshift(f(:,:,8),[1,1]);
    
    % add BC
    % top constant T
    g(1,:,4)=T1*ones(1,nx+1)-sum(g(1,:,[1 2 3 5]),3);
    %bottom constant T
    g(ny+1,:,2)=T2*ones(1,nx+1)-sum(g(ny+1,:,[1 3 4 5]),3);
    % left periodic BC
    f(domy,1,1)=f(domy,nx+1,1);
    f(domy,1,5)=f(domy,nx+1,5);
    f(domy,1,8)=f(domy,nx+1,8);
    g(domy,1,1)=g(domy,nx+1,1);
    % right periodic BC
    f(domy,nx+1,3)=f(domy,1,3);
    f(domy,nx+1,6)=f(domy,1,6);
    f(domy,nx+1,7)=f(domy,1,7);
    g(domy,nx+1,3)=g(domy,1,3); 
    % upper boundary
    % bounce back for f
    f(1,:,4)=f(1,:,2);
    f(1,:,7)=f(1,:,5);
    f(1,:,8)=f(1,:,6);
    % lower boundary
    % bounce back for f
    f(ny+1,:,2)=f(ny+1,:,4);
    f(ny+1,:,5)=f(ny+1,:,7); 
    f(ny+1,:,6)=f(ny+1,:,8);
 
    if rem(it,5)==0
        subplot(2,1,1)
        imagesc(T);
        fff=getframe;      
%         writeVideo(myobj,fff);
    end   
    
    
end
close(myobj);