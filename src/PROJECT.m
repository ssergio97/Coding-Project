%% Sergio Esperanza 
% 1404790 
% Explicit Discretization Method for 2 Dimensional Diffusion Equation 

clc
clear all 
close all 

% Domain of Interest Parameters and Boundary Conditions 

n = 10 ; % number of nodes 
ax = 0 ; % Boundary Conditions and given values. 
ay = 0 ; 
bx = 2*pi ;
by = 2*pi ; 


x = linspace(ax,bx,n) ; % x Axis 
y = linspace(ay,by,n) ; % Y Axis 
time = 1; % Simulation Time 
dt = .1 ; % Time Step 
dx = bx/n ; %spacial step in x
dy = by/n ; %spatial step in y
t = 0:dt:time; 

% defining fb(y) and gb(y)
for q = 1:n 
    fby(q) = y(q)*(by - y(q))^2 ;
    gby(q) = ((by - y(q))^2)*cos((pi*y(q))/by);
end
C = by*x; % Reducing for the u(x,y=ay) boundary condition gives this simplified expression

T1 = ones(n,n) ; %Initiliazing


dT = zeros(n,n); %Initiliazing

% explicit Discretization 
for z = 1:length(t) 
    for i = 2:n-1
   for j = 2:n-1
       dT(i,j) = (T1(i,j+1) - (2*(T1(i,j))) + T1(i,j-1))/(dx^2) +  (T1(i+1,j) - (2*(T1(i,j))) + T1(i-1,j))/(dy^2);
                 
   end 
      T1 =  (dT*dt); 
    end 
            T1(n,1:n) = C;  %Bottom Condition
            T1(1:n,1) = fby ; %LEFT Condition
            T1(1:n,n) = gby ;  %RIGHT Condition
            T1(i,1) = 0;       % Ghost Node Neumann Condition 
end

figure(3)
mesh(T1)

figure(4)
plot(T1,'DisplayName','T1') 

%% Implicit Time Integration for 2 Dimensional Diffusion Equation 
clc
clear all 
close all 


n = 10 ; % number of nodes 
ax = 0 ; % Boundary Conditions and given values. 
ay = 0 ; 
bx = 2*pi ;
by = 2*pi ; 


time = 1; % Simulation Time 
dt = .1 ; % Fixed Time Step 
dx = bx/n ; %spacial step in x
dy = by/n ; %spatial step in y
t = 0:dt:time; 
x = linspace(ax,bx,n) ;
y = linspace(ay,by,n) ;
b = dt/(dx*dx); % Parameter 
% Boundary Condition Equations
for w = 1:n 
    fby(w) = y(w)*(by - y(w))^2 ;
    gby(w) = ((by - y(w))^2)*cos((pi*y(w))/by);
end

C = by*x; % Reducing for the u(x,y=ay) boundary condition gives this simplified expression

%Setting up tridiagonal matrix
d1(1:n-2) = -b - b;
d2(1:n-1) = 2 + b+ b+ b+b;
d3(1:n-2) = -b - b; 
Diagonal = inv(diag(d2,0) + diag(d1,-1) + diag(d3,1));

%Implementing boundary conditions 
u(n,1:n) = C;  %Bottom Condition
u(1:n,1) = fby ; %LEFT Condition
u(1:n,n) = gby ;  %RIGHT Condition

%Implicit loop 
for z = 2:length(t)
    dtdt = u(2:n,z-1);
    u(2:n,z) = Diagonal*dtdt ;
      
end 

figure(3)
mesh(u(:,1:10))

figure(4)
plot(u(:,1:10),'DisplayName','u(:,1:10)')