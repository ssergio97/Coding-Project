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


x = linspace(ax,bx,n) ;
y = linspace(ay,by,n) ;
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

for z = 1:length(t) 
    for i = 2:n-1
   for j = 2:n-1
       dT(i,j) = (T1(i,j+1) - (2*(T1(i,j))) + T1(i,j-1))/(dx^2) +  (T1(i+1,j) - (2*(T1(i,j))) + T1(i-1,j))/(dy^2);
          
    
            
   end 
      T1 =  T1 + (dT*dt) 
    end 
            T1(n,1:n) = C;  %Bottom Condition
            T1(1:n,1) = fby ; %LEFT Condition
            T1(1:n,n) = gby ;  %RIGHT Condition
                             % Ghost Node Neumann Condition 
end


%% 
% Parameters to define the heat equation and the range in space and time
clc
clear all 
close all
n = 10 ; % number of nodes 
ax = 0 ; % Boundary Conditions and given values. 
ay = 0 ; 
bx = 2*pi ;
by = 2*pi ; 

x = linspace(ax,bx,n) ;
y = linspace(ay,by,n) ;
time = 10; % Simulation Time 
dt = .1 ; % Fixed Time Step 
dx = bx/n ; %spacial step in x
dy = by/n ; %spatial step in y
t = 0:dt:time; 
b = dt/(dx*dx); % Parameter of the method

d1(1:n-2) = -b - b;
d2(1:n-1) = 2 + 4*b;
d3(1:n-2) = -b - b; 
Diagonal = inv(diag(d2,0) + diag(d1,-1) + diag(d3,1));

for q = 1:t+1 
    for x = 1:t+1 
    u() 


for z = 2:t 
    dtdt = u(2:n,z-1)
    u(2:n,z) = Diagonal*dtdt 
end 
