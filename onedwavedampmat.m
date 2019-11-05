%1D wave equation simulation with damping 
%nice version 
clear all
close all

%parameters:
T = 1;          %time duration
Fs = 44100;     %sample rate in Hz
f0 = 10000;      %fundamental frequency 
k = 1/Fs;       %number of time steps
Ns = floor(T/k);%number of samples
L = 1;          %length of string in metres
c = 200;        %wave speed
%c = 2*f0*L;     %wave speed
h = 2*c*k;      %grid spacing (>ck=0.0023)
                %small h -> more accurate soln
N = floor(L/h); %number of grid segments
                %(there are N+1 grid points)
h = L/N;        %redefining h so that it matches with our N
b = 0.01;       %damping coefficient

%scheme matrices for loop:
b1 = 2-2*((c*k)/h)^2;
b2 = ((c*k)/h)^2;

A = (b/2*k + 1)*eye(N+1, N+1);
C = (b/2*k - 1)*eye(N+1, N+1);
AI = inv(A); 

Ones = ones(N+1, 1) ;
B = diag(b1 * Ones, 0) + diag(b2*Ones(1:N), -1) + diag(b2*Ones(1:N), 1);

%initialise 2 vectors of length N+1 (initial values):
u0 = zeros(N+1,1);
u1 = zeros(N+1,1);
%setting some initial displacement:
u0(round(N/2 - 20):round(N/2)) = 0.7; 
u1(round(N/2 - 20):round(N/2)) = 0.7; 

%create matrix of initial values instead?
%u = zeros(N+1,1);
out = zeros(Ns, 1); %output vector 
                    %(like the input of a mic)

%we need two loops, one for time steps
%and one for grid points

for n = 1:Ns %time steps 
    %imposing boundary conditions:
    u(1) = 0; 
    u(N+1) = 0;
    %computing new u
    u = A\(B*u1 + C*u0);     
    %updating initial conditions
    u0 = u1; 
    u1 = u;
    %taking one value of u as our output
    out(n) = u(round(N/2)); 
    
    plot(u);
    ylim([-1,1]); xlim([0, round(N)]);
    drawnow
end