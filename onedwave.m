%1D wave equation simulation
clear all
close all

%parameters:
T = 1;          %time duration
Fs = 44100;     %sample rate in Hz
k = 1/Fs;       %number of time steps
Ns = floor(T/k);%number of samples
c = 200;        %wave speed (orders of 100)
h = 2*c*k;      %grid spacing (>ck=0.0023)
                %small h -> more accurate soln
L = 1;          %length of string in metres
N = floor(L/h); %number of grid segments
                %(there are N+1 grid points)
h = L/N;        %redefining h so that it matches with our N
                
%scheme coefficients for loop:
c1 = 2-(2*(c*k)^2)/h^2;
c2 = ((c*k)/h)^2;

%initialise 2 vectors of length N+1 (initial values):
u0 = zeros(N+1,1);
u1 = zeros(N+1,1);
%setting some initial displacement:
u0(round(N/2):round(9+N/2)) = hann(10); 
u1(round(N/2):round(9+N/2)) = hann(10); 

%create matrix of initial values instead?
u = zeros(N+1,1);
out = zeros(Ns, 1); %output vector 
                    %(like the input of a mic)

%we need two loops, one for time steps
%and one for grid points

for n = 1:Ns %time steps 
    %imposing boundary conditions:
    u(1) = 0; 
    u(N+1) = 0;
    for l = 2:N %grid spacing 
        %u(l,n+1) = c1*u(l,n) - u(l,n-1) + c2*(u(l+1,n) + u(l-1,n));
        u(l) = c1*u1(l) - u0(l)+ c2*(u1(l+1)+u1(l-1));     
    end
    %updating initial conditions
    u0 = u1; 
    u1 = u;
    
    %taking one value of u as our output
    out(n) = u(round(N/2)); 
    
    plot(u);
    ylim([-1,1]); xlim([0, N+1]);
    drawnow
end