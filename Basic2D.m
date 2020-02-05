% Plate and Membrane Modelling

Fs = 20000;        % Sample rate
k = 1/Fs;          % Time step
c = 500;           % Wave speed
epsilon = 1;       % Aspect ratio of plate, 1 for a square
T = 1;             % End time
Lx = 1;            % Length of x axis
Ly = 1;            % Length of y axis

h = c*k;           
Nx = floor(Lx/h);  % Number of x grid points        
Ny = floor(Ly/h);  % Number of y grid points       
Ns = floor(T/k);    % Number of samples
N = (Nx+1)*(Ny+1);  % The size of our vectors, and dimensions of square matrices


lambda = 1/sqrt(2);  % Courant number
h = c*k/lambda;      % Grid spacing

% Excitation
point = [0.5,0.5];    % Excitation point
size = 0.5;           % Excitation size

% Creating raised cosine
[X, Y] = meshgrid([0:Nx]*h, [0:Ny]*h);
dist = sqrt((X-point(1)).^2 +(Y-point(2)).^2);
ind = sign(max(-dist+size/2,0)); 
Raised_Cosine = 0.5*ind'.*(1+cos(2*pi*dist'/size));

u1 = Raised_Cosine;
u2 = Raised_Cosine;


u = zeros(Nx+1, Ny+1);

for i=3:Ns-1
    
    u(2:Nx, 2:Ny) = 0.5*(u2(3:Nx+1, 2:Ny) + u2(1:Nx-1, 2:Ny)  ...
                       + u2(2:Nx, 3:Ny+1) + u2(2:Nx, 1:Ny-1)) ...
                       - u1(2:Nx, 2:Ny);
    
    u1 = u2;
    u2 = u;
    
    % Animating the plate
    mesh(u, 'FaceColor', 'interp');
    axis([-10 Nx+10 -10 Ny+20 -1 1])
    drawnow
    axis equal
    
    
end

 
