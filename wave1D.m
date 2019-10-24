clear all
close all

Fs=44100;  %Sample Rate
c=100;     %Wave Speed
h=1/400;   %Grid Spacing
L=1;       %Length of string
N=round(L/h);  %Number of grid segments
T=0.1;       %Length of simulation

u=zeros(round(T*Fs),N+1);  %Create matrix to store shape of string at each timestep 

u(1,:)=[hann(N+1)]'; %Give string some shape at first timestep

u(2,:)=u(1,:);        %Since 0 velocity, give timestep 1 same shape.

u(1,1)=0;
u(1,N+1)=0;
u(2,1)=0;
u(2,N+1)=0; %Fixed Boundary Conditions


for i=3:round(T*Fs)
   for j=2:N
    
      d=(((1/Fs)^2)*(c^2))/(h^2);   
      u(i,j)=2*(1-d)*u(i-1,j)+d*(u(i-1,j-1)+u(i-1,j+1))-u(i-2,j);
      u(i,1)=0;
      u(i,N+1)=0;
      
   end 
end

points=[0:N];

for i=1:round(T*Fs)
    plot(points,u(i,:));
    axis([0 N -1 1]);
    M(i)=getframe();
    
end

close all
