% h2d.m, solving u_t=Delta u by FFT 
n=16; dx=2*pi/n; dy=dx; x=0:dx:2*pi-dx; y=x; [X Y]=meshgrid(x,y); 
eta=1; h=0.01; mu=zeros(n,n); % holds F-multipliers 
K1=zeros(n,n); % holds Laplacian in F-space
mu2=zeros(n,n); % holds Laplacian in F-space
for j1=1:n; k1=j1-n/2-1; % fill multiplier matrix mu as if fft(u) was centered 
    for j2=1:n; k2=j2-n/2-1;ks=k1*k1+k2*k2; 
        % K1(j1,j2)=ks;
        % K2(j1,j2)=k2;
        mu(j1,j2)=(1-(1-eta)*ks*h)/(1+eta*ks*h);
        % mu2(j1,j2)=(1-(1-eta)*ks*h)/(1+eta*ks*h);
    end
end
% K1=fftshift(K1); 
% K2=fftshift(K2);
mu=fftshift(mu); % now adapt mu to true fft
% mu2=fftshift(mu2); % now adapt mu to true fft
u0=sin(X)+sin(2*Y)+cos(7*exp(-X)); umin=min(min(u0)); umax=max(max(u0));
% print u0
% K1'
% mu'-mu2'
u0'
clf;
surf(X,Y,u0); axis([0 2*pi 0 2*pi umin umax]); grid off;
more=1;u=u0; uf=fft2(u); more=askmore(more);
% uf2=uf;
while more==1 % integration loop
  uf=mu.*uf;  % actual integration, next line for plotting
  % uf2=mu2.*uf2;  % actual integration, next line for plotting
  u=ifft2(uf); 
  % u2=ifft2(uf2);
  surf(X,Y,real(u));  axis([0 2*pi 0 2*pi umin umax]); grid off;
  real(u')
  more=askmore(more);
end

