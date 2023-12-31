% ac2d.m, solving u_t=Delta u+u-u^3 by FFT 
n=16; dx=2*pi/n; dy=dx; x=0:dx:2*pi-dx; y=x; [X Y]=meshgrid(x,y); 
nu=1; h=0.01; mu=ones(n,n); % holds F-multipliers 
for j1=1:n; k1=j1-n/2-1; % fill multiplier matrix mu as if fft(u) was centered 
    for j2=1:n; k2=j2-n/2-1;ks=k1*k1+k2*k2; mu(j1,j2)=1/(1+nu*ks*h);end
end
mu=fftshift(mu); % now adapt mu to true fft
u0=zeros(n,n); nb=8; % blockwise random IC
% for j1=0:n/nb-1 
%     for j2=0:n/nb-1
%    u0(j1*nb+1:(j1+1)*nb, j2*nb+1:(j2+1)*nb)=0.5*(rand(1,1)-0.5)*ones(nb,nb); 
%     end
% end
u0=sin(X)+sin(2*Y)+cos(7*exp(-X));
u0'
%u0=0.5*(rand(n,n)-0.5);
% figure(1);clf; colormap cool; p1=surf(X,Y,u0);shading interp;
% axis([0 2*pi 0 2*pi -2 2]);view(-10,70);grid off; 
% tits=['t=0 '];title(tits);
more=1;u=u0; t=0; uf=fft2(u); more=askmore(more);
while more==1 % integration loop
  uf=mu.*(uf+h*fft2(u-u.^3));  u=ifft2(uf); t=t+h;
  %clf; pcolor(x,y,u);shading interp; colorbar;
  real(u')
  % clf; p1=surf(X,Y,real(u));
  % colormap cool; shading interp;
  % axis([0 2*pi 0 2*pi -1.1 1.1]);view(-10,70);grid off; 
  % tits=['t=', num2str(t)];title(tits);
  more=askmore(more);
end

