close all
for t=[10.45]
    xL=1;
uInf=1;
vInf=1;
omega=4*pi;
gamma=15;
r0=pi/omega;
x0=0.5+xL;
y0=0.5;
xL=3;
xc=@(x) mod(x-uInf*t-xL,1)-x0+xL;
yc=@(y) mod(y-vInf*t,1)-y0;

xx=linspace(xL,1+xL,50);
yy=linspace(0,1,60);
[XX,YY]= meshgrid(xx,yy);
HH=zeros(size(XX));
RR=zeros(size(XX));
for ii=1:size(XX,1)
    for jj=1:size(XX,2)
        r=sqrt(xc(XX(ii,jj))^2+yc(YY(ii,jj))^2);
        RR(ii,jj)=r;
        disp(r)
        if r<r0
            HH(ii,jj) = 10+(gamma/omega)^2*(hfunction(omega*r)-hfunction(pi));
        else
            HH(ii,jj)=10;
        end
    end
end
figure()
surf(XX,YY,HH)
end
        
function h=hfunction(x)
    h=2.*cos(x)+2.*x*sin(x)+1./8.*cos(2*x)+0.25*x*sin(2*x)+0.75*x^2;
end