%------------------Newton Method------------------%
function []=NewtonMethod(a, b)

%{
N=100;

x = linspace(-5,5);
y = linspace(5,-50);
[X,Y] = meshgrid(x,y);

F=10*(Y-X.^2).^2+(1-X).^2;
contour(X,Y,F)
%}
ff=@(x,y) 10*(y-x^2)^2+(1-x)^2;
gf=@(x,y) [-40*x*(y-x^2)-2*(1-x);20*(y-x^2)]; 
Hf=@(x,y) [80*y^2-40*(y-x^2)+2 -40*x;-40*x 20];
p = ((-inv(Hf(0,-1)))*gf(0,-1))
lamda = (-Hf(0,-1))-(gf(0,-1)/p)
norm(p)
norm([2/42;20/20])
i = inv(Hf(0,-1))
i*(-gf(0,-1))
%{
F =@(x,y) ff(x,y) + gf(x,y)*p + .5*Hf(x,y)*p

x = linspace(-5,5);
y = linspace(5,-50);
[X,Y] = meshgrid(x,y);

contour(X,Y,F(X,Y))


x0=[a;b];
%hold on;
%plot(x0(1),x0(2),'ok','MarkerFaceColor','k','MarkerSize',4)
w=0.8;
c=0.1;
for j=1:17
    pk=-(inv(Hf(x0(1),x0(2)))*gf(x0(1),x0(2)));
    alpha=1;
    while x0(2)+alpha*pk(2)<=0
        alpha=alpha*w;
    end
    nf=ff(x0(1)+alpha*pk(1),x0(2)+alpha*pk(2));
    dd=gf(x0(1),x0(2))'*pk;
    while (nf>ff(x0(1),x0(2))+c*dd*alpha)
        alpha=alpha*w;
        nf=ff(x0(1)+alpha*pk(1),x0(2)+alpha*pk(2));
    end
    Hf(x0(1),x0(2))
    alpha
    x1=x0+alpha*pk
    %plot(x1(1),x1(2),'ok','MarkerFaceColor','k','MarkerSize',4)
    %plot([x0(1);x1(1)],[x0(2);x1(2)],'k','LineWidth',2)
    x0=x1;
    pause
end
%}
%---------------------------------------------------%