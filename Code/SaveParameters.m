clear('prm'); %clear varaible if exists in workspace as not to add to it
%% Length Parameters
prm.r=1;
prm.L=2;
prm.H=3;

%% For Forward Kinematics
syms d1 d2 d3 x y theta L r H real
Eq1=(x-d1)^2+y^2==L^2;
Eq2=(x+r*cos(theta)-d2)^2+(y+r*sin(theta))^2==L^2;
Eq3=(x+r*cos(theta+pi/3)-d3)^2+(y+r*sin(theta+pi/3)-H)^2==L^2;

Eq12=collect(Eq1-Eq2,[x,y]);
Eq13=collect(Eq1-Eq3,[x,y]);

[A,b] = equationsToMatrix([Eq12;Eq13],[x;y]);

X=linsolve(A,b);
xt=expand(X(1)); %Expand breaks down cos(theta+pi/3) -> sin(theta)*cos(pi/3)+cos(theta)*sin(pi/3)
yt=expand(X(2));

syms T; %T=tan(theta/2)
sinT=2*T/(1+T^2);
cosT=(1-T^2)/(1+T^2);

xT=subs(xt,[sin(theta),cos(theta)],[sinT,cosT]);
yT=subs(yt,[sin(theta),cos(theta)],[sinT,cosT]);

prm.Eq1T=subs(subs(Eq1-L^2,[x,y],[xT,yT]),[r,L,H],[prm.r,prm.L,prm.H]);
prm.xT=subs(xT,[r,L,H],[prm.r,prm.L,prm.H]);
prm.yT=subs(yT,[r,L,H],[prm.r,prm.L,prm.H]);
%% Save to file
save('Parameters','prm');