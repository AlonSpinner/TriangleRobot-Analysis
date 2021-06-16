function Bool=checkSolutionValidity(x,y,theta,d1,d2,d3)
load('Parameters.mat');
r=prm.r;
H=prm.H;
L=prm.L;

D1=[d1,0];
D2=[d2,0];
D3=[d3,H];
A=[x,y];
B=A+r*[cos(theta),sin(theta)];
C=A+r*[cos(theta+pi/3),sin(theta+pi/3)];

dL1=abs(vecnorm(A-D1)-L);
dL2=abs(vecnorm(B-D2)-L);
dL3=abs(vecnorm(C-D3)-L);
dL=[dL1,dL2,dL3];

Bool=~any(dL>1e-10);
end