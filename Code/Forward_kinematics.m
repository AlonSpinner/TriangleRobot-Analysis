function [xVec,yVec,thetaVec]=Forward_kinematics(d1num,d2num,d3num)
load('Parameters.mat');
syms d1 d2 d3 T real

Eq1T=prm.Eq1T;
xT=prm.xT;
yT=prm.yT;
Eq1Tnum=subs(Eq1T,[d1,d2,d3],[d1num,d2num,d3num]);

ExprT=collect(lhs(Eq1Tnum),T);
[polyT,~]=numden(ExprT);
coffs=double(coeffs(polyT,'all'));
TVec=roots(coffs);

sVec = 2*TVec./(1+TVec.^2);
cVec = (1-TVec.^2)./(1+TVec.^2);
thetaVec = double(atan2(sVec,cVec));

xTNum=subs(xT,[d1,d2,d3],[d1num,d2num,d3num]);
yTNum=subs(yT,[d1,d2,d3],[d1num,d2num,d3num]);
xVec=double(subs(xTNum,T,TVec));
yVec=double(subs(yTNum,T,TVec));
end