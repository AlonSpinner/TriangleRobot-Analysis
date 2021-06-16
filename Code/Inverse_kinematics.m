function [d1,d2,d3]=Inverse_kinematics(x,y,theta,elbows,varargin)
method='numerical';
for ind	= 1:2:length(varargin)
	comm	= lower(varargin{ind});
	switch comm
		case 'method'
			method	= varargin{ind+1};
	end
end

switch method
    case 'numerical'
        load('Parameters.mat');
        r=prm.r;
        L=prm.L;
        H=prm.H;
    case 'symbolic'
        syms r L H real
end

d1=x+elbows(1).*sqrt(L^2-y.^2);
d2=x+r.*cos(theta)+elbows(2).*sqrt(L^2-(y+r.*sin(theta)).^2);
d3=x+r.*cos(theta+pi/3)+elbows(3).*sqrt(L.^2-(y+r.*sin(theta+pi/3)-H).^2);
end