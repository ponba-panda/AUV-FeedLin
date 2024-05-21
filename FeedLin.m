clc;
clear all;
ModelRexRov;

syms f1 f2 f3 f1d f2d f3d;
syms ux uy uz
tau = [f1;0;0;0;f2;f3];
taud = [f1d;0;0;0;f2d;f3d];
A = C(u,v,w,p,q,r) + D(u,v,w,p,q,r);
Minv = inv(M);
tm1 = Minv*(g(x,y,z,phi,theta,psi) + A*V);

%To compute Vdd
tm2 = (diff(A,u)*(Minv(1,:)*tau) + diff(A,v)*(Minv(2,:)*tau)+ diff(A,w)*(Minv(3,:)*tau) + diff(A,p)*(Minv(4,:)*tau) +  diff(A,q)*(Minv(5,:)*tau) + diff(A,r)*(Minv(6,:)*tau));

T1 = Tnb(1,:);
T2 = Tnb(2,:);
tm3 = diff(g(x,y,z,phi,theta,psi),phi)*(T1*v2) + diff(g(x,y,z,phi,theta,psi),theta)*(T2*v2);
tm4 = (diff(A,u)*(tm1(1,:)) + diff(A,v)*(tm1(2,:))+ diff(A,w)*(tm1(3,:)) + diff(A,p)*(tm1(4,:)) +  diff(A,q)*(tm1(5,:)) + diff(A,r)*(tm1(6,:)));

vdd = Minv*(taud - tm2*V - tm3 + tm4*V + A*tm1 - A*Minv*tau);

a = S(v2)*S(v2)*v1 + S(Minv(4:end,:)*tau)*v1 - S(tm1(4:end,:))*v1 + 2*S(v2)*(Minv(1:3,:)*tau - tm1(1:3,:)) + vdd(1:3,:);

fds = solve(a == [ux;uy;uz],[f1d,f2d,f3d]);

%Convert Everything to a Matlab Function
% fd1 = matlabFunction(fds.f1d,'File','fd1file');
% fd2 = matlabFunction(fds.f2d,'File','fd2file');
% fd3 = matlabFunction(fds.f3d,'File','fd3file');
% afile = matlabFunction(a,'File','afile');
