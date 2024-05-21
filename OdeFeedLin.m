clc;
clear all;

%Import Model
ModelRexRov

%Import Matlab Functions
Cm = @Cmfile;
Dm = @Dmfile;
fd1 = @fd1file;
fd2 = @fd2file;
fd3 = @fd3file;
gm = @gmfile;
Jm = @Jmfile;
Jminv = @Jminvfile;

%Ode Solver
Xinitial = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0].';% x y z phi theta psi, u, v, w, q, p, q, r, f1, f2, f3
[t Xfinal] = ode45(@(t,State) odefun(t,State,M,gm,Jm, Jminv, fd3,fd2,fd1,Dm,Cm,S),[0 10],Xinitial)

function Xdot = odefun(t,State,M,gm,Jm, Jminv, fd3,fd2,fd1,Dm,Cm,S)
Xdot = zeros(15,1);
%Input DEs
x = State(1);
y = State(2);
z = State(3);
phi = State(4);
theta = State(5);
if theta == pi/2
    theta = pi/2 +1e-3;
    State(5) = theta;
end
psi = State(6);
f1 = State(13);
f2 = State(14);
f3 = State(15);
tau = [f1;0;0;0;f2;f3];
Jinv = Jminv(State(1:6));
J = Jm(State(1:6));
V = State(7:12);
u = V(1);
v = V(2);
w = V(3);
p = V(4);
q = V(5);
r = V(6);
g1 = gm(State(1:6));
D1 = Dm(V);
C1 = Cm(V);
Xd = J*V;

Vd = M\(tau -g1 -D1*V-C1*V);
Xdd = J(1:3,1:3)*Vd(1:3) + J(1:3,1:3)*S(V(4:6))*V(1:3); 
Xdot(1:6) = Xd;
Xdot(7:12) = Vd;
ux1 = -0.01*1000*cos(10*t) + 100*(-0.01*100*sin(10*t) - Xdd(1)) + 100*(0.01*10*cos(10*t) - Xd(1)) + 1000*(0.01*sin(10*t) - x);
uy1 = 0 + 1*(0 - Xdd(2)) + 1*(0.01 - Xd(2)) + 1*(0.01*t - y);
uz1 = 0 + 1*(0 - Xdd(3)) + 1*(0.01 - Xd(3)) + 1*(0.01*t - z); 
Rnb = J(1:3,1:3);
Rbn = inv(Rnb);
us = Rbn*[ux1;uy1;uz1];
ux = us(1);
uy = us(2);
uz = us(3);
% u = 0;
% v = 0;
% w = 0;
% p = 0;
% q = 0;
% r = 0;
% theta = 0;
% phi = 0;
% ux = 0;
% uy = 0;
% uz = 0;
Xdot(13) = fd1(f1,f2,f3,p,phi,q,r,theta,u,ux,uy,uz,v,w);
Xdot(14) = fd2(f1,f2,f3,p,phi,q,r,theta,u,ux,uy,uz,v,w);
Xdot(15) = fd3(f1,f2,f3,p,phi,q,r,theta,u,ux,uy,uz,v,w);
end