function J = Jmfile(in1)
%Jmfile
%    J = Jmfile(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    28-Apr-2024 18:46:54

x = in1(1,:);
y = in1(2,:);
z = in1(3,:);
phi = in1(4,:);
theta = in1(5,:);
psi = in1(6,:);
t2 = cos(phi);
t3 = cos(psi);
t4 = cos(theta);
t5 = sin(phi);
t6 = sin(psi);
t7 = sin(theta);
t8 = tan(theta);
t9 = 1.0./t4;
J = reshape([t3.*t4,t4.*t6,t7.*-1.0,0.0,0.0,0.0,t2.*t6.*-1.0+t3.*t5.*t7,t2.*t3+t5.*t6.*t7,t4.*t5,0.0,0.0,0.0,t5.*t6+t2.*t3.*t7,t3.*t5.*-1.0+t2.*t6.*t7,t2.*t4,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,t5.*t8,t2,t5.*t9,0.0,0.0,0.0,t2.*t8,t5.*-1.0,t2.*t9],[6,6]);
end