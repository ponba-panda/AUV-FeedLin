function Jinv = Jminvfile(in1)
%Jminvfile
%    Jinv = Jminvfile(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    28-Apr-2024 18:46:55

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
t8 = t2.^2;
t9 = t3.^2;
t10 = t4.^2;
t11 = t5.^2;
t12 = t6.^2;
t13 = t7.^2;
t14 = t8.*t10;
t15 = t9.*t10;
t16 = t8+t11;
t17 = t8.*t13;
t18 = t10.*t11;
t19 = t9.*t13;
t20 = t10.*t12;
t21 = t11.*t13;
t22 = t12.*t13;
t23 = 1.0./t16;
t24 = t9.*t17;
t25 = t12.*t14;
t26 = t11.*t15;
t27 = t12.*t17;
t28 = t11.*t19;
t29 = t12.*t18;
t30 = t12.*t21;
t31 = t9.*t14;
t32 = t14+t17+t18+t21;
t33 = t15+t19+t20+t22;
t34 = 1.0./t32;
t35 = 1.0./t33;
t36 = t24+t25+t26+t27+t28+t29+t30+t31;
t37 = 1.0./t36;
Jinv = reshape([t3.*t4.*t35.*1.0,t37.*(t3.*t5.*t7.*-1.0+t2.*t6.*t10+t2.*t6.*t13).*-1.0,t37.*(t2.*t3.*t7+t5.*t6.*t10+t5.*t6.*t13).*1.0,0.0,0.0,0.0,t4.*t6.*t35.*1.0,t37.*(t2.*t3.*t10+t2.*t3.*t13+t5.*t6.*t7).*1.0,t37.*(-t2.*t6.*t7+t3.*t5.*t10.*1.0+t3.*t5.*t13.*1.0).*-1.0,0.0,0.0,0.0,(t7.*-1.0)./(t10+t13),t4.*t5.*t34.*1.0,t2.*t4.*t34.*1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,t2.*t23.*1.0,t5.*t23.*-1.0,0.0,0.0,0.0,t4.*tan(theta).*-1.0,t4.*t5.*t23.*1.0,t2.*t4.*t23.*1.0],[6,6]);
end
