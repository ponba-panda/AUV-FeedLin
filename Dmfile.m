function D = Dmfile(in1)
%Dmfile
%    D = Dmfile(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    28-Apr-2024 18:46:55

u = in1(1,:);
v = in1(2,:);
w = in1(3,:);
p = in1(4,:);
q = in1(5,:);
r = in1(6,:);
D = reshape([abs(u).*7.4822e+2+7.482e+1,0.0,0.0,0.0,0.0,0.0,0.0,abs(v).*9.9253e+2+6.948e+1,0.0,0.0,0.0,0.0,0.0,0.0,abs(w).*1.82101e+3+7.284e+2,0.0,0.0,0.0,0.0,0.0,0.0,abs(p).*6.72e+2+2.688e+2,0.0,0.0,0.0,0.0,0.0,0.0,abs(q).*7.7444e+2+3.0977e+2,0.0,0.0,0.0,0.0,0.0,0.0,abs(r).*5.2327e+2+1.05e+2],[6,6]);
end
