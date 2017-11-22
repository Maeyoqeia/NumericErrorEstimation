syms length1 length2 length3 positive
syms theta1 theta2 theta3 

syms t positive
syms th11 th12 th21 th22 th31 th32
% length1= 1
% length2 = 1
% length3 = 1

sol1 = (length1)*cos(theta1)+length2*cos(theta1+theta2)+length3*cos(theta1+theta2+theta3); %x
sol2 = (length1)*sin(theta1)+length2*sin(theta1+theta2)+length3*sin(theta1+theta2+theta3); %y
sol3 = (theta1+theta2+theta3); %groß theta
fkin_matrix = [sol1;sol2;sol3]; %ergibt die gelenkwinkel
%erfasse lineare interpolation dieser winkel - alternative
%trajektorienberechnung
pxtilde(t) = subs(sol1, theta1, ((1-t)*th11+t*th12))
pxtilde(t) = subs(pxtilde, theta2, ((1-t)*th21+t*th22))
pxtilde(t) = subs(pxtilde, theta3, ((1-t)*th31+t*th32))

pytilde(t) = subs(sol2, theta1, ((1-t)*th11+t*th12))
pytilde(t) = subs(pytilde, theta2, ((1-t)*th21+t*th22))
pytilde(t) = subs(pytilde, theta3, ((1-t)*th31+t*th32))

omegatilde(t) = subs(sol3, theta1, (1-t)*th11+t*th12)
omegatilde(t) = subs(omegatilde, theta2, (1-t)*th21+t*th22)
omegatilde(t) = subs(omegatilde, theta3, (1-t)*th31+t*th32)

J = [diff(sol1,theta1),diff(sol1,theta2),diff(sol1,theta3); 
    diff(sol2,theta1),diff(sol2,theta2),diff(sol2,theta3);
    diff(sol3,theta1),diff(sol3,theta2),diff(sol3,theta3)]

% px = (1-t)*px1+t*px2
% py = (1-t)*py1+t*py2
% 

pytilde_sub(t) = subs(pytilde, [th11, th12, th21, th22, th31,th32,length1, ...
    length2,length3], [0,pi/6, pi/2,2*pi/3, -pi/2, -5*pi/6,1,1,1])

% T = linspace(0,1,100)
% for i = 1:length(T)
% val = subs(pytilde_sub,t,T(i))
% py_tilde_sub_t(i) = val
% end
