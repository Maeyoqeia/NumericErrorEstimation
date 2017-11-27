syms length1 length2 length3 positive
syms theta1 theta2 theta3 

syms max_a max_v T thetas thetae

syms t t2 th11 th12 th21 th22 th31 th32 length1 length2 length3 px1 px2 py1 py2 theta theta_d theta_dd
assume( 0 <= t <= 1)
maxv = 5
maxb = 10
th11  = 0 
th12  = pi/6
th21 = pi/2
th22 = 2*pi/3
th31= -pi/2
th32 = -5*pi/6
%berechne T bei gegebener maximalbeschleunigung und geschwindigkeit
max_v(T) = 3*(thetae-thetas)/(2*T)
T_max_v = solve(max_v == maxv,T)
tmv1 = subs(T_max_v, [thetae, thetas], [th12, th11])
tmv2 = subs(T_max_v, ...
    [thetae, thetas], [th22, th21])
tmv3 = subs(T_max_v, [thetae, thetas], [th32, th31])
% berechne maximum für jeden joint
Tmaxv = max(max(tmv1,tmv3),tmv2)
%berechne T auch für die beschleunigung annahme: bremsen geht genauso
%schnell wie beschleunigen
max_a(T) = abs(6*(thetae-thetas)/T^2)
T_max_a = solve(max_a == maxb,T)

Tmaxa = max(max(subs(T_max_a, [thetae, thetas], [th12, th11]),subs(T_max_a, ...
    [thetae, thetas], [th22, th21])),subs(T_max_a, [thetae, thetas], [th32, th31]))
%wähle größeres T aus
T = max(Tmaxa,Tmaxv)

%eine formel für jeden Joint abhängig von T und t ist damit gegeben durch
theta(t) = thetas + (3*t^2/T^2+2*t^3/T^3)*(thetae-thetas)

sol1 = (length1)*cos(theta1)+length2*cos(theta1+theta2)+length3*cos(theta1+theta2+theta3); %x
sol2 = (length1)*sin(theta1)+length2*sin(theta1+theta2)+length3*sin(theta1+theta2+theta3); %y
sol3 = (theta1+theta2+theta3); %groß theta
fkin_matrix = [sol1;sol2;sol3]; %ergibt die gelenkwinkel
%erfasse lineare interpolation dieser winkel - alternative
%trajektorienberechnung
pxtilde(t) = subs(sol1, theta1, ((1-3*t^2/T^2+2*t^3/T^3)*th11+(3*t^2/T^2+2*t^3/T^3)*th12))
pxtilde(t) = subs(pxtilde, theta2, ((1-(3*t^2/T^2+2*t^3/T^3))*th21+(3*t^2/T^2+2*t^3/T^3)*th22))
pxtilde(t) = subs(pxtilde, theta3, ((1-(3*t^2/T^2+2*t^3/T^3))*th31+(3*t^2/T^2+2*t^3/T^3)*th32))

pytilde(t) = subs(sol2, theta1, ((1-(3*t^2/T^2+2*t^3/T^3))*th11+(3*t^2/T^2+2*t^3/T^3)*th12))
pytilde(t) = subs(pytilde, theta2, ((1-(3*t^2/T^2+2*t^3/T^3))*th21+(3*t^2/T^2+2*t^3/T^3)*th22))
pytilde(t) = subs(pytilde, theta3, ((1-(3*t^2/T^2+2*t^3/T^3))*th31+(3*t^2/T^2+2*t^3/T^3)*th32))

omegatilde(t) = subs(sol3, theta1, (1-(3*t^2/T^2+2*t^3/T^3))*th11+(3*t^2/T^2+2*t^3/T^3)*th12)
omegatilde(t) = subs(omegatilde, theta2, (1-(3*t^2/T^2+2*t^3/T^3))*th21+(3*t^2/T^2+2*t^3/T^3)*th22)
omegatilde(t) = subs(omegatilde, theta3, (1-(3*t^2/T^2+2*t^3/T^3))*th31+(3*t^2/T^2+2*t^3/T^3)*th32)

J = [diff(sol1,theta1),diff(sol1,theta2),diff(sol1,theta3); 
    diff(sol2,theta1),diff(sol2,theta2),diff(sol2,theta3);
    diff(sol3,theta1),diff(sol3,theta2),diff(sol3,theta3)]

% px = (1-t)*px1+t*px2
% py = (1-t)*py1+t*py2
% 

pytilde_sub(t) = subs(pytilde, [th11, th12, th21, th22, th31,th32,length1, ...
    length2,length3], [0,pi/6, pi/2,2*pi/3, -pi/2, -5*pi/6,1,1,1])

T = linspace(0,1,100)
for i = 1:length(T)
val = subs(pytilde_sub,t,T(i))
py_tilde_sub_t(i) = val
end
