
arms = [1 1 1 ];
testnumber = 1;
start_pt = [2;1];
end_pt = [1;1];
error_size = 0.01;
sSize1 = 2; % gibt an, wieviele Punkte der Linie im kartesischen Raum ausgewählt werden
sSize2 = 40; % gibt an wie viele Zwischenpunkte zwischen jedem Punkt eingefügt werden
angles_ori = zeros(1,sSize1);
angles = sample_multiple(angles_ori,0);
iteration = 1;

dir_path = sprintf('test__%d',testnumber);
mkdir(dir_path)
sampling = sample(start_pt,end_pt,-1,sSize1); %sample points from line
while(1)
transformed_pts_jsp = (transform(sampling, angles,arms)) %to joint space
    error_flag = 0;
    newsampling = [];
for i = 2:size(transformed_pts_jsp,2)

    newsampling = [newsampling sampling(:,i-1)]
    
    j1 = transformed_pts_jsp(:,i-1)
    j2 = transformed_pts_jsp(:,i)
syms length1 length2 length3 positive
syms theta1 theta2 theta3 

syms Ta

Ta = vmax/amax
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


pxtilde_sub(t) = subs(pxtilde, [th11, th12, th21, th22, th31,th32,length1, ...
    length2,length3], [j1(1), j2(1),j1(2), j2(2),j1(3),j2(3),arms(1),arms(2),arms(3)]);
pytilde_sub(t) = subs(pytilde, [th11, th12, th21, th22, th31,th32,length1, ...
    length2,length3], [j1(1), j2(1),j1(2), j2(2),j1(3),j2(3),arms(1),arms(2),arms(3)]);
%möglichkeit 1
max_t_y = solve(diff(py_sub-pytilde_sub)==0,t) %das t wo y die maximale abweichung hat
max_t_x = solve(diff(px_sub- pxtilde_sub)==0,t)

max_y = double(pytilde_sub(max_t_y)) %die maximale Abweichung an max_t_y
max_x = double(pxtilde_sub(max_t_x)) 
if isempty(max_y ) 
    if(pytilde_sub(0) < pytilde_sub(1))
        max_y = pytilde_sub(1)
        max_t_y = 1
    else
        max_y = pytilde_sub(0)
        max_t_y = 0
    end
end
if isempty(max_x)
    if(pxtilde_sub(0) < pxtilde_sub(1))
        max_x = pxtilde_sub(1)
        max_t_x = 1
    else
        max_x = pxtilde_sub(0)
        max_t_x = 0
    end
end
%möglichkeit 1
%schätze den fehler nach unten ab durch l<=sqrt(Ex^2+Ey^2)
punkty = (1-max_t_y)*sampling(2,i-1)+max_t_y*sampling(2,i)
punktx = (1-max_t_x)*sampling(1,i-1)+max_t_x*sampling(1,i)

Ex = abs(max_x - punktx)
Ey = abs(max_y-punkty)
exy = (Ex^2+Ey^2)
poss1 = [poss1 double(exy)]
if(exy > error_size)
    error_flag = 1
    newsampling = [newsampling sampling(:,i-1)*(0.5)+sampling(:,i)*0.5]
 end

%möglichkeit2 berechne den maximaltwert des vektors

l(t) = ((px_sub-pxtilde_sub)^2+(py_sub-pytilde_sub)^2)
t_err = solve(diff(l)==0,t)
max_err = l(t_err)
poss2 = [poss2,double(max_err)]
if(max_err > error_size)
%     error_flag = 1
%     newsampling = [newsampling sampling(:,i-1)*(0.5)+sampling(:,i)*0.5]
end

end

newsampling = [newsampling sampling(:,i)]
sampling = double(newsampling)
angles = zeros(1,size(sampling,2))
if(error_flag == 0) %nie ein fehler aufgetreten
    break
end

end
file_name = 'description.txt';
full_name = fullfile(dir_path, file_name);
fileID = fopen(full_name,'wt');
fprintf(fileID,'%s%f\n%s%d\n%s%d', 'Fehlergröße ',(error_size), 'Samplesize1 ', (sSize1), 'Samplesize2 ', (sSize2));
fclose(fileID);
