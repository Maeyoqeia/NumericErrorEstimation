
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

dir_path = sprintf('test_%d',testnumber);
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
% th11 = j1(1)
% th12 = j1(2)
% th21 =j1(3)
% th22 = j2(1)
% th31 = j2(2)
% th32 = j2(3)
syms t th11 th12 th21 th22 th31 th32 length1 length2 length3 px1 px2 py1 py2
assume( 0 <= t <= 1)

pytilde(t) = length3*sin(t*th12 + t*th22 + t*th32 - th11*(t - 1) - th21*(t - 1) - th31*(t - 1)) + length2*sin(t*th12 + t*th22 - th11*(t - 1) - th21*(t - 1)) + length1*sin(t*th12 - th11*(t - 1));

pxtilde(t) = length3*cos(t*th12 + t*th22 + t*th32 - th11*(t - 1)...
    - th21*(t - 1) - th31*(t - 1)) + length2*cos(t*th12 + t*th22 -...
    th11*(t - 1) - th21*(t - 1)) + length1*cos(t*th12 - th11*(t - 1));

pxtilde_sub(t) = subs(pxtilde, [th11, th12, th21, th22, th31,th32,length1, ...
    length2,length3], [j1(1), j2(1),j1(2), j2(2),j1(3),j2(3),arms(1),arms(2),arms(3)]);
pytilde_sub(t) = subs(pytilde, [th11, th12, th21, th22, th31,th32,length1, ...
    length2,length3], [j1(1), j2(1),j1(2), j2(2),j1(3),j2(3),arms(1),arms(2),arms(3)]);

max_t_y = solve(diff(pytilde_sub)==0,t) %das t wo y die maximale abweichung hat
max_t_x = solve(diff(pxtilde_sub)==0,t)

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
exy = sqrt(Ex^2+Ey^2)
if(exy > error_size)
    error_flag = 1
    newsampling = [newsampling sampling(:,i-1)*(0.5)+sampling(:,i)*0.5]
end
px(t) = (1-t)*px1+t*px2
py(t) = (1-t)*py1+t*py2
px_sub = subs(px, [px1, px2],[sampling(1,i-1), sampling(1,i)])
py_sub = subs(py, [py1, py2],[sampling(2,i-1), sampling(2,i)])

%möglichkeit2 berechne den maximaltwert des vektors

l(t) = (px_sub-pxtilde_sub)^2+(py_sub-pytilde_sub)^2
t_err = solve(diff(l)==0,t)
max_err = l(t_err)
% 
% %erste möglichkeit. interpoliere zwischen den ts, je nach höhe der
% %abweichung
% if isempty(max_y ) 
%     if(pytilde_sub(0) < pytilde_sub(1))
%         max_y = pytilde_sub(1)
%         max_t_y = 1
%     else
%         max_y = pytilde_sub(0)
%         max_t_y = 0
%     end
% end
% if isempty(max_x)
%     if(pxtilde_sub(0) < pxtilde_sub(1))
%         max_x = pxtilde_sub(1)
%         max_t_x = 1
%     else
%         max_x = pxtilde_sub(0)
%         max_t_x = 0
%     end
% end
% %berechne angestrebte punkte auf der geraden Linie genau an den t's
% %also zwischen sampling(:,i-1) und sampling(:,i)
% punkty = (1-max_t_y)*sampling(2,i-1)+max_t_y*sampling(2,i)
% punktx = (1-max_t_x)*sampling(1,i-1)+max_t_x*sampling(1,i)
% 
% error_y_max = (punkty-max_y)^2
% error_x_max = (punktx-max_x)^2
% error = error_y_max+error_x_max
% 
% max_error = max(error_y_max,error_x_max)
% error_y_weight = error_y_max/(error_y_max+error_x_max)
% t_error_weight = max_t_x*(1-error_y_weight)+max_t_y*error_y_weight
% 
% 
% %möglichkeit 2 wähle die achse mit der größeren abweichung
% max_tt = max(max_y,max_x) %%die größere maximale abweichung beider achsen x,y
%     if max_tt == max_y
%         fact_t = max_t_y   %das t, was zu der größeren max. Abweichung gehört
%     else
%         fact_t = max_t_x
%     end
% %t_error_weight = max_tt
%         %setze neuen punkt an t_error_weight
% s_point = sampling(:,i-1)*(1-t_error_weight)+sampling(:,i)*t_error_weight
% s_point_error = sum((s_point-[pxtilde_sub(t_error_weight);pytilde_sub(t_error_weight)]).^2)
% %teste ob der fehler zu groß ist
% if(error > error_size)
%     %fuege samplingpunkt hinzu 
%     error_flag = 1
%     
% end
%    newsampling = [newsampling s_point]
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
