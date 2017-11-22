
arms = [1 1 1 ]
testnumber = 1;
start_pt = [2;1];
end_pt = [1;1];
error = 0.1;
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
    error_flag = 0
    newsampling = []
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
syms t th11 th12 th21 th22 th31 th32 length1 length2 length3
assume( 0 <= t <= 1)

pytilde(t) = length3*sin(t*th12 + t*th22 + t*th32 - th11*(t - 1) - th21*(t - 1) - th31*(t - 1)) + length2*sin(t*th12 + t*th22 - th11*(t - 1) - th21*(t - 1)) + length1*sin(t*th12 - th11*(t - 1))

pxtilde(t) = length3*cos(t*th12 + t*th22 + t*th32 - th11*(t - 1)...
    - th21*(t - 1) - th31*(t - 1)) + length2*cos(t*th12 + t*th22 -...
    th11*(t - 1) - th21*(t - 1)) + length1*cos(t*th12 - th11*(t - 1))

pxtilde_sub(t) = subs(pxtilde, [th11, th12, th21, th22, th31,th32,length1, ...
    length2,length3], [j1(1), j2(1),j1(2), j2(2),j1(3),j2(3),arms(1),arms(2),arms(3)])
pytilde_sub(t) = subs(pytilde, [th11, th12, th21, th22, th31,th32,length1, ...
    length2,length3], [j1(1), j2(1),j1(2), j2(2),j1(3),j2(3),arms(1),arms(2),arms(3)])

max_t_y = solve(diff(pytilde_sub)==0,t)
max_t_x = solve(diff(pxtilde_sub)==0,t)

max_y = double(pytilde_sub(max_t_y))
max_x = double(pxtilde_sub(max_t_x)) 

max_tt = max(max_y,max_x)
    if max_tt == max_y
        fact_t = max_t_y
    else
        fact_t = max_t_x
%berechne fehler
punkt = (1-fact_t)*sampling(i-1)
if isempty(max_y ) 
    max_y = 0
end
if isempty(max_x)
    max_x = 0
end
if(max_y > error|| max_x > error)
    %fuege samplingpunkt hinzu 
    error_flag = 1
    
    end
    newsampling = [newsampling sampling(:,i-1).*(1-fact_t)+sampling(:,i).*fact_t]
end

end

if(error_flag == 0) %nie ein fehler aufgetreten
    break
end
newsampling = [newsampling sampling(:,i)]
sampling = double(newsampling)
angles = zeros(1,size(sampling,2))
end
file_name = 'description.txt';
full_name = fullfile(dir_path, file_name);
fileID = fopen(full_name,'wt');
fprintf(fileID,'%s%f\n%s%d\n%s%d', 'Fehlergröße ',(error), 'Samplesize1 ', (sSize1), 'Samplesize2 ', (sSize2));
fclose(fileID);
