
testnumber = 1;
arms = [1 1 1 ];
start_pt = [1;2];
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
iteration = 0
singularity = 0
while(1)
    iteration = iteration +1
    [transformed_pts_jsp,singularity] = (transform(sampling, angles,arms)); %to joint space
    if((transformed_pts_jsp(2,1)-transformed_pts_jsp(2,end))> pi)
        singularity = 1
    end
    error_flag = 0;
    newsampling = [];
    exy_array = [];
    l_array_max = [];
    l_array = {};
    px_array = {};
    py_array = {};
    for i = 2:size(transformed_pts_jsp,2)
        
        newsampling = [newsampling sampling(:,i-1)];
        
        j1 = transformed_pts_jsp(:,i-1)
        j2 = transformed_pts_jsp(:,i)
        % th11 = j1(1)
        % th12 = j1(2)
        % th21 =j1(3)
        % th22 = j2(1)
        % th31 = j2(2)
        % th32 = j2(3)
        syms t th11 th12 th21 th22 th31 th32 length1 length2 length3 px1 px2 py1 py2 real
        assume( 0 <= t <= 1)
        px(t) = (1-t)*px1+t*px2;
        py(t) = (1-t)*py1+t*py2;
        px_sub = subs(px, [px1, px2],[sampling(1,i-1), sampling(1,i)]);
        py_sub = subs(py, [py1, py2],[sampling(2,i-1), sampling(2,i)]);
        pytilde(t) = length3*sin(t*th12 + t*th22 + t*th32 - th11*(t - 1) - ...
            th21*(t - 1) - th31*(t - 1)) + length2*sin(t*th12 + t*th22 - ...
            th11*(t - 1) - th21*(t - 1)) + length1*sin(t*th12 - th11*(t - 1));
        
        
        pytilde(t) = length3*sin(t*th12 + t*th22 + t*th32 - th11*(t - 1) - ...
            th21*(t - 1) - th31*(t - 1)) + length2*sin(t*th12 + t*th22 -...
            th11*(t - 1) - th21*(t - 1)) + length1*sin(t*th12 - th11*(t - 1));
        
        pxtilde(t) = length3*cos(t*th12 + t*th22 + t*th32 - th11*(t - 1)...
            - th21*(t - 1) - th31*(t - 1)) + length2*cos(t*th12 + t*th22 -...
            th11*(t - 1) - th21*(t - 1)) + length1*cos(t*th12 - th11*(t - 1));
        
        pxtilde_sub(t) = subs(pxtilde, [th11, th12, th21, th22, th31,th32,length1, ...
            length2,length3], [j1(1), j2(1),j1(2), j2(2),j1(3),j2(3),arms(1),arms(2),arms(3)]);
        pytilde_sub(t) = subs(pytilde, [th11, th12, th21, th22, th31,th32,length1, ...
            length2,length3], [j1(1), j2(1),j1(2), j2(2),j1(3),j2(3),arms(1),arms(2),arms(3)]);
        
       % pxtilde_exp = rewrite(pxtilde_sub, 'exp');
      %  pytilde_exp = rewrite(pytilde_sub, 'exp');
        
        pxtilde_sub = real(pxtilde_exp);
        pytilde_sub = real(pytilde_exp);
        px_array{end+1} = (real(pxtilde_sub));
        py_array{end +1} = (real(pytilde_sub));
        max_t_y = real(solve(diff(py_sub-pytilde_sub)==0,t)); %das t wo y die maximale abweichung hat
        max_t_x = real(solve(diff(px_sub- pxtilde_sub)==0,t));
        %vpasolve 
        max_y = double(pytilde_sub(max_t_y)) ;%die maximale Abweichung an max_t_y
        max_x = double(pxtilde_sub(max_t_x));
        if isempty(max_y )
            string = 'empty';
            if(pytilde_sub(0) < pytilde_sub(1))
                max_y = pytilde_sub(1);
                max_t_y = 1;
            else
                max_y = pytilde_sub(0);
                max_t_y = 0;
            end
        end
        if isempty(max_x)
            if(pxtilde_sub(0) < pxtilde_sub(1))
                max_x = pxtilde_sub(1);
                max_t_x = 1;
            else
                max_x = pxtilde_sub(0);
                max_t_x = 0;
            end
        end
        %möglichkeit 1
        %schätze den fehler nach unten ab durch l<=sqrt(Ex^2+Ey^2)
        punkty = (1-max_t_y)*sampling(2,i-1)+max_t_y*sampling(2,i);
        punktx = (1-max_t_x)*sampling(1,i-1)+max_t_x*sampling(1,i);
        
        Ex = abs(max_x - punktx);
        Ey = abs(max_y-punkty);
        exy = (Ex^2+Ey^2);
        exy_array(end+1) = double(exy);
        if(exy > error_size)
            error_flag = 1
            newsampling = [newsampling sampling(:,i-1)*(0.5)+sampling(:,i)*0.5];
        end
        
        
        %möglichkeit2 berechne den maximaltwert des vektors
        
        l(t) = ((px_sub-pxtilde_sub)^2+(py_sub-pytilde_sub)^2);
        t_err = solve(diff(l,t)==0,t);
        max_err = double(real(l(t_err)));
        l_array_max = [l_array_max double(max_err)];
        l_array{end+1} = real(l);
        if(max_err > error_size)
            %     error_flag = 1
            %     newsampling = [newsampling sampling(:,i-1)*(0.5)+sampling(:,i)*0.5]
        end
        
    end
    
    newsampling = [newsampling sampling(:,i)];
    %plotting
    %originale linie
    line = sample(start_pt,end_pt,-1,sSize2); %sample points from line
    %markiere alle punkte vom newsampling
    
    figure%('visible','off');
    plot(line(1,:),line(2,:))
    hold on;
    plot(newsampling(1,:),newsampling(2,:), 'o')
    hold on;
    for p_index = 1:size(px_array,2)
        px_el = px_array(p_index);
        py_el = py_array(p_index);
        T = linspace(0,1,sSize2);
        for i = 1:length(T)
            val = subs(px_el,t,T(i));
            pxtilde_sub_i(i) = val;
        end
        for i = 1:length(T)
            val = subs(py_el,t,T(i));
            pytilde_sub_i(i) = val;
        end
        a = [pxtilde_sub_i;pytilde_sub_i]
        plot(a(1,:),a(2,:)); %plotte transformierte Punkte
        hold on;
    end
    title('Originale und transformierte Punkte')
    if(singularity == 1)
        title('Originale und transformierte Punkte, SINGULARITÄT')
    end
    str = sprintf('test_%d/iteration%d.jpg',testnumber,iteration);
    saveas(gcf,str)
    
    figure%('visible','off'); %plot errors
    
    %plotte l - den genauen Fehler
    for el= 1:size(l_array,2)
        l_curr = l_array(el);
        T = linspace(0,1,sSize2);
        for i = 1:length(T)
            val = subs(l_curr,t,T(i));
            l_curr_i(i) = double(val);
        end
        T_larray = linspace(el,el+1,sSize2)
        plot(T_larray, l_curr_i);
        hold on;
        plot([el, el+1], [exy_array(el), exy_array(el)])
        hold on;
    end
    difference = abs(l_array_max-exy_array);
    max_diff = max(difference);
    str_tit = sprintf('Quadr. Fehler und Abschätzung, max. Differenz = %f',max_diff);
    title(str_tit)
    str = sprintf('test_%d/fehler%d.jpg',testnumber,iteration);
    saveas(gcf,str)
    sampling = double(newsampling);
    angles = zeros(1,size(sampling,2));
    if(error_flag == 0) %nie ein fehler aufgetreten
        break
    end
    
end

file_name = 'description.txt';
full_name = fullfile(dir_path, file_name);
fileID = fopen(full_name,'wt');
fprintf(fileID,'%s%f\n%s%d\n%s%d\n%s%d%s%d%s\n%s%d%s%d%s', 'Fehlergröße: ',(error),...
    'Samplesize1: ', (sSize1), 'Samplesize2: ', (sSize2),...
    'Startpunkt: (', start_pt(1),',',start_pt(2),')',...
    'Endpunkt: (',end_pt(1),',',end_pt(2),')');fclose(fileID);
