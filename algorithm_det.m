ssize = 40
error = 0.1;
start_pt = [0;0];
end_pt = [0;1];
angle = 0
arms = [1,1,1];
%init

thetas = ikin(start_pt(1), start_pt(2), angle, arms);
maxd = arms(1)*arms(2);
nue = 0.5; %lernrate
dist = end_pt-start_pt;
dist_len = norm(dist);
step = 0;
points = [];
i = 1;
theta_list = zeros(3,1)
while(1)
theta_list(:,i) = thetas'
[j_0,d_0] = jac(arms(1),arms(2),arms(3),thetas(1),thetas(2),thetas(3));
t_i = abs(d_0)/max_d*nue + t_i_old

current_pt = start_pt + step*dist
points(1,i) = current_pt(1);
points(2,i) = current_pt(2);
points
i = i+1
thetas = ikin (current_pt(1),current_pt(2), angle, arms)

if(step >= 1)
    break
end
end

%berechne trajektorie
traj = trajectory1(theta_list, ssize)
path = retransform(arms, traj) %endgültiger pfad, gesampelt mit 40*lenght(traj)
%annahme: nur zwei gegebene punkte der geraden
ssize2 = length(theta_list)+(length(theta_list)-1)*ssize-2
  ori_pts = [];
  ori = [start_pt(1) end_pt(1)]
  ori2 = [start_pt(2) end_pt(2)]
  ori_pts(1,:) = sample_multiple(ori, ssize2);
  ori_pts(2,:) = sample_multiple(ori2,ssize2);

 trs_pts = path;
            trs_pts(3,:) = [];
            abstand = sqrt(sum((trs_pts - ori_pts).^2));
            abstand_quadriert = sum((trs_pts-ori_pts).^2);
            error = sum((trs_pts-ori_pts).^2);
            meanError = sum(error)/size(error,2);
            maxError = max(max(error));

            
            figure
            ax1 = subplot(2,1,1);
            l = length(error);
            five_p = ceil(l*0.05);
            [sorted_values,sortIndex] = sort(error(:),'descend');  
            max_index = sortIndex(1:five_p);  %# Get a linear index into A of the 5 largest values
        
            h(1) = plot(error);
            %line(xlim,[min_val,min_val])
            hold on;
            h(2) = plot(max_index, error(max_index), 'd');
            str = strcat('Fehler (Mean =',num2str(double(meanError)),')'); %double precision
            xlabel('Samplingpunkte')
            ylabel('Fehlergröße')
            legend(h, 'Fehler', '5% der Maximalwerte','Location','southeastoutside')
            title(ax1,str);

            ax2 = subplot(2,1,2);
            plot(trs_pts(1,:),trs_pts(2,:));
            hold on;
            plot(ori_pts(1,:),ori_pts(2,:));
            hold on;
            j(1) = plot(start_pt(1),start_pt(2), 'ob');
            hold on;
            j(2) = plot(end_pt(1),end_pt(2), 'og');
            xlabel('Weg in x-Richtung')
            ylabel('Weg in y-Richtung')
            legend(j, 'Startpunkt', 'Endpunkt','Location','southeastoutside')
            title(ax2,'Originale und transformierte Punkte');
            hold off;
           % str = sprintf('test%d/iteration%d.jpg',testnumber,iteration);
           % saveas(gcf,str)

