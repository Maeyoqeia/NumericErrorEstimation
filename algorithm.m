testnumber = 4;

error = 0.001;
sSize1 = 2; % gibt an, wieviele Punkte der Linie im kartesischen Raum ausgewählt werden
sSize2 = 40; % gibt an wie viele Zwischenpunkte zwischen jedem Punkt eingefügt werden
angles_ori = zeros(1,sSize1);
angles = sample_multiple(angles_ori,0);
iteration = 1;

dir_path = sprintf('test%d',testnumber);
mkdir(dir_path)
while(1)
t_obj = TransformationObject([0;0],[0;1],sSize1, sSize2, angles,[1,1,1],iteration,testnumber);
t_obj.toJSpace()
t_obj.trajGen()
t_obj.toTSpace()
t_obj.computeError();
t_obj.plot()
if(t_obj.maxError <= error)
    break
end
angles = sample_multiple(angles, 1);
sSize1 = size(angles,2);
iteration = iteration+1;
end
file_name = 'description.txt';
full_name = fullfile(dir_path, file_name);
fileID = fopen(full_name,'wt');
fprintf(fileID,'%s%f\n%s%d\n%s%d', 'Fehlergröße ',(error), 'Samplesize1 ', (sSize1), 'Samplesize2 ', (sSize2));
fclose(fileID);