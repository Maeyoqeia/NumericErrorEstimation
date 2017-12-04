testnumber = 2;

error = 0.01;
sSize1 = 2; % gibt an, wieviele Punkte der Linie im kartesischen Raum ausgewählt werden
sSize2 = 40; % gibt an wie viele Zwischenpunkte zwischen jedem Punkt eingefügt werden
angles_ori = [0, pi/2]
angles = sample_multiple(angles_ori,0);
iteration = 1;
startpunkt = [2;1];
endpunkt = [1;1];
startwinkel = 0;
endwinkel = 0;
angles = [startwinkel endwinkel]
dir_path = sprintf('test%d',testnumber);
mkdir(dir_path)
while(1)
t_obj = TransformationObject([2;1],[1;1],sSize1, sSize2, angles,[1,1,1],iteration,testnumber);
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
fprintf(fileID,'%s%f\n%s%d\n%s%d\n%s%d%s%d%s\n%s%d%s%d%s', 'Fehlergröße: ',(error),...
    'Samplesize1: ', (sSize1), 'Samplesize2: ', (sSize2),...
    'Startpunkt: (', startpunkt(1),',',startpunkt(2),')',...
    'Endpunkt: (',endpunkt(1),',',endpunkt(2),')');
fclose(fileID);