classdef TransformationObject < handle
    %TRANSFORMATIONOBJECT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sampleSize_1
        sampleSize_2
        start_pt
        end_pt
        angles
        arms
        original_pts
        transformed_pts_tsp
        transformed_pts_jsp
        error
        trajectory
        meanError
        maxError
        original_pts_polar_phi
        original_pts_polar_r
        dets2
        kond2
        iteration
        testnumber
        def
        singularity
    end
    
    methods
        function obj = TransformationObject(start_pt,end_pt, sampleSize_1, sampleSize_2,angles,arms,iteration,testnumber)
            obj.start_pt = start_pt;
            obj.end_pt = end_pt;
            obj.sampleSize_1 = sampleSize_1;
            obj.sampleSize_2 = sampleSize_2;
            obj.angles = angles;
            obj.arms = arms;
            obj.meanError = realmax;
            obj.maxError = realmax;
            obj.iteration = iteration;
            obj.testnumber = testnumber;
        end
        function obj = toJSpace(obj)
            obj.original_pts = sample(obj.start_pt,obj.end_pt,-1,obj.sampleSize_1); %sample points from line
            % [obj.original_pts_polar_phi,obj.original_pts_polar_r] = cart2pol(obj.original_pts(1,:),(obj.original_pts(2,:)));
            [obj.transformed_pts_jsp,singularity] = (transform(obj.original_pts, obj.angles,obj.arms)); %to joint space
            if((obj.transformed_pts_jsp(2,1)-obj.transformed_pts_jsp(2,end))> pi)
                obj.singularity = 1;
            end
            %berechne determinante von jacobi matrix mit werten(winkeln)
            %aus den spalten von pts_jsp
            
        end
        function obj = trajGen(obj)
            %build trajectory from transformed points
            obj.trajectory = (trajectory1(obj.transformed_pts_jsp,obj.sampleSize_2)); %in joint space
            for p = 1:size(obj.trajectory,2)
                point = obj.trajectory(:,p);
                obj.dets2(p) = obj.arms(1)*obj.arms(2)*sin(point(2));
                [j,d] = jac(obj.arms(1),obj.arms(2),obj.arms(3),point(1),point(2),point(3));
                if p==1
                    obj.def(p) = 0;
                else
                    obj.def(p) = obj.dets2(p)-obj.dets2(p-1);
                end
                obj.kond2(p) = cond(j);
            end
        end
        function obj = toTSpace(obj)
            obj.transformed_pts_tsp = retransform(obj.arms,obj.trajectory); %from joint space to
        end
        function obj = computeError(obj)
            ori_pts = [];
            ori_pts(1,:) = sample_multiple(obj.original_pts(1,:), obj.sampleSize_2);
            ori_pts(2,:) = sample_multiple(obj.original_pts(2,:),obj.sampleSize_2);
            
            trs_pts = obj.transformed_pts_tsp;
            trs_pts(3,:) = [];
            abstand = sqrt(sum((trs_pts - ori_pts).^2));
            abstand_quadriert = sum((trs_pts-ori_pts).^2);
            obj.error = sum((trs_pts-ori_pts).^2);
            obj.meanError = sum(obj.error)/size(obj.error,2);
            obj.maxError = max(max(obj.error));
        end
        function plot(obj)
            
            ori_pts = [];
            ori_pts(1,:) = sample_multiple(obj.original_pts(1,:), obj.sampleSize_2);
            ori_pts(2,:) = sample_multiple(obj.original_pts(2,:),obj.sampleSize_2);
            
            trs_pts = obj.transformed_pts_tsp;
            trs_pts(3,:) = [];
            
            figure('visible','off');
            ax1 = subplot(2,1,1);
            l = length(obj.error);
            five_p = ceil(l*0.05);
            [sorted_values,sortIndex] = sort(obj.error(:),'descend');
            max_index = sortIndex(1:five_p);  %# Get a linear index into A of the 5 largest values
            
            h(1) = plot(obj.error);
            %line(xlim,[min_val,min_val])
            hold on;
            h(2) = plot(max_index, obj.error(max_index), 'd');
            str = strcat('Fehler (Mean =',num2str(double(obj.meanError)),')'); %double precision
            xlabel('Samplingpunkte')
            ylabel('Fehlergröße')
            legend(h, 'Fehler', '5% der Maximalwerte','Location','southeastoutside')
            title(ax1,str);
            
            ax2 = subplot(2,1,2);
            plot(obj.transformed_pts_tsp(1,:),obj.transformed_pts_tsp(2,:));
            hold on;
            plot(obj.original_pts(1,:),obj.original_pts(2,:));
            hold on;
            j(1) = plot(obj.start_pt(1),obj.start_pt(2), 'ob');
            hold on;
            j(2) = plot(obj.end_pt(1),obj.end_pt(2), 'og');
            xlabel('Weg in x-Richtung')
            ylabel('Weg in y-Richtung')
            legend(j, 'Startpunkt', 'Endpunkt','Location','southeastoutside')
            tit = sprintf('Originale und transformierte Punkte');
            if(obj.singularity == 1)
                tit = sprintf('Originale und transformierte Punkte, SINGULARITÄT');
            end
            title(ax2,tit);
            hold off;
            str = sprintf('test%d/iteration%d.jpg',obj.testnumber,obj.iteration);
            saveas(gcf,str)
            %berechne determinante an max_values
            figure('visible','off');
            plot(obj.dets2, obj.error)
            xlabel('Fehlergröße')
            ylabel('Determinante')
            title('Determinante und Fehler')
            str = sprintf('test%d/determinante_und_fehler%d.jpg',obj.testnumber,obj.iteration);
            saveas(gcf,str)
            
            figure('visible','off');
            plot(obj.dets2)
            xlabel('Samplingpunkte')
            ylabel('Determinante')
            title('Determinante')
            str = sprintf('test%d/determinante%d.jpg',obj.testnumber,obj.iteration);
            saveas(gcf,str)
            
            figure('visible','off');
            plot(obj.kond2,obj.error)
            xlabel('Fehlergröße')
            ylabel('Konditionszahl')
            title('Konditionszahl und Fehler')
            str = sprintf('test%d/konditionszahl_und_fehler%d.jpg',obj.testnumber,obj.iteration);
            saveas(gcf,str)
            
            figure('visible','off');
            plot(obj.kond2);
            xlabel('Samplingpunkte')
            ylabel('Konditionszahl')
            title('Konditionszahl')
            str = sprintf('test%d/konditionszahl%d.jpg',obj.testnumber,obj.iteration);
            saveas(gcf,str);
            
            %plot trajectory in joint space
            figure('visible','off');
            plot3(obj.trajectory(1,:),obj.trajectory(2,:),obj.trajectory(3,:))
            axis([-pi pi -pi pi -pi pi])
        end
        
    end
end

