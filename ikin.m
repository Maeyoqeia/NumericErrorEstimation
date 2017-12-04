function [theta,theta2,singularity] = ikin(ox,oy,omega,length) %annahme - ox und oy sind an der base 0
%berechnet inverse kinematik nur für einen RRR-Roboter
%ox,oy, Koordinaten des Greifers
%omega-Rotation des Greifers, hier Winkel zur x-Achse
%length- armlängen
%output: beide Gelenkstellungen, die möglich sind (righty lefty solution)
% ox = sym(ox);
% oy = sym(oy);
% omega = sym(omega);
% length = sym(length);
%vpa increased variable precision (symbolic math toolbox) -default
%precision 32 digits
% ox = vpa(ox)
% oy = vpa(oy)
% omega = vpa(omega)
% length = vpa(length)
%an singularitäten steigt der fehler- wäre dann sinnvoll die rauszunehmen,
%oder? -> teste ob sin(theta2) = 0 ist, und wenn ja, nehme dann den anderen
%winkel
singularity = 0
o_wrist = [ox-length(3)*cos(omega);oy-length(3)*sin(omega);0];
l = sqrt(o_wrist(1)^2+o_wrist(2)^2);
stuff = (l^2-length(1)^2-length(2)^2)/(2*length(1)*length(2));
theta(2) = atan2(real(sqrt(1-stuff^2)),stuff);
theta2(2) = atan2(real(-sqrt(1-stuff^2)),stuff); 
IN = atan2(length(2)*sin(theta(2)),length(1)+length(2)*cos(theta(2))); 
OUT = atan2(o_wrist(2),o_wrist(1));
IN2 = atan2(length(2)*sin(theta2(2)),length(1)+length(2)*cos(theta2(2))); 
theta(1) = -IN+OUT ; %IN, IN2 negativ 
theta2(1) = -IN2+OUT;
if(abs(sin(theta(2))) < eps(0.75))
    singularity = 1
    %zweiten winkel nehmen?
end
theta(3) = omega-(theta(1)+theta(2));
theta2(3) = omega -(theta2(1)+theta2(2));
end