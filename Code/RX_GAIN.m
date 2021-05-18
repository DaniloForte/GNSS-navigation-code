function [RX_GAIN] = RX_GAIN(off_bore_angle)

% function to compute the receiver gain from the off-bore direction assuming an
% antenna pointing toward the Earth center and with a gain pattern as in: 
% "Use of GNSS for lunar missions and plans for lunar in-orbit development"

% PROTOTYPE:
   % [RX_GAIN] = RX_GAIN(off_bore_angle);
    
% INPUT:
   % off_bore_angle: angle between the Earth center direction and the
   % transmitter satellites        [deg]
% OUTPUT:
    % RX_GAIN: antenna gain in dB  [dB]


% VERSIONS:
    % 1/4/2021: First version

    
% Data available up to 12.2 deg from boresight

x_RX = [-12.2,-11:1:11, 12.2]; % Off-bore angle values

y_RX = [ 3 3.2 5 7 8.5 10 11.1 12 12.8 13.3 13.7 13.9 14 ...
         13.9 13.7 13.3 12.8 12 11.1 10 8.5 7 5 3.2 3]; % Gain values
     
     
if off_bore_angle <= 12.2 % This is always the case for altitude above 20 Earth radius 
    
    RX_GAIN = interp1(x_RX,y_RX,off_bore_angle);
    
else % If the angle is greater than the defined interval return a very low gain (arbitrary)
    
    RX_GAIN = -6; 
    
end

end

