function [D_shift] = doppler_shift(f_0,v_r,c)

% Function to compute the doppler shift due to the relative motion
% of two object.
% v_r should be positive if the two object are moving away from each other

% PROTOTYPE:
    % [D_shift] = doppler_shift(f_0,v_r,c);
    
% INPUT:
   % f_0: emitted signal frequency          [Hz]
   % v_r: transmitter and receiver antenna center redial relative velocity [km/s]
   % c: speed of light in vacuum            [km/s]
% OUTPUT:
   % D_shift: Doppler shift                 [Hz]


% VERSIONS:
   % 1/4/2021: First version

   
D_shift = - f_0 * (v_r/c);


end

