function [A_d] = FSP_Loss(d,f,c)

% Simple function to compute the Free Space Path(FSP) loss function of the path
% length, the signal frequency  and the speed of light in vacuum.

% PROTOTYPE:
    % [A_d] = FSP_Loss(d,f,c);
    
% INPUT:
   % d: transmitter and receiver antenna center distance      [km]
   % f: signal frequency                                      [Hz]
   % c: speed of light in vacuum                              [km/s]
% OUTPUT:
   % A_d: FSP loss in decibels                                [dB] 


% VERSIONS:
    % 1/4/2021: First version
    
    
% FSP losses in dB

A_d = 20*log10( (4 * pi)/c * d * f);


end

