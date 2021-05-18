function [GPS,GALILEO,MULTI_GNSS] = ephemeridis_visibility(RX_Trajectory,GPS,GALILEO,MULTI_GNSS,aq_treshold,signalGPS,signalGALILEO)

% Function to perform ephemeridis visibility analysis for the GPS and
% Galielo constellations. The simulation data must have a 10s timestep.
% The demodulation thresholds and time are described in the thesis.

% PROTOTYPE:
   % [GPS,GALILEO,MULTI_GNSS] = ephemeridis_visibility('kernels\spk\Halo_LUMIO.bsp',GPS,GALILEO,MULTI_GNSS,15,'L1','E1');
    
% INPUT:
   % RX_Trajectory: path of the receiver trajectory SPK
   % GPS,GALILEO,MULTI_GNSS: structures containing the simulation
   % results from GNSS_visibility_function. NOTE(!!!) a 10s timestep is
   % needed for the simulation data.
   % aq_treshold: receiver acquistion/tracking threshold  [dB-Hz]
   % signalGPS,signalGALILEO: strings containing the selected GPS and Galileo signal band (should be
   % consistent with the simulation reults used)
% OUTPUT:
   % GPS,GALILEO,MULTI_GNSS: the imput structures at which the results of
   % the ephemerides visibility have been added


% VERSIONS:
   % 1/4/2021: First version


% Loading the required kernels.  We need the leapseconds kernel to convert input
% UTC time strings into ET. We also need SPK files for the GNSS constellations 
% and the the solar system planets and satelittes.

METAKR = 'kernels\GNSS_visibility.tm';
cspice_furnsh ( METAKR );

% We also need the spk of the receiver trajectory 
cspice_furnsh (RX_Trajectory); 

ID = cspice_spkobj(RX_Trajectory,1);
ID_RX = char(string(ID));


%% GPS 

% L1 signal

if isequal(signalGPS,'L1')
    
    GPS.ephemeridis_tracking = [];
    GPS.almanc_demodulated_flag = [];

    % For each SV

    for k = 1 : 30

        n = 0; % consecutive step above demodulation treshold counter
        m = 0; % ephemeridis avaiability counter 
        
        % For each timestep

        for r = 1 : numel(GPS.epoch)
            
            % if the signal in above the demodulation trehsold

            if   GPS.C_N0(r,k) > 26.5 

                n = n + 1;
                
            % if the signal in below the demodulation treshold

            elseif  GPS.C_N0(r,k) < 26.5  

                n = 0;

            end
            
            % if 5 consecutive steps (50s) are above the treshold the
            % ephemeridis are now avaible and are valid for 4h

            if n == 5

               m = 1;

            end
            
            % if 120 consecutive steps (1200s) are above the treshold the
            % almanac is now avaible and are valid for 12days            

            if n == 120

               GPS.almanc_demodulated_flag(r,k) = 1;

            end
            
            % the ephemeridis are valid for 4h (14400s) after they are
            % demodulated
            
            if m > 0 && m < 1440

                m = m + 1;

            elseif m == 1440

                m = 0;

            end
            
            % checking ephemeridis avaibility + tracking for each time step
            
            if m > 0 &&  GPS.C_N0(r,k) > aq_treshold

                GPS.ephemeridis_tracking(r,k) = 1;
            else
                
                GPS.ephemeridis_tracking(r,k) = 0;
                 
            end
        end
    end
    
% L5 signal

elseif isequal(signalGPS,'L5')
    
    GPS.ephemeridis_tracking = [];
    GPS.almanc_demodulated_flag = [];

    % For each SV

    for k = 1 : 30

        n = 0; % consecutive step above demodulation treshold counter
        m = 0; % ephemeridis avaiability counter 
        
        % For each timestep

        for r = 1 : numel(GPS.epoch)
            
            % if the signal in above the demodulation trehsold

            if   GPS.C_N0(r,k) > 26.1

                n = n + 1;
                
            % if the signal in below the demodulation treshold

            elseif  GPS.C_N0(r,k) < 26.1  

                n = 0;

            end
            
            % if 3 consecutive steps (24s) are above the treshold the
            % ephemeridis are now avaible and are valid for 4h

            if n == 3

               m = 1;

            end
            
            % if 60 consecutive steps (600s) are above the treshold the
            % almanac is now avaible and are valid for 12days            

            if n == 60

               GPS.almanc_demodulated_flag(r,k) = 1;

            end
            
            % the ephemeridis are valid for 4h (14400s) after they are
            % demodulated
            
            if m > 0 && m < 1440

                m = m + 1;

            elseif m == 1440

                m = 0;

            end
            
            % checking ephemeridis avaibility + tracking for each time step
            
            if m > 0 &&  GPS.C_N0(r,k) > aq_treshold

                GPS.ephemeridis_tracking(r,k) = 1;
            else
                
                GPS.ephemeridis_tracking(r,k) = 0;
                 
            end
        end
    end

end

%% GALILEO 

% E1 signal

if isequal(signalGALILEO,'E1')
    
    GALILEO.ephemeridis_tracking = [];
    GALILEO.almanc_demodulated_flag = [];

    % For each SV

    for k = 1 : 24  

        n = 0; % consecutive step above demodulation treshold counter
        m = 0; % ephemeridis avaiability counter 
        
        % For each timestep

        for r = 1 : numel(GALILEO.epoch)
            
            % if the signal in above the demodulation trehsold

            if   GALILEO.C_N0(r,k) > 27.7 

                n = n + 1;
                
            % if the signal in below the demodulation treshold

            elseif  GPS.C_N0(r,k) < 27.7  

                n = 0;

            end
            
            % if 3 consecutive steps (30s) are above the treshold the
            % ephemeridis are now avaible and are valid for 4h

            if n == 3

               m = 1;

            end
            
            % if 72 consecutive steps (720s) are above the treshold the
            % almanac is now avaible and are valid for 12days            

            if n == 72

               GALILEO.almanc_demodulated_flag(r,k) = 1;

            end
            
            % the ephemeridis are valid for 4h (14400s) after they are
            % demodulated
            
            if m > 0 && m < 1440

                m = m + 1;

            elseif m == 1440

                m = 0;

            end
            
            % checking ephemeridis avaibility + tracking for each time step
            
            if m > 0 &&  GALILEO.C_N0(r,k) > aq_treshold

                GALILEO.ephemeridis_tracking(r,k) = 1;
                
            else
                
                 GALILEO.ephemeridis_tracking(r,k) = 0;
                 
            end
        end
    end
    
% L5 signal

elseif isequal(signalGALILEO,'E5a')
    
    GALILEO.ephemeridis_tracking = [];
    GALILEO.almanc_demodulated_flag = [];

    % For each SV

    for k = 1 : 24

        n = 0; % consecutive step above demodulation treshold counter
        m = 0; % ephemeridis avaiability counter 
        
        % For each timestep

        for r = 1 : numel(GALILEO.epoch)
            
            % if the signal in above the demodulation trehsold

            if   GALILEO.C_N0(r,k) > 20.7

                n = n + 1;
                
            % if the signal in below the demodulation treshold

            elseif  GALILEO.C_N0(r,k) < 20.7  

                n = 0;

            end
            
            % if 5 consecutive steps (50s) are above the treshold the
            % ephemeridis are now avaible and are valid for 4h

            if n == 3

               m = 1;

            end
            
            % if 60 consecutive steps (600s) are above the treshold the
            % almanac is now avaible and are valid for 12days            

            if n == 60

               GALILEO.almanc_demodulated_flag(r,k) = 1;

            end
            
            % the ephemeridis are valid for 4h (14400s) after they are
            % demodulated
            
            if m > 0 && m < 1440

                m = m + 1;

            elseif m == 1440

                m = 0;

            end
            
            % checking ephemeridis avaibility + tracking for each time step
            
            if m > 0 &&  GALILEO.C_N0(r,k) > aq_treshold

                GALILEO.ephemeridis_tracking(r,k) = 1;
            else
                
                GALILEO.ephemeridis_tracking(r,k) = 0;
            end
        end
    end

end


%% GDOP

% The GDOP is now computed considering only the SV tracked for which the
% ephemeridis are avaible 


%% GPS
% At each time step 

for r = 1 : numel(GPS.epoch)

et = GPS.epoch(r);

% Number of acquired and tracked SV at the time considered

SV_tracked = sum(GPS.ephemeridis_tracking(r,:));

GPS.ephemeridis_n_tracked(r,1) = SV_tracked;


H = zeros(SV_tracked,4);

n = 1;

% For each SV of the constellation

    for k = 1 : 30

        % if the SV is being tracked

        if GPS.ephemeridis_tracking(r,k) > 0 && SV_tracked > 3

            GPSID = -(600 + k);
            GPSID = char(string(GPSID));  

            % Position vector of the tracked SV as seen from the TX

            a_n = cspice_spkpos ( GPSID, et,'J2000', 'XLT + S', ID_RX);

            % Normalizzation of the position vector

            a_n = a_n' / norm(a_n);

            H(n,:) = [a_n 1];

            n = n + 1;

        end

    end

        if SV_tracked > 3

            D = inv(H' * H);

            GDOP = sqrt(trace(D));

            GPS.ephemeridis_GDOP(r,1) = GDOP;

        else

            GPS.ephemeridis_GDOP(r,1) = nan;

        end
end

%% GALILEO

% At each time step 

for r = 1 : numel(GALILEO.epoch)

et = GALILEO.epoch(r);

% Number of acquired and tracked SV at the time considered

SV_tracked = sum(GALILEO.ephemeridis_tracking(r,:));

GALILEO.ephemeridis_n_tracked(r,1) = SV_tracked;


H = zeros(SV_tracked,4);

n = 1;

% For each SV of the constellation

    for k = 1 : 24

        % if the SV is being tracked

        if GALILEO.ephemeridis_tracking(r,k) > 0 && SV_tracked > 3

            GALILEOID = -(650 + k);
            GALILEOID = char(string(GALILEOID));  

            % Position vector of the tracked SV as seen from the TX

            a_n = cspice_spkpos ( GALILEOID, et,'J2000', 'XLT + S', ID_RX);

            % Normalizzation of the position vector

            a_n = a_n' / norm(a_n);

            H(n,:) = [a_n 1];

            n = n + 1;

        end

    end

        if SV_tracked > 3

            D = inv(H' * H);

            GDOP = sqrt(trace(D));

            GALILEO.ephemeridis_GDOP(r,1) = GDOP;

        else

            GALILEO.ephemeridis_GDOP(r,1) = nan;

        end
end


%% GPS + GALILEO

% At each time step 

for r = 1 : numel(GPS.epoch)
    
et = GPS.epoch(r);

% Number of acquired and tracked SV at the timeconsidered

SV_tracked_GPS = sum(GPS.ephemeridis_tracking(r,:));
SV_tracked_GALILEO = sum(GALILEO.ephemeridis_tracking(r,:));

SV_tracked = (SV_tracked_GPS + SV_tracked_GALILEO);

H = zeros(SV_tracked,4);

n = 1;

% For each SV of the GPS constellation

    for k = 1 : 30

        % if the SV is being tracked

        if GPS.ephemeridis_tracking(r,k) > 0 && SV_tracked > 3

            GPSID = -(600 + k);
            GPSID = char(string(GPSID));  

            % Position vector of the tracked SV as seen from the TX

            a_n = cspice_spkpos ( GPSID, et,'J2000', 'XLT + S', ID_RX);

            % Normalizzation of the position vector

            a_n = a_n' / norm(a_n);

            H(n,:) = [a_n 1];

            n = n + 1;

        end
        
    end
    
    % For each SV of the GALILEO constellation

    for k = 1 : 24

        % if the SV is being tracked

        if GALILEO.ephemeridis_tracking(r,k) > 0 && SV_tracked > 3

            GALILEOID = -(650 + k);
            GALILEOID = char(string(GALILEOID));  

            % Position vector of the tracked SV as seen from the TX

            a_n = cspice_spkpos ( GALILEOID, et,'J2000', 'XLT + S', ID_RX);

            % Normalizzation of the position vector

            a_n = a_n' / norm(a_n);

            H(n,:) = [a_n 1];

            n = n + 1;

        end
        
    end
    
    if SV_tracked > 3
        
        D = inv(H' * H);

        GDOP = sqrt(trace(D));

        MULTI_GNSS.ephemeridis_GPSandGAL_GDOP(r,1) = GDOP;
        
     else
        
        MULTI_GNSS.ephemeridis_GPSandGAL_GDOP(r,1) =  nan;
    
    end
end


%% clearing the kernel pool

cspice_kclear

end

