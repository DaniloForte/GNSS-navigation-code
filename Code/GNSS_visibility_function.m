function [GPS,GALILEO,GLONASS,BEIDOU,LNNSS_6SV,LNNSS_8SV,ID_RX,MULTI_GNSS] = GNSS_visibility_function(RX_Trajectory,tstep,aq_treshold,signalGPS,signalGALILEO,GLO_flag,BEI_flag,LNNSS_flag,LNNSS,Simulation_time)

% Main function to compute the time window in which the GNSS satellites
% are in view (no occultation from Earth and Moon), the distance, relative
% velocity and angular separation in order to compute the C/N0 of the link
% budget. Also the GDOP is computed for each constellation and some
% combination. Also the LNNSS is included

% PROTOTYPE:
   % [GPS,GALILEO,GLONASS,BEIDOU,LNNSS_6SV,LNNSS_8SV,ID_RX,MULTI_GNSS] = GNSS_visibility_function('kernels\spk\MTO.bsp',600,15,'L1','E1',1,1,1,'H');
    
% INPUT:
   % RX_Trajectory: path of the receiver trajectory SPK
   % tstep: simulation timestep                           [s]
   % aq_treshold: receiver acquistion/tracking threshold  [dB-Hz]
   % signalGPS,signalGALILEO: strings containing the selected GPS and Galileo signal band
   % GLO_flag,BEI_flag: flag to simulate also GLONASS and BEIDOU Constellations (1:YES, 0:NO)
   % LNNSS_flag:flag to simulate also the LNNSS Constellations (1:YES, 0:NO)
   % LNNSS: string to specify which LNNSS architecture to simulate (L = Light, H = heavy)
   % Simulation_time: if specified select a simulation timespan   [s]
   
% OUTPUT:
   % GPS,GALILEO,GLONASS,BEIDOU,LNNSS_6SV,LNNSS_8SV,MULTI_GNSS: structures containting all the simuations results
   % ID_RX: receiver SPK ID


% VERSIONS:
   % 1/4/2021: First version



%% Local parameters
c = cspice_clight; % Speed of light [km/s]
dpr = cspice_dpr; % degree per radiant

METAKR = 'kernels\GNSS_visibility.tm';

% Load the kernels that this program requires.  We need a leapseconds kernel to convert input
% UTC time strings into ET.  We also need SPK files for the GNSS constellations and the the solar system planets and satellites.

cspice_furnsh ( METAKR );

% We also need the spk of the receiver trajectory 

cspice_furnsh (RX_Trajectory); 

% Recovering the id and the time coverage of the receiver trajectory

ID = cspice_spkobj(RX_Trajectory,1);


cover = cspice_spkcov(RX_Trajectory,ID,1);

timstr = cspice_timout( cover(:,1)','YYYY MON DD HR:MN:SC.### (TDB) ::TDB' );

% Returning the coverage in the command window

 fprintf( '========================================\n')
 fprintf( 'Coverage for object %d\n', ID )

fprintf('   Start: %s\n'  ,timstr(1,:) )
fprintf('   Stop:  %s\n\n' , timstr(2,:) )


if ID == -515 % For MMS ( half an orbit)

    et_span = [cspice_str2et('2020 DEC 10 12:00:01.000  TDB'),cspice_str2et('2020 DEC 12 06:01:31.339  TDB')];
    LNNSS = 0;
 
elseif ID == -516 % For Moon transfer orbit ( One orbit (about 10 days))
    
     et_span = [cover(1) + 1 , cover(2)]; 
     
     if nargin == 10 % If the Integration time has been specified in the imput 
        et_span = [et_span(1), et_span(1) + Simulation_time]; 
     end
    

elseif ID == -517 % For Lunar pathfinder (one E-M synodic period 29.48 days)
    
    et_span = [cspice_str2et('2022-DEC-1, 12:00:00  TDB'),cspice_str2et('2022-DEC-1, 12:00:00  TDB') + 29.48*86400];
    
     if nargin == 10 % If the Integration time has been specified in the imput 
        et_span = [et_span(1), et_span(1) + Simulation_time]; 
     end
    
elseif ID == -100099 % For LUMIO (one E-M synodic period 29.48 days)
    
    et_span = [cspice_str2et('2024-MAR-21, 12:00:00  TDB'),cspice_str2et('2024-MAR-21, 12:00:00  TDB') + 29.48*86400 ];
    
     if nargin == 10 % If the Integration time has been specified in the imput 
        et_span = [et_span(1), et_span(1) + Simulation_time]; 
     end   
     
else  % Other RX trajectory spk
    
         et_span = [cover(1) + 1 ,cover(2)]; 
         
      if nargin == 10 % If the Integration time has been specified in the imput 
        et_span = [et_span(1), et_span(1) + Simulation_time]; 
     end
end

  
% To reduce the time space to explore for each GNSS SV, we first compute
% the time window where the SV as seen from the RX is not occulted by the
% Earth

% Computation of occultation of the GNSS satellites by the EARTH as seen
% from the TX

% Any occultations lasting less than 10 minutes or less than a timestep will be ignored .

step_OCC = max(600,tstep);
MAXWIN  =  10000;

ID_RX = char(string(ID));

cnfine = cspice_wninsd( et_span(1), et_span(2));

% Time span discretization 

tspan = et_span(1) : tstep : et_span(2);


%% GPS

GPS = struct('epoch',[],'RX_TXdistance',[],'RX_TXv_radial',[],'RX_offbore',[],'TX_offbore',[],'TX_azimuth',[],'visibility_flag',[]); % Struct containing the results

for k = 1:30
    
    GPSID = -(600 + k);
    GPSID = char(string(GPSID));  
    
    
    % Time windows of Earth occultation
    
    occtyp  = 'any';
    front   = 'earth';
    fshape  = 'ellipsoid';
    fframe  = 'ITRF93';
    back    = GPSID;
    bshape  = 'point';
    bframe  = '-';
    obsrvr  = ID_RX;
    abcorr  = 'LT+S';

   
    EARTH_occultation = cspice_gfoclt( occtyp, front, fshape, fframe, ...
                          back, bshape,bframe,          ...
                          abcorr, obsrvr, step_OCC, cnfine,  ...
                          MAXWIN);
                     
   
                      
    % List the beginning and ending times in each interval
    % if result contains data.

   
    left_E = zeros(1,numel( EARTH_occultation)/2 + 1);
    right_E = zeros(1,numel( EARTH_occultation)/2 + 1);
    
    for i = 1:numel(EARTH_occultation)/2

     [left_E(i), right_E(i)] = cspice_wnfetd( EARTH_occultation, i );

    end
    
    
    % Adding to the windows the endtime of the simulation, needed to
    % compute after the last occultation window or if the window is empty
    
    left_E(end) = tspan(end);
    right_E(end) = tspan(end);
    
    
    % Time windows of Moon occultation
    
    occtyp  = 'any';
    front   = 'moon';
    fshape  = 'ellipsoid';
    fframe  = 'moon_PA';
    back    = GPSID;
    bshape  = 'point';
    bframe  = '-';
    obsrvr  = ID_RX;
    abcorr  = 'LT+S';

   
    MOON_occultation = cspice_gfoclt( occtyp, front, fshape, fframe, ...
                          back, bshape,bframe,          ...
                          abcorr, obsrvr, step_OCC, cnfine,  ...
                          MAXWIN);
                     
   
                      
    % List the beginning and ending times in each interval
    % if result contains data.

   
    
    left_M = zeros(1,numel( MOON_occultation)/2 + 1);
    right_M = zeros(1,numel( MOON_occultation)/2 + 1);
    
    for i = 1:numel(MOON_occultation)/2

     [left_M(i), right_M(i)] = cspice_wnfetd( MOON_occultation, i );

    end
    
    % Adding to the windows the endtime of the simulation, needed to
    % compute after the last occultation window or if the window is empty
    
    left_M(end) = tspan(end);
    right_M(end) = tspan(end);
     
    
    
    % For each time step where the TX is in view of the SV, the ideal
    % range, relative velocity and angular separation with respect to the
    % SV nadir in evaluated
    
    % initializzation of the Earth and Moon occultation windows counter
    
     n = 1; % Earth occultation windows counter
     m = 1; % Moon occultation windows counter
    
     
    for r = 1 : numel(tspan) % for each time step
        
         et = tspan(r);
         
         % saving the epoch in the output structure
         
         GPS.epoch(r,1) = et ;
         
         
         % If the timestep is not inside an occultation window the
         % computation are performed
         
         if et <= left_E(n) && et <= left_M(m) 
             
             % Computing all the relevant position and state vector at the
             % current timestep 
           
             % State of the RX as seen from the TX SV
             
             [state,ltime] = cspice_spkezr (ID_RX, et,'J2000', 'XLT', GPSID);

             TX2RX_pos_J2000 = state(1:3,1); % Position vector
             TX2RX_vel_J2000 = state(4:6,1); % Velocity vector
             
             % Position of Earth as seen from the TX SV ( an its nadir
             % pointing)
             
             TX2EA_pos_J2000 = cspice_spkpos ( 'EARTH', et  ,'J2000', 'LT', GPSID);
             
             % Position of the SUN as seen from the TX SV 
             
             TX2SUN_pos_J2000 = cspice_spkpos ( 'SUN', et  ,'J2000', 'LT', GPSID); 
             
             % Only the Sun direction is needed

             TX2SUN_pos_J2000 = cspice_vhat(TX2SUN_pos_J2000);
             
             % Position of Earth as seen from the RX ( an its high gain
             % antenna pointing)
             
             RX2EA_pos_J2000 = cspice_spkpos ( 'EARTH', et ,'J2000', 'LT', ID_RX);

             % We can compute the distance between the TX and RX as the norm of the position vector    
             
             d = norm(TX2RX_pos_J2000);
             
             range = ltime * c; % equivalently is the light time multiplied by the speed of light

             % Radial velocity from the SV to the RX
             
             v_radial = dot(TX2RX_vel_J2000,(TX2RX_pos_J2000/d));
             
             % To compute the off-bore angle and the azimuth of RX wrt 3D
             % antenna pattern of each SV a simple yaw-steering (YS) attitude model is
             % implemented 
             
             % For the GPS constellation two different body reference frame
             % are in use depending on the SV block.
            
             % For IIR and IIR-M block the x_YS axis is pointed opposite to
             % the Sun direction.
             % For IIF and III block the x_YS axis is pointed toward the Sun direction.
             
             % From GPSID -601 to - 608 --> BIIR   SV 
             % From GPSID -609 to - 615 --> BIIR-M SV
             % From GPSID -616 to - 627 --> BII-F  SV
             % From GPSID -628 to - 630 --> BIII   SV
             
             if k < 16 % IIR/IIR-M block
                 
                 % In nominal operation the z-axis is pointing the nadir
                 % direcion toward Earth, the x-axis is opposite to the Sun
                 % direction and the y-axis complete the othonormal basis

                 k_YS = cspice_vhat(TX2EA_pos_J2000);

                 j_YS = cross(-TX2SUN_pos_J2000,TX2EA_pos_J2000)/norm(cross(-TX2SUN_pos_J2000,TX2EA_pos_J2000));

                 i_YS = cross(j_YS,k_YS);
                 
                 % A rotation matrix can be compute to rotate a vector in
                 % the SV YS fram ( ideally equivalent to the body frame)

                 R_J20002YS = [i_YS'; j_YS'; k_YS'];
                 
                 % The direction vector of the RX in the TX body frame can
                 % be computed rotating the TX2RX position vector from the
                 % J2000 frame
                 
                 TX2RX_pos_YS = R_J20002YS * cspice_vhat(TX2RX_pos_J2000); 
                 
                 % The angle can be computed as a transformation into
                 % spherical coordinates where:
                 
                 % Theta is the angle between the z-axis ( nadir pointing)
                 % and the RX direction in the YS_frame. This is the
                 % off-bore angle. It is bounded in [0, pi];
                 
                 % Phi is the angle between the +x_YS axis and the
                 % projection of the RX direction into the xy plane. It is
                 % consistend with the convention for the azimuth of the
                 % directivity pattern and is the azimuth. It is bounded in [-pi, pi];
           
                 
                 [~,theta,phi] = cspice_recsph(TX2RX_pos_YS);
                 
                 TX2RX_offbore = cspice_convrt ( theta,'RADIANS', 'DEGREES');
                 TX2RX_azimuth = cspice_convrt ( phi,'RADIANS', 'DEGREES');
                 
             else % IIF/III block 
                 
                 % In nominal operation the z-axis is pointing the nadir
                 % direcion toward Earth, the x-axis is toward the Sun
                 % direction and the y-axis complete the othonormal basis
                 % (then the computation are the same)
       

                 k_YS = cspice_vhat(TX2EA_pos_J2000);

                 j_YS = cross(TX2SUN_pos_J2000,TX2EA_pos_J2000)/norm(cross(TX2SUN_pos_J2000,TX2EA_pos_J2000));

                 i_YS = cross(j_YS,k_YS);

                 R_J20002YS = [i_YS'; j_YS'; k_YS'];
                                 
                 TX2RX_pos_YS = R_J20002YS * cspice_vhat(TX2RX_pos_J2000); 
                               
                 [~,theta,phi] = cspice_recsph(TX2RX_pos_YS);
                 
                 TX2RX_offbore = cspice_convrt ( theta,'RADIANS', 'DEGREES');
                 TX2RX_azimuth = cspice_convrt ( phi,'RADIANS', 'DEGREES');

             
             end
             
             % The Gain pattern data has the azimuth defined in the range
             % [0 2pi]
             
             if  TX2RX_azimuth < 0
                 
                  TX2RX_azimuth = 2*pi*dpr + TX2RX_azimuth;
                  
             end
             
             % Angular separation between the Earth pointing antenna of the
             % receiver SPC and the position of the SV ( off-bore angle)


             RX2TX_offbore = cspice_convrt ( cspice_vsep(RX2EA_pos_J2000, -TX2RX_pos_J2000),'RADIANS', 'DEGREES'); 


             GPS.RX_TXdistance(r,k) = d ;
             GPS.RX_TXv_radial(r,k) = v_radial ;
             GPS.TX_offbore(r,k) = TX2RX_offbore ;
             GPS.TX_azimuth(r,k) = TX2RX_azimuth ;
             GPS.RX_offbore(r,k) = RX2TX_offbore ;
             GPS.visibility_flag(r,k) = 1 ;
                 
         else
             
             % else if the RX is in occultation or by Earth or Moon the SV
             % is not visible
             
             GPS.visibility_flag(r,k) = 0 ;

         end
         
         % Once the time step surpass a time window the counter is increased to analyze the next one
         
         if et >= right_E(n) && n < length(right_E)
             n = n + 1;
         end
         
         if et >= right_M(m) && m < length(right_M)
             m = m + 1;
         end
         

    end
     
    
    % In the first step the occultation flag is always 1 even if there is
    % no visibility due to the way the logic check for the spice
    % windows is perfomed, so an extra check for the first row is needed for
    % each SV
    
    if (left_E(1) - GPS.epoch(1,1)) == 0
        
        GPS.visibility_flag(1,k) = 0 ;
        
    end
        
        
         
                      
end

%% GALILEO

GALILEO = struct('epoch',[],'RX_TXdistance',[],'RX_TXv_radial',[],'RX_offbore',[],'TX_offbore',[],'TX_azimuth',[],'visibility_flag',[]); % Struct containing the results

for k = 1:24
    
    GALILEOID = -(650 + k);
    GALILEOID = char(string(GALILEOID));  
    
    
    % Time windows of Earth occultation
    
    occtyp  = 'any';
    front   = 'earth';
    fshape  = 'ellipsoid';
    fframe  = 'ITRF93';
    back    = GALILEOID;
    bshape  = 'point';
    bframe  = '-';
    obsrvr  = ID_RX;
    abcorr  = 'LT+S';

   
    EARTH_occultation = cspice_gfoclt( occtyp, front, fshape, fframe, ...
                          back, bshape,bframe,          ...
                          abcorr, obsrvr, step_OCC, cnfine,  ...
                          MAXWIN);
                     
   
                      
    % List the beginning and ending times in each interval
    % if result contains data.

   
    left_E = zeros(1,numel( EARTH_occultation)/2 + 1);
    right_E = zeros(1,numel( EARTH_occultation)/2 + 1);
    
    for i = 1:numel(EARTH_occultation)/2

     [left_E(i), right_E(i)] = cspice_wnfetd( EARTH_occultation, i );

    end
    
    
    % Adding to the windows the endtime of the simulation, needed to
    % compute after the last occultation window or if the window is empty
    
    left_E(end) = tspan(end);
    right_E(end) = tspan(end);
    
    
    % Time windows of Moon occultation
    
    occtyp  = 'any';
    front   = 'moon';
    fshape  = 'ellipsoid';
    fframe  = 'moon_PA';
    back    = GALILEOID;
    bshape  = 'point';
    bframe  = '-';
    obsrvr  = ID_RX;
    abcorr  = 'LT+S';

   
    MOON_occultation = cspice_gfoclt( occtyp, front, fshape, fframe, ...
                          back, bshape,bframe,          ...
                          abcorr, obsrvr, step_OCC, cnfine,  ...
                          MAXWIN);
                     
   
                      
    % List the beginning and ending times in each interval
    % if result contains data.

   
    
    left_M = zeros(1,numel( MOON_occultation)/2 + 1);
    right_M = zeros(1,numel( MOON_occultation)/2 + 1);
    
    for i = 1:numel(MOON_occultation)/2

     [left_M(i), right_M(i)] = cspice_wnfetd( MOON_occultation, i );

    end
    
    % Adding to the windows the endtime of the simulation, needed to
    % compute after the last occultation window or if the window is empty
    
    left_M(end) = tspan(end);
    right_M(end) = tspan(end);
     
    
    
    % For each time step where the TX is in view of the SV, the ideal
    % range, relative velocity and angular separation with respect to the
    % SV nadir in evaluated
    
    % initializzation of the Earth and Moon occultation windows counter
    
     n = 1; % Earth occultation windows counter
     m = 1; % Moon occultation windows counter
    
     
    for r = 1 : numel(tspan) % for each time step
        
         et = tspan(r);
         
         % saving the epoch in the output structure
         
         GALILEO.epoch(r,1) = et ;
         
         
         % If the timestep is not inside an occultation window the
         % computation are performed
         
         if et <= left_E(n) && et <= left_M(m)
             
             % Computing all the relevant position and state vector at the
             % current timestep 
           
             % State of the RX as seen from the TX SV
             
             [state,ltime] = cspice_spkezr (ID_RX, et,'J2000', 'XLT', GALILEOID);

             TX2RX_pos_J2000 = state(1:3,1); % Position vector
             TX2RX_vel_J2000 = state(4:6,1); % Velocity vector
             
             % Position of Earth as seen from the TX SV ( an its nadir
             % pointing)
             
             TX2EA_pos_J2000 = cspice_spkpos ( 'EARTH', et  ,'J2000', 'LT', GALILEOID);
             
             % Position of the SUN as seen from the TX SV 
             
             TX2SUN_pos_J2000 = cspice_spkpos ( 'SUN', et  ,'J2000', 'LT', GALILEOID); 
             
             % Only the Sun direction is needed 

             TX2SUN_pos_J2000 = cspice_vhat(TX2SUN_pos_J2000); % vector normalizzation 
             
             % Position of Earth as seen from the RX ( an its high gain
             % antenna pointing)
             
             RX2EA_pos_J2000 = cspice_spkpos ( 'EARTH', et ,'J2000', 'LT', ID_RX);

             % We can compute the distance between the TX and RX as the norm of the position vector    
             
             d = norm(TX2RX_pos_J2000);
             
             range = ltime * c; % equivalently is the light time multiplied by the speed of light

             % Radial velocity from the SV to the RX
             
             v_radial = dot(TX2RX_vel_J2000,(TX2RX_pos_J2000/d));
             
             % To compute the off-bore angle and the azimuth of RX wrt 3D
             % antenna pattern of each SV a simple yaw-steering (YS) attitude model is
             % implemented 
             
             % For the GALILEO constellation as no precise information were found, each Sv
             % is assumed as a BIIF/BIII' GPS SV .
            
             % For IIF and III block the x_YS axis is pointed toward the Sun direction.
                              
             % In nominal operation the z-axis is pointing the nadir
             % direcion toward Earth, the x-axis is toward the Sun
             % direction and the y-axis complete the othonormal basis
            

             k_YS = cspice_vhat(TX2EA_pos_J2000);

             j_YS = cross(TX2SUN_pos_J2000,TX2EA_pos_J2000)/norm(cross(TX2SUN_pos_J2000,TX2EA_pos_J2000));

             i_YS = cross(j_YS,k_YS);
             
             % A rotation matrix can be compute to rotate a vector in
             % the SV YS fram ( ideally equivalent to the body frame)

             R_J20002YS = [i_YS'; j_YS'; k_YS'];
             
             % The direction vector of the RX in the TX body frame can
             % be computed rotating the RX2RX position vector from the
             % J2000 frame 

             TX2RX_pos_YS = R_J20002YS * cspice_vhat(TX2RX_pos_J2000); 
             
             % The angle can be computed as a transformation into
             % spherical coordinates where:

             % Theta is the angle between the z-axis ( nadir pointing)
             % and the RX direction in the YS_frame. This is the
             % off-bore angle. It is bounded in [0, pi];

             % Phi is the angle between the +x_YS axis and the
             % projection of the RX direction into the xy plane. It is
             % consistend with the convention for the azimuth of the
             % directivity pattern and is the azimuth. It is bounded in [-pi, pi];

             [~,theta,phi] = cspice_recsph(TX2RX_pos_YS);

             TX2RX_offbore = cspice_convrt ( theta,'RADIANS', 'DEGREES');
             TX2RX_azimuth = cspice_convrt ( phi,'RADIANS', 'DEGREES');

             
             
             % The Gain pattern data has the azimuth defined in the range
             % [0 2pi]
             
             if  TX2RX_azimuth < 0
                 
                  TX2RX_azimuth = 2*pi*dpr + TX2RX_azimuth;
                  
             end
             
             % Angular separation between the Earth pointing antenna of the
             % receiver SPC and the position of the SV ( off-bore angle)

             RX2TX_offbore = cspice_convrt ( cspice_vsep(RX2EA_pos_J2000, -TX2RX_pos_J2000),'RADIANS', 'DEGREES'); 


             GALILEO.RX_TXdistance(r,k) = d ;
             GALILEO.RX_TXv_radial(r,k) = v_radial ;
             GALILEO.TX_offbore(r,k) = TX2RX_offbore ;
             GALILEO.TX_azimuth(r,k) = TX2RX_azimuth ;
             GALILEO.RX_offbore(r,k) = RX2TX_offbore ;
             GALILEO.visibility_flag(r,k) = 1 ;
  
             
         else
             
             % else if the RX is in occultation or by Earth or Moon the SV
             % is not visible
             
             GALILEO.visibility_flag(r,k) = 0 ;

         end
         
         % Once the time step surpass a time window the counter is increased to analyze the next one
         
         if et >= right_E(n) && n < length(right_E)
             n = n + 1;
         end
         
         if et >= right_M(m) && m < length(right_M)
             m = m + 1;
         end
         

    end
     
    % In the first step the occultation flag is always 1 even if there is
    % no visibility due to the way the logic check for the spice
    % windows is perfomed, so an extra check for the first row is needed for
    % each SV
    
    if (left_E(1) - GALILEO.epoch(1,1)) == 0
        
        GALILEO.visibility_flag(1,k) = 0 ;
        
    end

       
                      
end


%% GLONASS



GLONASS = struct('epoch',[],'RX_TXdistance',[],'RX_TXv_radial',[],'RX_offbore',[],'TX_offbore',[],'TX_azimuth',[],'visibility_flag',[]); % Struct containing the results

if GLO_flag == 1
    
for k = 1:24
    
    GLONASSID = -(700 + k);
    GLONASSID = char(string(GLONASSID));  
    
    
    % Time windows of Earth occultation
    
    occtyp  = 'any';
    front   = 'earth';
    fshape  = 'ellipsoid';
    fframe  = 'ITRF93';
    back    = GLONASSID;
    bshape  = 'point';
    bframe  = '-';
    obsrvr  = ID_RX;
    abcorr  = 'LT+S';

   
    EARTH_occultation = cspice_gfoclt( occtyp, front, fshape, fframe, ...
                          back, bshape,bframe,          ...
                          abcorr, obsrvr, step_OCC, cnfine,  ...
                          MAXWIN);
                     
   
                      
    % List the beginning and ending times in each interval
    % if result contains data.

   
    left_E = zeros(1,numel( EARTH_occultation)/2 + 1);
    right_E = zeros(1,numel( EARTH_occultation)/2 + 1);
    
    for i = 1:numel(EARTH_occultation)/2

     [left_E(i), right_E(i)] = cspice_wnfetd( EARTH_occultation, i );

    end
    
    
    % Adding to the windows the endtime of the simulation, needed to
    % compute after the last occultation window or if the window is empty
    
    left_E(end) = tspan(end);
    right_E(end) = tspan(end);
    
    
    % Time windows of Moon occultation
    
    occtyp  = 'any';
    front   = 'moon';
    fshape  = 'ellipsoid';
    fframe  = 'moon_PA';
    back    = GLONASSID;
    bshape  = 'point';
    bframe  = '-';
    obsrvr  = ID_RX;
    abcorr  = 'LT+S';

   
    MOON_occultation = cspice_gfoclt( occtyp, front, fshape, fframe, ...
                          back, bshape,bframe,          ...
                          abcorr, obsrvr, step_OCC, cnfine,  ...
                          MAXWIN);
                     
   
                      
    % List the beginning and ending times in each interval
    % if result contains data.

   
    
    left_M = zeros(1,numel( MOON_occultation)/2 + 1);
    right_M = zeros(1,numel( MOON_occultation)/2 + 1);
    
    for i = 1:numel(MOON_occultation)/2

     [left_M(i), right_M(i)] = cspice_wnfetd( MOON_occultation, i );

    end
    
    % Adding to the windows the endtime of the simulation, needed to
    % compute after the last occultation window or if the window is empty
    
    left_M(end) = tspan(end);
    right_M(end) = tspan(end);
     
    
    
    % For each time step where the TX is in view of the SV, the ideal
    % range, relative velocity and angular separation with respect to the
    % SV nadir in evaluated
    
    % initializzation of the Earth and Moon occultation windows counter
    
     n = 1; % Earth occultation windows counter
     m = 1; % Moon occultation windows counter
    
     
    for r = 1 : numel(tspan) % for each time step
        
         et = tspan(r);
         
         % saving the epoch in the output structure
         
         GLONASS.epoch(r,1) = et ;
         
         
         % If the timestep is not inside an occultation window the
         % computation are performed
         
         if et <= left_E(n) && et <= left_M(m)

             % Computing all the relevant position and state vector at the
             % current timestep 
           
             % State of the RX as seen from the TX SV
             
             [state,ltime] = cspice_spkezr (ID_RX, et,'J2000', 'XLT', GLONASSID);

             TX2RX_pos_J2000 = state(1:3,1); % Position vector
             TX2RX_vel_J2000 = state(4:6,1); % Velocity vector
             
             % Position of Earth as seen from the TX SV ( an its nadir
             % pointing)
             
             TX2EA_pos_J2000 = cspice_spkpos ( 'EARTH', et  ,'J2000', 'LT', GLONASSID);
             
             % Position of the SUN as seen from the TX SV 
             
             TX2SUN_pos_J2000 = cspice_spkpos ( 'SUN', et  ,'J2000', 'LT', GLONASSID); 
             
             % Only the Sun direction is needed 

             TX2SUN_pos_J2000 = cspice_vhat(TX2SUN_pos_J2000); % vector normalizzation 
             
             % Position of Earth as seen from the RX ( an its high gain
             % antenna pointing)
             
             RX2EA_pos_J2000 = cspice_spkpos ( 'EARTH', et ,'J2000', 'LT', ID_RX);

             % We can compute the distance between the TX and RX as the norm of the position vector    
             
             d = norm(TX2RX_pos_J2000);
             
             range = ltime * c; % equivalently is the light time multiplied by the speed of light

             % Radial velocity from the SV to the RX
             
             v_radial = dot(TX2RX_vel_J2000,(TX2RX_pos_J2000/d));
             
             % To compute the off-bore angle and the azimuth of RX wrt 3D
             % antenna pattern of each SV a simple yaw-steering (YS) attitude model is
             % implemented 
             
             % For the GLONASS constellation as no precise information were found, each Sv
             % is assumed as a BIIF/BIII' GPS SV .
            
             % For IIF and III block the x_YS axis is pointed toward the Sun direction.
                              
             % In nominal operation the z-axis is pointing the nadir
             % direcion toward Earth, the x-axis is toward the Sun
             % direction and the y-axis complete the othonormal basis
            

             k_YS = cspice_vhat(TX2EA_pos_J2000);

             j_YS = cross(TX2SUN_pos_J2000,TX2EA_pos_J2000)/norm(cross(TX2SUN_pos_J2000,TX2EA_pos_J2000));

             i_YS = cross(j_YS,k_YS);
             
             % A rotation matrix can be compute to rotate a vector in
             % the SV YS fram ( ideally equivalent to the body frame)

             R_J20002YS = [i_YS'; j_YS'; k_YS'];
             
             % The direction vector of the RX in the TX body frame can
             % be computed rotating the RX2RX position vector from the
             % J2000 frame 

             TX2RX_pos_YS = R_J20002YS * cspice_vhat(TX2RX_pos_J2000); 
             
             % The angle can be computed as a transformation into
             % spherical coordinates where:

             % Theta is the angle between the z-axis ( nadir pointing)
             % and the RX direction in the YS_frame. This is the
             % off-bore angle. It is bounded in [0, pi];

             % Phi is the angle between the +x_YS axis and the
             % projection of the RX direction into the xy plane. It is
             % consistend with the convention for the azimuth of the
             % directivity pattern and is the azimuth. It is bounded in [-pi, pi];

             [~,theta,phi] = cspice_recsph(TX2RX_pos_YS);

             TX2RX_offbore = cspice_convrt ( theta,'RADIANS', 'DEGREES');
             TX2RX_azimuth = cspice_convrt ( phi,'RADIANS', 'DEGREES');

             
             
             % The Gain pattern data has the azimuth defined in the range
             % [0 2pi]
             
             if  TX2RX_azimuth < 0
                 
                  TX2RX_azimuth = 2*pi*dpr + TX2RX_azimuth;
                  
             end
             
             % Angular separation between the Earth pointing antenna of the
             % receiver SPC and the position of the SV ( off-bore angle)

             RX2TX_offbore = cspice_convrt ( cspice_vsep(RX2EA_pos_J2000, -TX2RX_pos_J2000),'RADIANS', 'DEGREES'); 


             GLONASS.RX_TXdistance(r,k) = d ;
             GLONASS.RX_TXv_radial(r,k) = v_radial ;
             GLONASS.TX_offbore(r,k) = TX2RX_offbore ;
             GLONASS.TX_azimuth(r,k) = TX2RX_azimuth ;
             GLONASS.RX_offbore(r,k) = RX2TX_offbore ;
             GLONASS.visibility_flag(r,k) = 1 ;
                 
         else
             
             % else if the RX is in occultation or by Earth or Moon the SV
             % is not visible
             
             GLONASS.visibility_flag(r,k) = 0 ;

         end
         
         % Once the time step surpass a time window the counter is increased to analyze the next one
         
         if et >= right_E(n) && n < length(right_E)
             n = n + 1;
         end
         
         if et >= right_M(m) && m < length(right_M)
             m = m + 1;
         end
         

    end
     
    % In the first step the occultation flag is always 1 even if there is
    % no visibility due to the way spice the logic check for the spice
    % windows is perfomed, so an extra check for the first row need for
    % each SV
    
    if (left_E(1) - GLONASS.epoch(1,1)) == 0
        
        GLONASS.visibility_flag(1,k) = 0 ;
        
    end

       
                      
end
end

%% BEIDOU-3

BEIDOU = struct('epoch',[],'RX_TXdistance',[],'RX_TXv_radial',[],'RX_offbore',[],'TX_offbore',[],'TX_azimuth',[],'visibility_flag',[]); % Struct containing the results

if BEI_flag == 1
    
for k = 1:30
    
    BEIDOUID = -(750 + k);
    BEIDOUID = char(string(BEIDOUID));  
    
    
    % Time windows of Earth occultation
    
    occtyp  = 'any';
    front   = 'earth';
    fshape  = 'ellipsoid';
    fframe  = 'ITRF93';
    back    = BEIDOUID;
    bshape  = 'point';
    bframe  = '-';
    obsrvr  = ID_RX;
    abcorr  = 'LT+S';

   
    EARTH_occultation = cspice_gfoclt( occtyp, front, fshape, fframe, ...
                          back, bshape,bframe,          ...
                          abcorr, obsrvr, step_OCC, cnfine,  ...
                          MAXWIN);
                     
   
                      
    % List the beginning and ending times in each interval
    % if result contains data.

   
    left_E = zeros(1,numel( EARTH_occultation)/2 + 1);
    right_E = zeros(1,numel( EARTH_occultation)/2 + 1);
    
    for i = 1:numel(EARTH_occultation)/2

     [left_E(i), right_E(i)] = cspice_wnfetd( EARTH_occultation, i );

    end
    
    
    % Adding to the windows the endtime of the simulation, needed to
    % compute after the last occultation window or if the window is empty
    
    left_E(end) = tspan(end);
    right_E(end) = tspan(end);
    
    
    % Time windows of Moon occultation
    
    occtyp  = 'any';
    front   = 'moon';
    fshape  = 'ellipsoid';
    fframe  = 'moon_PA';
    back    = BEIDOUID;
    bshape  = 'point';
    bframe  = '-';
    obsrvr  = ID_RX;
    abcorr  = 'LT+S';

   
    MOON_occultation = cspice_gfoclt( occtyp, front, fshape, fframe, ...
                          back, bshape,bframe,          ...
                          abcorr, obsrvr, step_OCC, cnfine,  ...
                          MAXWIN);
                     
   
                      
    % List the beginning and ending times in each interval
    % if result contains data.

   
    
    left_M = zeros(1,numel( MOON_occultation)/2 + 1);
    right_M = zeros(1,numel( MOON_occultation)/2 + 1);
    
    for i = 1:numel(MOON_occultation)/2

     [left_M(i), right_M(i)] = cspice_wnfetd( MOON_occultation, i );

    end
    
    % Adding to the windows the endtime of the simulation, needed to
    % compute after the last occultation window or if the window is empty
    
    left_M(end) = tspan(end);
    right_M(end) = tspan(end);
     
    
    
    % For each time step where the TX is in view of the SV, the ideal
    % range, relative velocity and angular separation with respect to the
    % SV nadir in evaluated
    
    % initializzation of the Earth and Moon occultation windows counter
    
     n = 1; % Earth occultation windows counter
     m = 1; % Moon occultation windows counter
    
     
    for r = 1 : numel(tspan) % for each time step
        
         et = tspan(r);
         
         % saving the epoch in the output structure
         
         BEIDOU.epoch(r,1) = et ;
         
         
         % If the timestep is not inside an occultation window the
         % computation are performed
         
         if et <= left_E(n) && et <= left_M(m)

             % Computing all the relevant position and state vector at the
             % current timestep 
           
             % State of the RX as seen from the TX SV
             
             [state,ltime] = cspice_spkezr (ID_RX, et,'J2000', 'XLT', BEIDOUID);

             TX2RX_pos_J2000 = state(1:3,1); % Position vector
             TX2RX_vel_J2000 = state(4:6,1); % Velocity vector
             
             % Position of Earth as seen from the TX SV ( an its nadir
             % pointing)
             
             TX2EA_pos_J2000 = cspice_spkpos ( 'EARTH', et  ,'J2000', 'LT', BEIDOUID);
             
             % Position of the SUN as seen from the TX SV 
             
             TX2SUN_pos_J2000 = cspice_spkpos ( 'SUN', et  ,'J2000', 'LT', BEIDOUID); 
             
             % Only the Sun direction is needed 

             TX2SUN_pos_J2000 = cspice_vhat(TX2SUN_pos_J2000); % vector normalizzation 
             
             % Position of Earth as seen from the RX ( an its high gain
             % antenna pointing)
             
             RX2EA_pos_J2000 = cspice_spkpos ( 'EARTH', et ,'J2000', 'LT', ID_RX);

             % We can compute the distance between the TX and RX as the norm of the position vector    
             
             d = norm(TX2RX_pos_J2000);
             
             range = ltime * c; % equivalently is the light time multiplied by the speed of light

             % Radial velocity from the SV to the RX
             
             v_radial = dot(TX2RX_vel_J2000,(TX2RX_pos_J2000/d));
             
             % To compute the off-bore angle and the azimuth of RX wrt 3D
             % antenna pattern of each SV a simple yaw-steering (YS) attitude model is
             % implemented 
             
             % For the BEIDOU constellation as no precise information were found, each Sv
             % is assumed as a BIIF/BIII' GPS SV .
            
             % For IIF and III block the x_YS axis is pointed toward the Sun direction.
                              
             % In nominal operation the z-axis is pointing the nadir
             % direcion toward Earth, the x-axis is toward the Sun
             % direction and the y-axis complete the othonormal basis
            

             k_YS = cspice_vhat(TX2EA_pos_J2000);

             j_YS = cross(TX2SUN_pos_J2000,TX2EA_pos_J2000)/norm(cross(TX2SUN_pos_J2000,TX2EA_pos_J2000));

             i_YS = cross(j_YS,k_YS);
             
             % A rotation matrix can be compute to rotate a vector in
             % the SV YS fram ( ideally equivalent to the body frame)

             R_J20002YS = [i_YS'; j_YS'; k_YS'];
             
             % The direction vector of the RX in the TX body frame can
             % be computed rotating the RX2RX position vector from the
             % J2000 frame 

             TX2RX_pos_YS = R_J20002YS * cspice_vhat(TX2RX_pos_J2000); 
             
             % The angle can be computed as a transformation into
             % spherical coordinates where:

             % Theta is the angle between the z-axis ( nadir pointing)
             % and the RX direction in the YS_frame. This is the
             % off-bore angle. It is bounded in [0, pi];

             % Phi is the angle between the +x_YS axis and the
             % projection of the RX direction into the xy plane. It is
             % consistend with the convention for the azimuth of the
             % directivity pattern and is the azimuth. It is bounded in [-pi, pi];

             [~,theta,phi] = cspice_recsph(TX2RX_pos_YS);

             TX2RX_offbore = cspice_convrt ( theta,'RADIANS', 'DEGREES');
             TX2RX_azimuth = cspice_convrt ( phi,'RADIANS', 'DEGREES');

             
             
             % The Gain pattern data has the azimuth defined in the range
             % [0 2pi]
             
             if  TX2RX_azimuth < 0
                 
                  TX2RX_azimuth = 2*pi*dpr + TX2RX_azimuth;
                  
             end
             
             % Angular separation between the Earth pointing antenna of the
             % receiver SPC and the position of the SV ( off-bore angle)

             RX2TX_offbore = cspice_convrt ( cspice_vsep(RX2EA_pos_J2000, -TX2RX_pos_J2000),'RADIANS', 'DEGREES'); 


             BEIDOU.RX_TXdistance(r,k) = d ;
             BEIDOU.RX_TXv_radial(r,k) = v_radial ;
             BEIDOU.TX_offbore(r,k) = TX2RX_offbore ;
             BEIDOU.TX_azimuth(r,k) = TX2RX_azimuth ;
             BEIDOU.RX_offbore(r,k) = RX2TX_offbore ;
             BEIDOU.visibility_flag(r,k) = 1 ;
                 
         else
             
             % else if the RX is in occultation or by Earth or Moon the SV
             % is not visible
             
             BEIDOU.visibility_flag(r,k) = 0 ;

         end
         
         % Once the time step surpass a time window the counter is increased to analyze the next one
         
         if et >= right_E(n) && n < length(right_E)
             n = n + 1;
         end
         
         if et >= right_M(m) && m < length(right_M)
             m = m + 1;
         end
         

    end
     
    % In the first step the occultation flag is always 1 even if there is
    % no visibility due to the way spice the logic check for the spice
    % windows is perfomed, so an extra check for the first row need for
    % each SV
    
    if (left_E(1) - BEIDOU.epoch(1,1)) == 0
        
        BEIDOU.visibility_flag(1,k) = 0 ;
        
    end

       
                      
end
end

%% Evaluation of the link budget as Carrier over noise ( C/N0 [dB-Hz]), doppler shift ( D_shift [Hz]) value for each SV at each time step and GDOP


% Constants and parameters

% Implementation losses 

L_B = 2; % tipical value

k_b = 1.381*1e-23; % Boltzmann constant

T_DS = 13.571; % [K] temperature of the cosmic radiation for the L5 frequency (worst case)(source ITU-R P.372-14)

T_E = 290; % [K] temperature of Earth radiation for the L1 frequency (source "Antenna pattern" and ITU-R P.372-14)

T_M = 280; % [K] temperature of Moon radiation for the L1 frequency at full Moon (worst case scenario)(source ITU-R P.372-14)

% Earth radius

R_E = cspice_bodvrd ( 'EARTH', 'RADII',3);
R_E = mean(R_E);

% Moon radius

R_M = cspice_bodvrd ( 'MOON', 'RADII',3);
R_M = mean(R_M);

% Antenna efficiency 

ANT_e = 0.75; % tipical value

% Low Noise Amplifier (LAN) noise figure 

N_f = 2; 

% Extrapolation of the gain patterns data for the GPS SV's antenna for the
% three block avaible. BIII block has a similar gain pattern to the BIIF
% block.

[BIIR_GAIN,BIIRM_GAIN,BIIF_GAIN] = GPS_GAIN_data(0);

% Definition of the grid required to interpolate the gain in 2D 

[X_grid_IIF,Y_grid_IIF] = meshgrid(0:1:360,0:1:90);
[X_grid_IIR,Y_grid_IIR] = meshgrid(0:10:360,-90:2:90);

% Evaluation of the link for each SV and time step divided by constellation

%% GPS constellations

if isequal(signalGPS,'L1')

    % L1 signal

    f = 1575.42 * 1e6; % [Hz]

    % Power of the signal in the considered bandwidth

    p_sharing = 0.25;

    % For the GPS constellation specific information of the antenna Gain and trasmitted power
    % for each generation block are avaiable

    % From GPSID -601 to - 608 --> BIIR   SV 
    % From GPSID -609 to - 615 --> BIIR-M SV
    % From GPSID -616 to - 627 --> BII-F  SV
    % From GPSID -628 to - 630 --> BIII   SV


    % For each SV 

    for k = 1:30

        % Mean emitted power of the TX for GPS SV (source "GNSS satellite
        % transmit power and its impact on orbit % determination")

        if k <= 8  % BIIR SV

            P_TX =  10*log10(p_sharing * 60); % 60W mean transmitted power in dB considering the power in the bandwidth
            G = BIIR_GAIN;
            X_grid = X_grid_IIR;
            Y_grid = Y_grid_IIR;

        elseif k <= 15 % BIIR-M SV

            P_TX = 10*log10(p_sharing * 145); % 145W mean transmitted power in dB considering the power in the bandwidth
            G = BIIRM_GAIN;
            X_grid = X_grid_IIR;
            Y_grid = Y_grid_IIR;

        elseif k <= 30 % BII-F and BIII SV

            P_TX = 10*log10(p_sharing * 240); % 240W mean transmitted power in dB considering the power in the bandwidth
            G = BIIF_GAIN;
            X_grid = X_grid_IIF;
            Y_grid = Y_grid_IIF;

        end

        % For each time step

        for r = 1 : numel(tspan)

            % If the SV is visible and with a separation less than 90 deg

            TX_offbore = GPS.TX_offbore(r,k);
            TX_azimuth = GPS.TX_azimuth(r,k);

            RX2TX_offbore = GPS.RX_offbore(r,k);

            if  GPS.visibility_flag(r,k) == 1 && TX_offbore <= 90

                % Computation of the free space path losses for C/A L1 signal

                d = GPS.RX_TXdistance(r,k); % distance of the RX wrt the TX

                A_d = FSP_Loss(d,f,c); % free space path losses [dB]


                % Evaluation of the RX gain, function of the off-boresigth angle

                 G_RX = RX_GAIN(RX2TX_offbore);

                % Evaluation of the TX gain, function of the off-boresigth
                % angle and the SV attitude

                % Block BIIR and BIIR-M

                if  k < 16

                    if TX_azimuth > pi*dpr 

                        G_TX = interp2(X_grid,Y_grid,G,TX_azimuth,-TX_offbore);

                    else

                        G_TX = interp2(X_grid,Y_grid,G,TX_azimuth, TX_offbore);

                    end

                else

                % Block BIIF and BIII

                G_TX = interp2(X_grid,Y_grid,G,TX_azimuth, TX_offbore);

                end

                % For an isotropic antenna the T_a is a weighted mean of the temperature
                % for a unitary sphere integration


                et = GPS.epoch(r);

                pos = cspice_spkpos ( 'EARTH', et,'J2000', 'LT', ID_RX);
                d_E = norm(pos);

                pos = cspice_spkpos ( 'MOON', et,'J2000', 'LT', ID_RX);
                d_M = norm(pos);

                % Earth solid angle seen from the RX 

                theta = atan(R_E/d_E);
                OM_E = 2*pi * ( 1 - cos(theta));

                % Moon solid angle seen from the RX 

                theta = atan(R_M/d_M);
                OM_M = 2*pi * ( 1 - cos(theta));

                % Antenna temperature 

                T_ant = (1/(4*pi)) * (OM_E * (25 * T_E) + OM_M * T_M + (4*pi - (OM_E + OM_M)) * T_DS);

                % The antennta noise temperature including thermal inefficiencies

                T_ant =  T_ant   + 290 * (1/ANT_e -1);

                % Low Noise Amplifier (LNA) noise temperature

                T_amp = 290*(10^(N_f/10) - 1); % For an assumed 290 K  temperature 

                % Total system temperature

                T_sys = T_ant + T_amp;

                % Carrier signal [dBW]

                C = P_TX - A_d + G_TX + G_RX - L_B;


                % Noise [dBW - Hz]

                N_0 = 10*log10( k_b * ( T_sys));

                C_N0 = C - N_0; % dB - Hz

                GPS.C_N0(r,k) = C_N0;

                GPS.acquired(r,k) = C_N0 > aq_treshold;


                % Doppler shift computation

                GPS.D_shift(r,k) = doppler_shift(f,GPS.RX_TXv_radial(r,k),c);

            end


        end

    end



elseif isequal(signalGPS,'L5')
    

    % L5 signal 

    f = 1176.45 * 1e6; % [Hz]

    % Power of the signal in the considered bandwidth

    p_sharing = 0.5;

    % For the GPS constellation specific information of the antenna Gain and trasmitted power
    % for each generation block are avaiable

    % From GPSID -601 to - 608 --> BIIR   SV 
    % From GPSID -609 to - 615 --> BIIR-M SV
    % From GPSID -616 to - 627 --> BII-F  SV
    % From GPSID -628 to - 630 --> BIII   SV


    % Only B-IIF and B-III SVs broacast the L5 signal

    for k = 16:30

        % Mean emitted power of the TX for GPS SV (source "GNSS satellite
        % transmit power and its impact on orbit % determination")

        P_TX = 10*log10(p_sharing * 240); % 240W mean transmitted power in dB considering the power in the bandwidth 
        
        G =  BIIF_GAIN; 
        X_grid = X_grid_IIF;
        Y_grid = Y_grid_IIF;

        % For each time step

        for r = 1 : numel(tspan)

            % If the SV is visible and with a separation less than 90 deg

            TX_offbore = GPS.TX_offbore(r,k);
            TX_azimuth = GPS.TX_azimuth(r,k);

            RX2TX_offbore = GPS.RX_offbore(r,k);

            if  GPS.visibility_flag(r,k) == 1 && TX_offbore <= 90

                % Computation of the free space path losses for C/A L1 signal

                d = GPS.RX_TXdistance(r,k); % distance of the RX wrt the TX

                A_d = FSP_Loss(d,f,c); % free space path losses [dB]


                % Evaluation of the RX gain, function of the off-boresigth angle

                 G_RX = RX_GAIN(RX2TX_offbore);

                % Evaluation of the TX gain, function of the off-boresigth
                % angle and the SV attitude

                % Block BIIF and BIII

                G_TX = interp2(X_grid,Y_grid,G,TX_azimuth, TX_offbore/(1.1));

                % For an isotropic antenna the T_a is a weighted mean of the temperature
                % for a unitary sphere integration


                et = GPS.epoch(r);

                pos = cspice_spkpos ( 'EARTH', et,'J2000', 'LT', ID_RX);
                d_E = norm(pos);

                pos = cspice_spkpos ( 'MOON', et,'J2000', 'LT', ID_RX);
                d_M = norm(pos);

                % Earth solid angle seen from the RX 

                theta = atan(R_E/d_E);
                OM_E = 2*pi * ( 1 - cos(theta));

                % Moon solid angle seen from the RX 

                theta = atan(R_M/d_M);
                OM_M = 2*pi * ( 1 - cos(theta));

                % Antenna temperature 

                T_ant = (1/(4*pi)) * (OM_E * (25 * T_E) + OM_M * T_M + (4*pi - (OM_E + OM_M)) * T_DS);

                % The antennta noise temperature including thermal inefficiencies

                T_ant =  T_ant   + 290 * (1/ANT_e -1);

                % Low Noise Amplifier (LNA) noise temperature

                T_amp = 290*(10^(N_f/10) - 1); % For an assumed 290 K  temperature 

                % Total system temperature

                T_sys = T_ant + T_amp;

                % Carrier signal [dBW]

                C = P_TX - A_d + G_TX + G_RX - L_B;


                % Noise [dBW - Hz]

                N_0 = 10*log10( k_b * ( T_sys));

                C_N0 = C - N_0; % dB - Hz

                GPS.C_N0(r,k) = C_N0;

                GPS.acquired(r,k) = C_N0 > aq_treshold;


                % Doppler shift computation

                GPS.D_shift(r,k) = doppler_shift(f,GPS.RX_TXv_radial(r,k),c);

            end


        end

    end

end

%% GDOP

% We can compute under a certain set of hypothesis the GDOP as in "GPS
% principles and applications" (Kaplan)

% At each time step 

for r = 1 : numel(GPS.epoch)

et = GPS.epoch(r);

% Number of acquired and tracked SV at the time considered

SV_tracked = sum(GPS.acquired(r,:));

GPS.n_tracked(r,1) = SV_tracked;


H = zeros(SV_tracked,4);

n = 1;

% For each SV of the constellation

    for k = 1 : 30

        % if the SV is being tracked

        if GPS.acquired(r,k) > 0 && SV_tracked > 3

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

        GPS.GDOP(r,1) = GDOP;

    else

        GPS.GDOP(r,1) = nan;

    end
end



%% GALILEO constellations

if isequal(signalGALILEO,'E1')

    % E1 signal ( Carrier frequency equal to L1)

    f = 1575.42 * 1e6; % [Hz]

    % Power of the signal in the considered bandwidth

    p_sharing = 0.25;

    % For the GALILEO constellation partial and not official data both for the
    % gain and the transmitted power are present 

    G = BIIF_GAIN; % Using the GPS IIF/III-A gain pattern as reference for Galileo constellation

    X_grid = X_grid_IIF;
    Y_grid = Y_grid_IIF;

    % For each SV 

    for k = 1:24

        % For each time step

        % Mean emitted power of the TX for GALILEO GSAT SV (source "GNSS satellite
        % transmit power and its impact on orbit % determination")

        if k <= 2

            P_TX = 10*log10(p_sharing * 135); % 135W mean transmitted power in dB considering the power in the bandwidth

        elseif k == 3

            P_TX = 10*log10(p_sharing * 95); % 95W mean transmitted power in dB considering the power in the bandwidth

        else

            P_TX = 10*log10(p_sharing * 265); % 255W mean transmitted power in dB considering the power in the bandwidth

        end

        for r = 1 : numel(tspan)

             % If the SV is visible and with a separation less than 90 deg

            TX_offbore = GALILEO.TX_offbore(r,k);
            TX_azimuth = GALILEO.TX_azimuth(r,k);

            RX2TX_offbore = GALILEO.RX_offbore(r,k);

            if  GALILEO.visibility_flag(r,k) == 1 && TX_offbore <= 90

                % Computation of the free space path losses for C/A L1 signal

                d = GALILEO.RX_TXdistance(r,k); % distance of the RX wrt the TX

                A_d = FSP_Loss(d,f,c); % free space path losses [dB]


                % Evaluation of the RX gain, function of the off-boresigth angle

                 G_RX = RX_GAIN(RX2TX_offbore);

                % Evaluation of the TX gain, function of the off-boresigth
                % angle and the SV attitude


                G_TX = interp2(X_grid,Y_grid,G,TX_azimuth, TX_offbore);


                % For an isotropic antenna the T_a is a weighted mean of the temperature
                % for a unitary sphere integration


                et = GALILEO.epoch(r);

                pos = cspice_spkpos ( 'EARTH', et,'J2000', 'LT', ID_RX);
                d_E = norm(pos);

                pos = cspice_spkpos ( 'MOON', et,'J2000', 'LT', ID_RX);
                d_M = norm(pos);

                % Earth solid angle seen from the RX 

                theta = atan(R_E/d_E);
                OM_E = 2*pi * ( 1 - cos(theta));

                % Moon solid angle seen from the RX 

                theta = atan(R_M/d_M);
                OM_M = 2*pi * ( 1 - cos(theta));

                % Antenna temperature 

                T_ant = (1/(4*pi)) * (OM_E * (25 * T_E) + OM_M * T_M + (4*pi - (OM_E + OM_M)) * T_DS);

                % The antennta noise temperature including thermal inefficiencies

                T_ant =  T_ant   + 290 * (1/ANT_e -1);

                % Low Noise Amplifier (LNA) noise temperature

                T_amp = 290*(10^(N_f/10) - 1); % For an assumed 290 K  temperature 

                % Total system temperature

                T_sys = T_ant + T_amp;

                % Carrier signal [dBW]

                C = P_TX - A_d + G_TX + G_RX - L_B;


                % Noise [dBW - Hz]

                N_0 = 10*log10( k_b * ( T_sys));

                C_N0 = C - N_0; % dB - Hz

                GALILEO.C_N0(r,k) = C_N0;

                GALILEO.acquired(r,k) = C_N0 > aq_treshold;


                % Doppler shift computation

                GALILEO.D_shift(r,k) = doppler_shift(f,GALILEO.RX_TXv_radial(r,k),c);

            end


        end

    end


elseif isequal(signalGALILEO,'E5a')

    % E5a signal ( Carrier frequency equal to L1)

    f = 1176.45 * 1e6; % [Hz]

    % Power of the signal in the considered bandwidth

    p_sharing = 0.5;

    % For the GALILEO constellation partial and not official data both for the
    % gain and the transmitted power are present 

    G = BIIF_GAIN; % Using the GPS IIF/III-A gain pattern as reference for Galileo constellation 

    X_grid = X_grid_IIF;
    Y_grid = Y_grid_IIF;

    % For each SV 

    for k = 1:24

        % For each time step

        % Mean emitted power of the TX for GALILEO GSAT SV (source "GNSS satellite
        % transmit power and its impact on orbit % determination")

        if k <= 2

            P_TX = 10*log10(p_sharing * 135);  % 135W mean transmitted power in dB considering the power in the bandwidth 
            
        elseif k == 3

            P_TX = 10*log10(p_sharing * 95);   % 95W mean transmitted power in dB considering the power in the bandwidth 
            
        else

            P_TX = 10*log10(p_sharing * 265); % 265W mean transmitted power in dB considering the power in the bandwidth 
            
        end

        for r = 1 : numel(tspan)

             % If the SV is visible and with a separation less than 90 deg

            TX_offbore = GALILEO.TX_offbore(r,k);
            TX_azimuth = GALILEO.TX_azimuth(r,k);

            RX2TX_offbore = GALILEO.RX_offbore(r,k);

            if  GALILEO.visibility_flag(r,k) == 1 && TX_offbore <= 90

                % Computation of the free space path losses for C/A L1 signal

                d = GALILEO.RX_TXdistance(r,k); % distance of the RX wrt the TX

                A_d = FSP_Loss(d,f,c); % free space path losses [dB]


                % Evaluation of the RX gain, function of the off-boresigth angle

                G_RX = RX_GAIN(RX2TX_offbore);

                % Evaluation of the TX gain, function of the off-boresigth
                % angle and the SV attitude


                G_TX = interp2(X_grid,Y_grid,G,TX_azimuth, TX_offbore/(1.1));


                % For an isotropic antenna the T_a is a weighted mean of the temperature
                % for a unitary sphere integration


                et = GALILEO.epoch(r);

                pos = cspice_spkpos ( 'EARTH', et,'J2000', 'LT', ID_RX);
                d_E = norm(pos);

                pos = cspice_spkpos ( 'MOON', et,'J2000', 'LT', ID_RX);
                d_M = norm(pos);

                % Earth solid angle seen from the RX 

                theta = atan(R_E/d_E);
                OM_E = 2*pi * ( 1 - cos(theta));

                % Moon solid angle seen from the RX 

                theta = atan(R_M/d_M);
                OM_M = 2*pi * ( 1 - cos(theta));

                % Antenna temperature 

                T_ant = (1/(4*pi)) * (OM_E * (25 * T_E) + OM_M * T_M + (4*pi - (OM_E + OM_M)) * T_DS);

                % The antennta noise temperature including thermal inefficiencies

                T_ant =  T_ant   + 290 * (1/ANT_e -1);

                % Low Noise Amplifier (LNA) noise temperature

                T_amp = 290*(10^(N_f/10) - 1); % For an assumed 290 K  temperature 

                % Total system temperature

                T_sys = T_ant + T_amp;

                % Carrier signal [dBW]

                C = P_TX - A_d + G_TX + G_RX - L_B;


                % Noise [dBW - Hz]

                N_0 = 10*log10( k_b * ( T_sys));

                C_N0 = C - N_0; % dB - Hz

                GALILEO.C_N0(r,k) = C_N0;

                GALILEO.acquired(r,k) = C_N0 > aq_treshold;


                % Doppler shift computation

                GALILEO.D_shift(r,k) = doppler_shift(f,GALILEO.RX_TXv_radial(r,k),c);

            end


        end

    end
    
end
%% GDOP

% We can compute under a certain set of hypothesis the GDOP as in "GPS
% principles and applications" (Kaplan)

% At each time step 

for r = 1 : numel(GALILEO.epoch)
    
et = GALILEO.epoch(r);

% Number of acquired and tracked SV at the time considered

SV_tracked = sum(GALILEO.acquired(r,:));

GALILEO.n_tracked(r,1) = SV_tracked;

H = zeros(SV_tracked,4);

n = 1;

% For each SV of the constellation

    for k = 1 : 24

        % if the SV is being tracked

        if GALILEO.acquired(r,k) > 0 && SV_tracked > 3

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

        GALILEO.GDOP(r,1) = GDOP;
        
   else
       
        GALILEO.GDOP(r,1) =  nan;
    
   end
end



%% GLONASS constellations

if GLO_flag == 1
    
% G1/L1 signal

f = 1602 * 1e6; % [Hz]

% For the GLONASS COSMOS-M SV data about the power transmitted in the band
% are available

% Power of the signal in the considered bandwidth

p_sharing = 0.25;

P_TX = 10*log10(p_sharing*50); % 50W mean transmitted power in dB 

% For the GLONASS constellation data for the
% gain pattern was not found

G = BIIF_GAIN; % Using the GPS IIF/III-A gain pattern as reference also for GLONASS constellation

% For each SV 

for k = 1:24

    % For each time step
    
    for r = 1 : numel(tspan)
       
        % If the SV is visible and with a separation less than 90 deg
        
        TX_offbore = GLONASS.TX_offbore(r,k);
        TX_azimuth = GLONASS.TX_azimuth(r,k);
        
        RX2TX_offbore = GLONASS.RX_offbore(r,k);
        
        if  GLONASS.visibility_flag(r,k) == 1 && TX_offbore <= 90
        
            % Computation of the free space path losses for C/A L1 signal
               
            d = GLONASS.RX_TXdistance(r,k); % distance of the RX wrt the TX
       
            A_d = FSP_Loss(d,f,c); % free space path losses [dB]
            
            
            % Evaluation of the RX gain, function of the off-boresigth angle
                       
             G_RX = RX_GAIN(RX2TX_offbore);

            % Evaluation of the TX gain, function of the off-boresigth
            % angle and the SV attitude

            
            G_TX = interp2(X_grid,Y_grid,G,TX_azimuth, TX_offbore);


            % For an isotropic antenna the T_a is a weighted mean of the temperature
            % for a unitary sphere integration


            et = GLONASS.epoch(r);

            pos = cspice_spkpos ( 'EARTH', et,'J2000', 'LT', ID_RX);
            d_E = norm(pos);

            pos = cspice_spkpos ( 'MOON', et,'J2000', 'LT', ID_RX);
            d_M = norm(pos);

            % Earth solid angle seen from the RX 

            theta = atan(R_E/d_E);
            OM_E = 2*pi * ( 1 - cos(theta));

            % Moon solid angle seen from the RX 

            theta = atan(R_M/d_M);
            OM_M = 2*pi * ( 1 - cos(theta));

            % Antenna temperature 

            T_ant = (1/(4*pi)) * (OM_E * (25 * T_E) + OM_M * T_M + (4*pi - (OM_E + OM_M)) * T_DS);

            % The antennta noise temperature including thermal inefficiencies

            T_ant =  T_ant   + 290 * (1/ANT_e -1);

            % Low Noise Amplifier (LNA) noise temperature

            T_amp = 290*(10^(N_f/10) - 1); % For an assumed 290 K  temperature 

            % Total system temperature

            T_sys = T_ant + T_amp;

            % Carrier signal [dBW]

            C = P_TX - A_d + G_TX + G_RX - L_B;


            % Noise [dBW - Hz]

            N_0 = 10*log10( k_b * ( T_sys));

            C_N0 = C - N_0; % dB - Hz

            GLONASS.C_N0(r,k) = C_N0;
            
            GLONASS.acquired(r,k) = C_N0 > aq_treshold;
            
           
            % Doppler shift computation
            
            GLONASS.D_shift(r,k) = doppler_shift(f,GLONASS.RX_TXv_radial(r,k),c);
            
        end

        
    end
    
end


%% GDOP

% We can compute under a certain set of hypothesis the GDOP as in "GPS
% principles and applications" (Kaplan)

% At each time step 

for r = 1 : numel(GLONASS.epoch)
    
et = GLONASS.epoch(r);

% Number of acquired and tracked SV at the time considered

SV_tracked = sum(GLONASS.acquired(r,:));

GLONASS.n_tracked(r,1) = SV_tracked;


H = zeros(SV_tracked,4);

n = 1;

% For each SV of the constellation

    for k = 1 : 24

        % if the SV is being tracked

        if GLONASS.acquired(r,k) > 0 && SV_tracked > 3

            GLONASSID = -(700 + k);
            GLONASSID = char(string(GLONASSID));  

            % Position vector of the tracked SV as seen from the TX

            a_n = cspice_spkpos ( GLONASSID, et,'J2000', 'XLT + S', ID_RX);

            % Normalizzation of the position vector

            a_n = a_n' / norm(a_n);

            H(n,:) = [a_n 1];

            n = n + 1;

        end
        
    end
    
    if SV_tracked > 3
        
        D = inv(H' * H);

        GDOP = sqrt(trace(D));

        GLONASS.GDOP(r,1) = GDOP;
        
    else
        
        GLONASS.GDOP(r,1) =  nan;
    
    end
end
end

%% BEIDOU-3 constellations

if BEI_flag == 1

% B1 signal

f = 1575.42 * 1e6; % [Hz]

% Power of the signal in the considered bandwidth

p_sharing = 0.25;

% For the BEIDOU III constellation data for the
% gain pattern was not found

G = BIIF_GAIN; % Using the GPS IIF/III-A gain pattern as reference also for BEIDOU III constellation

% For each SV 

for k = 1:30

    % Mean emitted power of the TX for BEUDOU II SV used as values for BEIDOU III
    %(source "GNSS satellite transmit power and its impact on orbit determination")
    
    % SV in MEO
    
    if k <= 24

        P_TX = 10*log10(p_sharing * 130); % 130W mean transmitted power in dB considering the power in the bandwidth
    
    else 
        
        % SV IN IGSO and GEO
        
        P_TX = 10*log10(p_sharing * 185); % 1855W mean transmitted power in dB considering the power in the bandwidth
        
    end
    
    % For each time step
    
    for r = 1 : numel(tspan)
       
        % If the SV is visible and with a separation less than 90 deg
        
        TX_offbore = BEIDOU.TX_offbore(r,k);
        TX_azimuth = BEIDOU.TX_azimuth(r,k);
        
        RX2TX_offbore = BEIDOU.RX_offbore(r,k);
        
        if  BEIDOU.visibility_flag(r,k) == 1 && TX_offbore <= 90
        
            % Computation of the free space path losses for C/A L1 signal
               
            d = BEIDOU.RX_TXdistance(r,k); % distance of the RX wrt the TX
       
            A_d = FSP_Loss(d,f,c); % free space path losses [dB]
            
            
            % Evaluation of the RX gain, function of the off-boresigth angle
                       
             G_RX = RX_GAIN(RX2TX_offbore);

            % Evaluation of the TX gain, function of the off-boresigth
            % angle and the SV attitude

            
            G_TX = interp2(X_grid,Y_grid,G,TX_azimuth, TX_offbore);


            % For an isotropic antenna the T_a is a weighted mean of the temperature
            % for a unitary sphere integration


            et = BEIDOU.epoch(r);

            pos = cspice_spkpos ( 'EARTH', et,'J2000', 'LT', ID_RX);
            d_E = norm(pos);

            pos = cspice_spkpos ( 'MOON', et,'J2000', 'LT', ID_RX);
            d_M = norm(pos);

            % Earth solid angle seen from the RX 

            theta = atan(R_E/d_E);
            OM_E = 2*pi * ( 1 - cos(theta));

            % Moon solid angle seen from the RX 

            theta = atan(R_M/d_M);
            OM_M = 2*pi * ( 1 - cos(theta));

            % Antenna temperature 

            T_ant = (1/(4*pi)) * (OM_E * (25 * T_E) + OM_M * T_M + (4*pi - (OM_E + OM_M)) * T_DS);

            % The antennta noise temperature including thermal inefficiencies

            T_ant =  T_ant   + 290 * (1/ANT_e -1);

            % Low Noise Amplifier (LNA) noise temperature

            T_amp = 290*(10^(N_f/10) - 1); % For an assumed 290 K  temperature 

            % Total system temperature

            T_sys = T_ant + T_amp;

            % Carrier signal [dBW]

            C = P_TX - A_d + G_TX + G_RX - L_B;


            % Noise [dBW - Hz]

            N_0 = 10*log10( k_b * ( T_sys));

            C_N0 = C - N_0; % dB - Hz

            BEIDOU.C_N0(r,k) = C_N0;
            
            BEIDOU.acquired(r,k) = C_N0 > aq_treshold;
            
           
            % Doppler shift computation
            
            BEIDOU.D_shift(r,k) = doppler_shift(f,BEIDOU.RX_TXv_radial(r,k),c);
            
        end

        
    end
    
end


%% GDOP

% We can compute under a certain set of hypothesis the GDOP as in "GPS
% principles and applications" (Kaplan)

% At each time step 

for r = 1 : numel(BEIDOU.epoch)
    
et = BEIDOU.epoch(r);

% Number of acquired and tracked SV at the time considered

SV_tracked = sum(BEIDOU.acquired(r,:));

BEIDOU.n_tracked(r,1) = SV_tracked;


H = zeros(SV_tracked,4);

n = 1;

% For each SV of the constellation

    for k = 1 : 30

        % if the SV is being tracked

        if BEIDOU.acquired(r,k) > 0 && SV_tracked > 3

            BEIDOUID = -(750 + k);
            BEIDOUID = char(string(BEIDOUID));  

            % Position vector of the tracked SV as seen from the TX

            a_n = cspice_spkpos ( BEIDOUID, et,'J2000', 'XLT + S', ID_RX);

            % Normalizzation of the position vector

            a_n = a_n' / norm(a_n);

            H(n,:) = [a_n 1];

            n = n + 1;

        end
        
    end
    
    if SV_tracked > 3
        
        D = inv(H' * H);

        GDOP = sqrt(trace(D));

        BEIDOU.GDOP(r,1) = GDOP;
        
    else
        
        BEIDOU.GDOP(r,1) =  nan;
    
    end
end
end

%% GDOP for multiple constellations

MULTI_GNSS = struct('GPSandGAL_GDOP',[],'ALL_GNSS_GDOP',[]);


%% GPS + GALILEO

% At each time step 

for r = 1 : numel(GPS.epoch)
    
et = GPS.epoch(r);

% Number of acquired and tracked SV at the timeconsidered

SV_tracked_GPS = sum(GPS.acquired(r,:));
SV_tracked_GALILEO = sum(GALILEO.acquired(r,:));

SV_tracked = (SV_tracked_GPS + SV_tracked_GALILEO);

H = zeros(SV_tracked,4);

n = 1;

% For each SV of the GPS constellation

    for k = 1 : 30

        % if the SV is being tracked

        if GPS.acquired(r,k) > 0 && SV_tracked > 3

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

        if GALILEO.acquired(r,k) > 0 && SV_tracked > 3

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

        MULTI_GNSS.GPSandGAL_GDOP(r,1) = GDOP;
        
     else
        
        MULTI_GNSS.GPSandGAL_GDOP(r,1) =  nan;
    
    end
end

%% ALL GNSS

if GLO_flag == 1 && BEI_flag == 1

% At each time step 

for r = 1 : numel(GPS.epoch)
    
et = GPS.epoch(r);

% Number of acquired and tracked SV at the timeconsidered

SV_tracked_GPS = sum(GPS.acquired(r,:));
SV_tracked_GALILEO = sum(GALILEO.acquired(r,:));
SV_tracked_GLONASS = sum(GLONASS.acquired(r,:));
SV_tracked_BEIDOU = sum(BEIDOU.acquired(r,:));

SV_tracked = (SV_tracked_GPS + SV_tracked_GALILEO + SV_tracked_GLONASS + SV_tracked_BEIDOU);

H = zeros(SV_tracked,4);

n = 1;

% For each SV of the GPS constellation

    for k = 1 : 30

        % if the SV is being tracked

        if GPS.acquired(r,k) > 0 && SV_tracked > 3

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

        if GALILEO.acquired(r,k) > 0 && SV_tracked > 3

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
        
    % For each SV of the GLONASS constellation

    for k = 1 : 24

        % if the SV is being tracked

        if GLONASS.acquired(r,k) > 0 && SV_tracked > 3

            GLONASSID = -(700 + k);
            GLONASSID = char(string(GLONASSID));  

            % Position vector of the tracked SV as seen from the TX

            a_n = cspice_spkpos ( GLONASSID, et,'J2000', 'XLT + S', ID_RX);

            % Normalizzation of the position vector

            a_n = a_n' / norm(a_n);

            H(n,:) = [a_n 1];

            n = n + 1;

        end
        
    end
    
    % For each SV of the BEIDOU-3 constellation

    for k = 1 : 30

        % if the SV is being tracked

        if BEIDOU.acquired(r,k) > 0 && SV_tracked > 3

            BEIDOUID = -(750 + k);
            BEIDOUID = char(string(BEIDOUID));  

            % Position vector of the tracked SV as seen from the TX

            a_n = cspice_spkpos ( BEIDOUID, et,'J2000', 'XLT + S', ID_RX);

            % Normalizzation of the position vector

            a_n = a_n' / norm(a_n);

            H(n,:) = [a_n 1];

            n = n + 1;

        end
        
    end
    
    if SV_tracked > 3
        
        D = inv(H' * H);

        GDOP = sqrt(trace(D));

        MULTI_GNSS.ALL_GNSS_GDOP(r,1) = GDOP;
        
     else
        
        MULTI_GNSS.ALL_GNSS_GDOP(r,1) =  nan;
    
    end
end
end 


%% LNNSS

% Similiar visibility and link budget computation are performed for a
% theoretical LNNSS constellation of 6 and 8 SV called "LNNSS" in
% two orbital plane inclined in order to be orthogonal to the mean orbital plane of the MOON in the
% time span considered (5 years starting 10-DEC-2020) and orthogonal to
% each other and the antenna array pointing to the Moon center instead of
% the nadir direction

LNNSS_6SV = struct('epoch',[],'RX_TXdistance',[],'RX_TXv_radial',[],'RX_offbore',[],'TX_offbore',[],'TX_azimuth',[],'visibility_flag',[]); % Struct containing the results

LNNSS_8SV = struct('epoch',[],'RX_TXdistance',[],'RX_TXv_radial',[],'RX_offbore',[],'TX_offbore',[],'TX_azimuth',[],'visibility_flag',[]); % Struct containing the results

    if LNNSS_flag == 1

    if isequal(LNNSS,'H')


        % Heavy 6 SV constellation
        cspice_furnsh('kernels/spk/6SV_H_LNNSS.bsp');

    elseif isequal(LNNSS,'L')


        % Light 6 SV constellation
        cspice_furnsh('kernels/spk/6SV_L_LNNSS.bsp');

    end


        %% LNNSS constellation

        for k = 1:6

            LNNSSID = -(800 + k);
            LNNSSID = char(string(LNNSSID));  


            % Time windows of Earth occultation

            occtyp  = 'any';
            front   = 'earth';
            fshape  = 'ellipsoid';
            fframe  = 'ITRF93';
            back    =  LNNSSID;
            bshape  = 'point';
            bframe  = '-';
            obsrvr  = ID_RX;
            abcorr  = 'LT+S';


            EARTH_occultation = cspice_gfoclt( occtyp, front, fshape, fframe, ...
                                  back, bshape,bframe,          ...
                                  abcorr, obsrvr, step_OCC, cnfine,  ...
                                  MAXWIN);



            % List the beginning and ending times in each interval
            % if result contains data.


            left_E = zeros(1,numel( EARTH_occultation)/2 + 1);
            right_E = zeros(1,numel( EARTH_occultation)/2 + 1);

            for i = 1:numel(EARTH_occultation)/2

             [left_E(i), right_E(i)] = cspice_wnfetd( EARTH_occultation, i );

            end


            % Adding to the windows the endtime of the simulation, needed to
            % compute after the last occultation window or if the window is empty

            left_E(end) = tspan(end);
            right_E(end) = tspan(end);


            % Time windows of Moon occultation

            occtyp  = 'any';
            front   = 'moon';
            fshape  = 'ellipsoid';
            fframe  = 'moon_PA';
            back    =  LNNSSID;
            bshape  = 'point';
            bframe  = '-';
            obsrvr  = ID_RX;
            abcorr  = 'LT+S';


            MOON_occultation = cspice_gfoclt( occtyp, front, fshape, fframe, ...
                                  back, bshape,bframe,          ...
                                  abcorr, obsrvr, step_OCC, cnfine,  ...
                                  MAXWIN);



            % List the beginning and ending times in each interval
            % if result contains data.



            left_M = zeros(1,numel( MOON_occultation)/2 + 1);
            right_M = zeros(1,numel( MOON_occultation)/2 + 1);

            for i = 1:numel(MOON_occultation)/2

             [left_M(i), right_M(i)] = cspice_wnfetd( MOON_occultation, i );

            end

            % Adding to the windows the endtime of the simulation, needed to
            % compute after the last occultation window or if the window is empty

            left_M(end) = tspan(end);
            right_M(end) = tspan(end);



            % For each time step where the TX is in view of the SV, the ideal
            % range, relative velocity and angular separation with respect to the
            % SV nadir in evaluated

            % initializzation of the Earth and Moon occultation windows counter

             n = 1; % Earth occultation windows counter
             m = 1; % Moon occultation windows counter


            for r = 1 : numel(tspan) % for each time step

                 et = tspan(r);

                 % saving the epoch in the output structure

                  LNNSS_6SV.epoch(r,1) = et ;


                 % If the timestep is not inside an occultation window the
                 % computation are performed

                 if et <= left_E(n) && et <= left_M(m)

                     % Computing all the relevant position and state vector at the
                     % current timestep 

                     % State of the RX as seen from the TX SV

                     [state,ltime] = cspice_spkezr (ID_RX, et,'J2000', 'XLT',  LNNSSID);

                     TX2RX_pos_J2000 = state(1:3,1); % Position vector
                     TX2RX_vel_J2000 = state(4:6,1); % Velocity vector

                     % Position of MOON as seen from the TX SV ( an its array pointing)

                     TX2MO_pos_J2000 = cspice_spkpos ( 'MOON', et  ,'J2000', 'LT',  LNNSSID);

                     % Position of the SUN as seen from the TX SV 

                     TX2SUN_pos_J2000 = cspice_spkpos ( 'SUN', et  ,'J2000', 'LT',  LNNSSID); 

                     % Only the Sun direction is needed 

                     TX2SUN_pos_J2000 = cspice_vhat(TX2SUN_pos_J2000); % vector normalizzation 

                     % Position of Earth as seen from the RX ( an its high gain
                     % antenna pointing)

                     RX2EA_pos_J2000 = cspice_spkpos ( 'EARTH', et ,'J2000', 'LT', ID_RX);

                     % We can compute the distance between the TX and RX as the norm of the position vector    

                     d = norm(TX2RX_pos_J2000);

                     range = ltime * c; % equivalently is the light time multiplied by the speed of light

                     % Radial velocity from the SV to the RX

                     v_radial = dot(TX2RX_vel_J2000,(TX2RX_pos_J2000/d));

                     % To compute the off-bore angle and the azimuth of RX wrt 3D
                     % antenna pattern of each SV a simple yaw-steering (YS) attitude model is
                     % implemented 

                     % For the LNNSS constellation is in similarity to the
                     % GALILEO constellation only the SVs are poiting the array
                     % toward the Moon center.


                     % In nominal operation the z-axis is pointing the Moon center,
                     % the x-axis is toward the Sun direction and the y-axis complete the orthonormal basis


                     k_YS = cspice_vhat(TX2MO_pos_J2000);

                     j_YS = cross(TX2SUN_pos_J2000,TX2MO_pos_J2000)/norm(cross(TX2SUN_pos_J2000,TX2MO_pos_J2000));

                     i_YS = cross(j_YS,k_YS);

                     % A rotation matrix can be compute to rotate a vector in
                     % the SV YS fram ( ideally equivalent to the body frame)

                     R_J20002YS = [i_YS'; j_YS'; k_YS'];

                     % The direction vector of the RX in the TX body frame can
                     % be computed rotating the RX2RX position vector from the
                     % J2000 frame 

                     TX2RX_pos_YS = R_J20002YS * cspice_vhat(TX2RX_pos_J2000); 

                     % The angle can be computed as a transformation into
                     % spherical coordinates where:

                     % Theta is the angle between the z-axis ( nadir pointing)
                     % and the RX direction in the YS_frame. This is the
                     % off-bore angle. It is bounded in [0, pi];

                     % Phi is the angle between the +x_YS axis and the
                     % projection of the RX direction into the xy plane. It is
                     % consistend with the convention for the azimuth of the
                     % directivity pattern and is the azimuth. It is bounded in [-pi, pi];

                     [~,theta,phi] = cspice_recsph(TX2RX_pos_YS);

                     TX2RX_offbore = cspice_convrt ( theta,'RADIANS', 'DEGREES');
                     TX2RX_azimuth = cspice_convrt ( phi,'RADIANS', 'DEGREES');



                     % The Gain pattern data has the azimuth defined in the range
                     % [0 2pi]

                     if  TX2RX_azimuth < 0

                          TX2RX_azimuth = 2*pi*dpr + TX2RX_azimuth;

                     end

                     % Angular separation between the Earth pointing antenna of the
                     % receiver SPC and the position of the SV ( off-bore angle)

                     RX2TX_offbore = cspice_convrt ( cspice_vsep(RX2EA_pos_J2000, -TX2RX_pos_J2000),'RADIANS', 'DEGREES'); 


                     LNNSS_6SV.RX_TXdistance(r,k) = d ;
                     LNNSS_6SV.RX_TXv_radial(r,k) = v_radial ;
                     LNNSS_6SV.TX_offbore(r,k) = TX2RX_offbore ;
                     LNNSS_6SV.TX_azimuth(r,k) = TX2RX_azimuth ;
                     LNNSS_6SV.RX_offbore(r,k) = RX2TX_offbore ;
                     LNNSS_6SV.visibility_flag(r,k) = 1 ;


                 else

                     % else if the RX is in occultation or by Earth or Moon the SV
                     % is not visible

                     LNNSS_6SV.visibility_flag(r,k) = 0 ;

                 end

                 % Once the time step surpass a time window the counter is increased to analyze the next one

                 if et >= right_E(n) && n < length(right_E)
                     n = n + 1;
                 end

                 if et >= right_M(m) && m < length(right_M)
                     m = m + 1;
                 end


            end

            % In the first step the occultation flag is always 1 even if there is
            % no visibility due to the way the logic check for the spice
            % windows is perfomed, so an extra check for the first row is needed for
            % each SV

            if (left_E(1) - LNNSS_6SV.epoch(1,1)) == 0

                LNNSS_6SV.visibility_flag(1,k) = 0 ;

            end



        end


        %% LNNSS link budget

        % E5a signal (Carrier frequency equal to L5)

        f = 1176.45 * 1e6; % [Hz]

        % Power of the signal in the considered bandwidth

        p_sharing = 0.5;

        % For the GALILEO constellation partial and not official data both for the
        % gain and the transmitted power are present 

        G = BIIF_GAIN; % Using the GPS IIF/III-A gain pattern as reference for LNNSS constellation 

        X_grid = X_grid_IIF;
        Y_grid = Y_grid_IIF;

        % For each SV 

        for k = 1:6

            % For each time step

            % Mean emitted power of the TX for GALILEO GSAT SV (source "GNSS satellite
            % transmit power and its impact on orbit % determination")

            P_TX = 10*log10(p_sharing * 265/4); %  66.25W mean transmitted power in dB considering the power in the bandwidth 


            for r = 1 : numel(tspan)

                % If the SV is visible and with a separation less than 90 deg

                TX_offbore = LNNSS_6SV.TX_offbore(r,k);
                TX_azimuth = LNNSS_6SV.TX_azimuth(r,k);

                RX2TX_offbore = LNNSS_6SV.RX_offbore(r,k);

                if  LNNSS_6SV.visibility_flag(r,k) == 1 && TX_offbore <= 90

                    % Computation of the free space path losses for C/A L1 signal

                    d = LNNSS_6SV.RX_TXdistance(r,k); % distance of the RX wrt the TX

                    A_d = FSP_Loss(d,f,c); % free space path losses [dB]


                    % Evaluation of the RX gain, function of the off-boresigth angle

                     G_RX = RX_GAIN(RX2TX_offbore);

                    % Evaluation of the TX gain, function of the off-boresigth
                    % angle and the SV attitude


                    G_TX = interp2(X_grid,Y_grid,G,TX_azimuth, TX_offbore/(1.1));


                    % For an isotropic antenna the T_a is a weighted mean of the temperature
                    % for a unitary sphere integration


                    et = LNNSS_6SV.epoch(r);

                    pos = cspice_spkpos ( 'EARTH', et,'J2000', 'LT', ID_RX);
                    d_E = norm(pos);

                    pos = cspice_spkpos ( 'MOON', et,'J2000', 'LT', ID_RX);
                    d_M = norm(pos);

                    % Earth solid angle seen from the RX 

                    theta = atan(R_E/d_E);
                    OM_E = 2*pi * ( 1 - cos(theta));

                    % Moon solid angle seen from the RX 

                    theta = atan(R_M/d_M);
                    OM_M = 2*pi * ( 1 - cos(theta));

                    % Antenna temperature 

                    T_ant = (1/(4*pi)) * (OM_E * (25 * T_E) + OM_M * T_M + (4*pi - (OM_E + OM_M)) * T_DS);

                    % The antennta noise temperature including thermal inefficiencies

                    T_ant =  T_ant   + 290 * (1/ANT_e -1);

                    % Low Noise Amplifier (LNA) noise temperature

                    T_amp = 290*(10^(N_f/10) - 1); % For an assumed 290 K  temperature 

                    % Total system temperature

                    T_sys = T_ant + T_amp;

                    % Carrier signal [dBW]

                    C = P_TX - A_d + G_TX + G_RX - L_B;


                    % Noise [dBW - Hz]

                    N_0 = 10*log10( k_b * ( T_sys));

                    C_N0 = C - N_0; % dB - Hz

                    LNNSS_6SV.C_N0(r,k) = C_N0;

                    LNNSS_6SV.acquired(r,k) = C_N0 > aq_treshold;


                    % Doppler shift computation

                    LNNSS_6SV.D_shift(r,k) = doppler_shift(f,LNNSS_6SV.RX_TXv_radial(r,k),c);

                end


            end

        end


        %% GDOP

        % We can compute under a certain set of hypothesis the GDOP as in "GPS
        % principles and applications" (Kaplan)

        % At each time step 

        for r = 1 : numel(LNNSS_6SV.epoch)

        et = LNNSS_6SV.epoch(r);

        % Number of acquired and tracked SV at the time considered

        SV_tracked = sum(LNNSS_6SV.acquired(r,:));

        LNNSS_6SV.n_tracked(r,1) = SV_tracked;

        H = zeros(SV_tracked,4);

        n = 1;

        % For each SV of the constellation

            for k = 1 : 6

                % if the SV is being tracked

                if LNNSS_6SV.acquired(r,k) > 0 && SV_tracked > 3

                    LNNSSID = -(800 + k);
                    LNNSSID = char(string(LNNSSID));  

                    % Position vector of the tracked SV as seen from the TX

                    a_n = cspice_spkpos ( LNNSSID, et,'J2000', 'XLT + S', ID_RX);

                    % Normalizzation of the position vector

                    a_n = a_n' / norm(a_n);

                    H(n,:) = [a_n 1];

                    n = n + 1;

                end

            end


           if SV_tracked > 3

                D = inv(H' * H);

                GDOP = sqrt(trace(D));

                LNNSS_6SV.GDOP(r,1) = GDOP;

           else

               LNNSS_6SV.GDOP(r,1) =  nan;

           end
        end



        % 8 SVs constellation

        if isequal(LNNSS,'H')

        cspice_unload('kernels/spk/6SV_H_LNNSS.bsp');
        cspice_furnsh('kernels/spk/8SV_H_LNNSS.bsp');

        elseif isequal(LNNSS,'L')

        cspice_unload('kernels/spk/6SV_L_LNNSS.bsp');
        cspice_furnsh('kernels/spk/8SV_L_LNNSS.bsp');

        end


        %% LNNSS constellation

        for k = 1:8

            LNNSSID = -(800 + k);
            LNNSSID = char(string(LNNSSID));  


            % Time windows of Earth occultation

            occtyp  = 'any';
            front   = 'earth';
            fshape  = 'ellipsoid';
            fframe  = 'ITRF93';
            back    =  LNNSSID;
            bshape  = 'point';
            bframe  = '-';
            obsrvr  = ID_RX;
            abcorr  = 'LT+S';


            EARTH_occultation = cspice_gfoclt( occtyp, front, fshape, fframe, ...
                                  back, bshape,bframe,          ...
                                  abcorr, obsrvr, step_OCC, cnfine,  ...
                                  MAXWIN);



            % List the beginning and ending times in each interval
            % if result contains data.


            left_E = zeros(1,numel( EARTH_occultation)/2 + 1);
            right_E = zeros(1,numel( EARTH_occultation)/2 + 1);

            for i = 1:numel(EARTH_occultation)/2

             [left_E(i), right_E(i)] = cspice_wnfetd( EARTH_occultation, i );

            end


            % Adding to the windows the endtime of the simulation, needed to
            % compute after the last occultation window or if the window is empty

            left_E(end) = tspan(end);
            right_E(end) = tspan(end);


            % Time windows of Moon occultation

            occtyp  = 'any';
            front   = 'moon';
            fshape  = 'ellipsoid';
            fframe  = 'moon_PA';
            back    =  LNNSSID;
            bshape  = 'point';
            bframe  = '-';
            obsrvr  = ID_RX;
            abcorr  = 'LT+S';


            MOON_occultation = cspice_gfoclt( occtyp, front, fshape, fframe, ...
                                  back, bshape,bframe,          ...
                                  abcorr, obsrvr, step_OCC, cnfine,  ...
                                  MAXWIN);



            % List the beginning and ending times in each interval
            % if result contains data.



            left_M = zeros(1,numel( MOON_occultation)/2 + 1);
            right_M = zeros(1,numel( MOON_occultation)/2 + 1);

            for i = 1:numel(MOON_occultation)/2

             [left_M(i), right_M(i)] = cspice_wnfetd( MOON_occultation, i );

            end

            % Adding to the windows the endtime of the simulation, needed to
            % compute after the last occultation window or if the window is empty

            left_M(end) = tspan(end);
            right_M(end) = tspan(end);



            % For each time step where the TX is in view of the SV, the ideal
            % range, relative velocity and angular separation with respect to the
            % SV nadir in evaluated

            % initializzation of the Earth and Moon occultation windows counter

             n = 1; % Earth occultation windows counter
             m = 1; % Moon occultation windows counter


            for r = 1 : numel(tspan) % for each time step

                 et = tspan(r);

                 % saving the epoch in the output structure

                  LNNSS_8SV.epoch(r,1) = et ;


                 % If the timestep is not inside an occultation window the
                 % computation are performed

                 if et <= left_E(n) && et <= left_M(m)

                     % Computing all the relevant position and state vector at the
                     % current timestep 

                     % State of the RX as seen from the TX SV

                     [state,ltime] = cspice_spkezr (ID_RX, et,'J2000', 'XLT',  LNNSSID);

                     TX2RX_pos_J2000 = state(1:3,1); % Position vector
                     TX2RX_vel_J2000 = state(4:6,1); % Velocity vector

                     % Position of MOON as seen from the TX SV ( an its array pointing)

                     TX2MO_pos_J2000 = cspice_spkpos ( 'MOON', et  ,'J2000', 'LT',  LNNSSID);

                     % Position of the SUN as seen from the TX SV 

                     TX2SUN_pos_J2000 = cspice_spkpos ( 'SUN', et  ,'J2000', 'LT',  LNNSSID); 

                     % Only the Sun direction is needed 

                     TX2SUN_pos_J2000 = cspice_vhat(TX2SUN_pos_J2000); % vector normalizzation 

                     % Position of Earth as seen from the RX ( an its high gain
                     % antenna pointing)

                     RX2EA_pos_J2000 = cspice_spkpos ( 'EARTH', et ,'J2000', 'LT', ID_RX);

                     % We can compute the distance between the TX and RX as the norm of the position vector    

                     d = norm(TX2RX_pos_J2000);

                     range = ltime * c; % equivalently is the light time multiplied by the speed of light

                     % Radial velocity from the SV to the RX

                     v_radial = dot(TX2RX_vel_J2000,(TX2RX_pos_J2000/d));

                     % To compute the off-bore angle and the azimuth of RX wrt 3D
                     % antenna pattern of each SV a simple yaw-steering (YS) attitude model is
                     % implemented 

                     % For the LNNSS constellation is in similarity to the
                     % GALILEO constellation only the SVs are poiting the array
                     % toward the Moon center.


                     % In nominal operation the z-axis is pointing the Moon center,
                     % the x-axis is toward the Sun direction and the y-axis complete the orthonormal basis


                     k_YS = cspice_vhat(TX2MO_pos_J2000);

                     j_YS = cross(TX2SUN_pos_J2000,TX2MO_pos_J2000)/norm(cross(TX2SUN_pos_J2000,TX2MO_pos_J2000));

                     i_YS = cross(j_YS,k_YS);

                     % A rotation matrix can be compute to rotate a vector in
                     % the SV YS fram ( ideally equivalent to the body frame)

                     R_J20002YS = [i_YS'; j_YS'; k_YS'];

                     % The direction vector of the RX in the TX body frame can
                     % be computed rotating the RX2RX position vector from the
                     % J2000 frame 

                     TX2RX_pos_YS = R_J20002YS * cspice_vhat(TX2RX_pos_J2000); 

                     % The angle can be computed as a transformation into
                     % spherical coordinates where:

                     % Theta is the angle between the z-axis ( nadir pointing)
                     % and the RX direction in the YS_frame. This is the
                     % off-bore angle. It is bounded in [0, pi];

                     % Phi is the angle between the +x_YS axis and the
                     % projection of the RX direction into the xy plane. It is
                     % consistend with the convention for the azimuth of the
                     % directivity pattern and is the azimuth. It is bounded in [-pi, pi];

                     [~,theta,phi] = cspice_recsph(TX2RX_pos_YS);

                     TX2RX_offbore = cspice_convrt ( theta,'RADIANS', 'DEGREES');
                     TX2RX_azimuth = cspice_convrt ( phi,'RADIANS', 'DEGREES');



                     % The Gain pattern data has the azimuth defined in the range
                     % [0 2pi]

                     if  TX2RX_azimuth < 0

                          TX2RX_azimuth = 2*pi*dpr + TX2RX_azimuth;

                     end

                     % Angular separation between the Earth pointing antenna of the
                     % receiver SPC and the position of the SV ( off-bore angle)

                     RX2TX_offbore = cspice_convrt ( cspice_vsep(RX2EA_pos_J2000, -TX2RX_pos_J2000),'RADIANS', 'DEGREES'); 


                     LNNSS_8SV.RX_TXdistance(r,k) = d ;
                     LNNSS_8SV.RX_TXv_radial(r,k) = v_radial ;
                     LNNSS_8SV.TX_offbore(r,k) = TX2RX_offbore ;
                     LNNSS_8SV.TX_azimuth(r,k) = TX2RX_azimuth ;
                     LNNSS_8SV.RX_offbore(r,k) = RX2TX_offbore ;
                     LNNSS_8SV.visibility_flag(r,k) = 1 ;


                 else

                     % else if the RX is in occultation or by Earth or Moon the SV
                     % is not visible

                     LNNSS_8SV.visibility_flag(r,k) = 0 ;

                 end

                 % Once the time step surpass a time window the counter is increased to analyze the next one

                 if et >= right_E(n) && n < length(right_E)
                     n = n + 1;
                 end

                 if et >= right_M(m) && m < length(right_M)
                     m = m + 1;
                 end


            end

            % In the first step the occultation flag is always 1 even if there is
            % no visibility due to the way the logic check for the spice
            % windows is perfomed, so an extra check for the first row is needed for
            % each SV

            if (left_E(1) - LNNSS_8SV.epoch(1,1)) == 0

                LNNSS_8SV.visibility_flag(1,k) = 0 ;

            end



        end


        %% LNNSS link budget


        % E5a signal (Carrier frequency equal to L5)

        f = 1176.45 * 1e6; % [Hz]

        % Power of the signal in the considered bandwidth

        p_sharing = 0.5;

        % For the GALILEO constellation partial and not official data both for the
        % gain and the transmitted power are present 

        G = BIIF_GAIN; % Using the GPS IIF/III-A gain pattern as reference for Galileo constellation 

        X_grid = X_grid_IIF;
        Y_grid = Y_grid_IIF;

        % For each SV 

        for k = 1:8

            % For each time step

            % Mean emitted power of the TX for GALILEO GSAT SV (source "GNSS satellite
            % transmit power and its impact on orbit % determination")


            P_TX = 10*log10(p_sharing * 265/4); % 66.25W mean transmitted power in dB considering the power in the bandwidth 


            for r = 1 : numel(tspan)

                 % If the SV is visible and with a separation less than 90 deg

                TX_offbore = LNNSS_8SV.TX_offbore(r,k);
                TX_azimuth = LNNSS_8SV.TX_azimuth(r,k);

                RX2TX_offbore = LNNSS_8SV.RX_offbore(r,k);

                if  LNNSS_8SV.visibility_flag(r,k) == 1 && TX_offbore <= 90

                    % Computation of the free space path losses for C/A L1 signal

                    d = LNNSS_8SV.RX_TXdistance(r,k); % distance of the RX wrt the TX

                    A_d = FSP_Loss(d,f,c); % free space path losses [dB]


                    % Evaluation of the RX gain, function of the off-boresigth angle

                     G_RX = RX_GAIN(RX2TX_offbore);

                    % Evaluation of the TX gain, function of the off-boresigth
                    % angle and the SV attitude


                    G_TX = interp2(X_grid,Y_grid,G,TX_azimuth, TX_offbore/(1.1));


                    % For an isotropic antenna the T_a is a weighted mean of the temperature
                    % for a unitary sphere integration


                    et = LNNSS_8SV.epoch(r);

                    pos = cspice_spkpos ( 'EARTH', et,'J2000', 'LT', ID_RX);
                    d_E = norm(pos);

                    pos = cspice_spkpos ( 'MOON', et,'J2000', 'LT', ID_RX);
                    d_M = norm(pos);

                    % Earth solid angle seen from the RX 

                    theta = atan(R_E/d_E);
                    OM_E = 2*pi * ( 1 - cos(theta));

                    % Moon solid angle seen from the RX 

                    theta = atan(R_M/d_M);
                    OM_M = 2*pi * ( 1 - cos(theta));

                    % Antenna temperature 

                    T_ant = (1/(4*pi)) * (OM_E * (25 * T_E) + OM_M * T_M + (4*pi - (OM_E + OM_M)) * T_DS);

                    % The antennta noise temperature including thermal inefficiencies

                    T_ant =  T_ant   + 290 * (1/ANT_e -1);

                    % Low Noise Amplifier (LNA) noise temperature

                    T_amp = 290*(10^(N_f/10) - 1); % For an assumed 290 K  temperature 

                    % Total system temperature

                    T_sys = T_ant + T_amp;

                    % Carrier signal [dBW]

                    C = P_TX - A_d + G_TX + G_RX - L_B;


                    % Noise [dBW - Hz]

                    N_0 = 10*log10( k_b * ( T_sys));

                    C_N0 = C - N_0; % dB - Hz

                    LNNSS_8SV.C_N0(r,k) = C_N0;

                    LNNSS_8SV.acquired(r,k) = C_N0 > aq_treshold;


                    % Doppler shift computation

                    LNNSS_8SV.D_shift(r,k) = doppler_shift(f,LNNSS_8SV.RX_TXv_radial(r,k),c);

                end


            end

        end


        %% GDOP

        % We can compute under a certain set of hypothesis the GDOP as in "GPS
        % principles and applications" (Kaplan)

        % At each time step 

        for r = 1 : numel(LNNSS_8SV.epoch)

        et = LNNSS_8SV.epoch(r);

        % Number of acquired and tracked SV at the time considered

        SV_tracked = sum(LNNSS_8SV.acquired(r,:));

        LNNSS_8SV.n_tracked(r,1) = SV_tracked;

        H = zeros(SV_tracked,4);

        n = 1;

        % For each SV of the constellation

            for k = 1 : 8

                % if the SV is being tracked

                if LNNSS_8SV.acquired(r,k) > 0 && SV_tracked > 3

                    LNNSSID = -(800 + k);
                    LNNSSID = char(string(LNNSSID));  

                    % Position vector of the tracked SV as seen from the TX

                    a_n = cspice_spkpos ( LNNSSID, et,'J2000', 'XLT + S', ID_RX);

                    % Normalizzation of the position vector

                    a_n = a_n' / norm(a_n);

                    H(n,:) = [a_n 1];

                    n = n + 1;

                end

            end


           if SV_tracked > 3

                D = inv(H' * H);

                GDOP = sqrt(trace(D));

                LNNSS_8SV.GDOP(r,1) = GDOP;

           else

               LNNSS_8SV.GDOP(r,1) =  nan;

           end
        end



    % GPS + GALILEO + LNNSS

        %% 6 SV LNNSS

        % At each time step 

        if isequal(LNNSS,'H')

        cspice_unload('kernels/spk/8SV_H_LNNSS.bsp');
        cspice_furnsh('kernels/spk/6SV_H_LNNSS.bsp');

        elseif isequal(LNNSS,'L')

        cspice_unload('kernels/spk/8SV_L_LNNSS.bsp');
        cspice_furnsh('kernels/spk/6SV_L_LNNSS.bsp');

        end


        for r = 1 : numel(GPS.epoch)

        et = GPS.epoch(r);

        % Number of acquired and tracked SV at the timeconsidered

        SV_tracked_GPS = sum(GPS.acquired(r,:));
        SV_tracked_GALILEO = sum(GALILEO.acquired(r,:));
        SV_tracked_LNNSS = sum(LNNSS_6SV.acquired(r,:));
        % SV_tracked_MOON_GALILEO = sum(MOON_GALILEO_8SV.acquired(r,:));

        SV_tracked = (SV_tracked_GPS + SV_tracked_GALILEO + SV_tracked_LNNSS);

        H = zeros(SV_tracked,4);

        n = 1;

        % For each SV of the GPS constellation

            for k = 1 : 30

                % if the SV is being tracked

                if GPS.acquired(r,k) > 0 && SV_tracked > 3

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

                if GALILEO.acquired(r,k) > 0 && SV_tracked > 3

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


            for k = 1 : 6

                % if the SV is being tracked

                if LNNSS_6SV.acquired(r,k) > 0 && SV_tracked > 3

                    LNNSSID = -(800 + k);
                    LNNSSID = char(string(LNNSSID));  

                    % Position vector of the tracked SV as seen from the TX

                    a_n = cspice_spkpos ( LNNSSID, et,'J2000', 'XLT + S', ID_RX);

                    % Normalizzation of the position vector

                    a_n = a_n' / norm(a_n);

                    H(n,:) = [a_n 1];

                    n = n + 1;

                end

            end

            if SV_tracked > 3

                D = inv(H' * H);

                GDOP = sqrt(trace(D));

                MULTI_GNSS.GPS_GAL_LNNSS6SV_GDOP(r,1) = GDOP;

             else

                MULTI_GNSS.GPS_GAL_LNNSS6SV_GDOP(r,1) =  nan;

            end
        end


        %% 8 SVs LNNSS

        if isequal(LNNSS,'H')

        cspice_unload('kernels/spk/6SV_H_LNNSS.bsp');
        cspice_furnsh('kernels/spk/8SV_H_LNNSS.bsp');

        elseif isequal(LNNSS,'L')

        cspice_unload('kernels/spk/6SV_L_LNNSS.bsp');
        cspice_furnsh('kernels/spk/8SV_L_LNNSS.bsp');

        end


        % At each time step

        for r = 1 : numel(GPS.epoch)

        et = GPS.epoch(r);

        % Number of acquired and tracked SV at the timeconsidered

        SV_tracked_GPS = sum(GPS.acquired(r,:));
        SV_tracked_GALILEO = sum(GALILEO.acquired(r,:));
        SV_tracked_LNNSS = sum(LNNSS_8SV.acquired(r,:));

        SV_tracked = (SV_tracked_GPS + SV_tracked_GALILEO + SV_tracked_LNNSS);

        H = zeros(SV_tracked,4);

        n = 1;

        % For each SV of the GPS constellation

            for k = 1 : 30

                % if the SV is being tracked

                if GPS.acquired(r,k) > 0 && SV_tracked > 3

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

                if GALILEO.acquired(r,k) > 0 && SV_tracked > 3

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


            for k = 1 : 8

                % if the SV is being tracked

                if LNNSS_8SV.acquired(r,k) > 0 && SV_tracked > 3

                    LNNSSID = -(800 + k);
                    LNNSSID = char(string(LNNSSID));  

                    % Position vector of the tracked SV as seen from the TX

                    a_n = cspice_spkpos ( LNNSSID, et,'J2000', 'XLT + S', ID_RX);

                    % Normalizzation of the position vector

                    a_n = a_n' / norm(a_n);

                    H(n,:) = [a_n 1];

                    n = n + 1;

                end

            end

            if SV_tracked > 3

                D = inv(H' * H);

                GDOP = sqrt(trace(D));

                MULTI_GNSS.GPS_GAL_LNNSS8SV_GDOP(r,1) = GDOP;

             else

                MULTI_GNSS.GPS_GAL_LNNSS8SV_GDOP(r,1) =  nan;

            end
        end

    end

%% clearing the kernel pool

cspice_kclear;

end

