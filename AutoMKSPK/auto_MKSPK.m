% Script to automate the mkpsk call in the cmd prompt in order to automate
% the generation of the spk of the GNSS Constallation

% The excell contains the path direction of the setup files and should be
% loaded manually as a string

%% TLE from celestrak
% GPS

for k = 1:30
    
    system(GNSSspk(k,1));
    
end

% Merging the generated spk in a single one
system(GNSSspk(k + 1,1));

system('brief -t c:\users\gobli\desktop\TLE_GPS.bsp &');


% GALILEO

for k = 1:24
    
    system(GNSSspk(k,2));
    
end

% Merging the generated spk in a single one
system(GNSSspk(k + 1,2));

system('brief -t c:\users\gobli\desktop\TLE_GALILEO.bsp &');


% GLONASS

for k = 1:24
    
    system(GNSSspk(k,3));
    
end

% Merging the generated spk in a single one
system(GNSSspk(k + 1,3));

system('brief -t c:\users\gobli\desktop\TLE_GLONASS.bsp &');


% BEIDOU III

for k = 1:30
    
    system(GNSSspk(k,4));
    
end

% Merging the generated spk in a single one
system(GNSSspk(k + 1,4));

system('brief -t c:\users\gobli\desktop\TLE_BEIDOU.bsp &');


%% Propagated states

% SPK Generated using the TLE from celestrak
% GPS

for k = 1:30
    
    system(STATESGNSSspk(k,1));
    
end

% Merging the generated spk in a single one
system(STATESGNSSspk(k + 1,1));

system('brief -t c:\users\gobli\desktop\2K25_GPS.bsp &');


% SPK Generated using the TLE from celestrak
% GALILEO

for k = 1:24
    
    system(STATESGNSSspk(k,2));
    
end

% Merging the generated spk in a single one
system(STATESGNSSspk(k + 1,2));

system('brief -t c:\users\gobli\desktop\2K25_GALILEO.bsp &');

% SPK Generated using the TLE from celestrak
% GLONASS

for k = 1:24
    
    system(STATESGNSSspk(k,3));
    
end

% Merging the generated spk in a single one
system(STATESGNSSspk(k + 1,3));

system('brief -t c:\users\gobli\desktop\2K25_GLONASS.bsp &');


% SPK Generated using the TLE from celestrak
% BEIDOU III

for k = 1:30
    
    system(STATESGNSSspk(k,4));
    
end

% Merging the generated spk in a single one
system(STATESGNSSspk(k + 1,4));

system('brief -t c:\users\gobli\desktop\2K25_BEIDOU.bsp &');




