function [BIIR_GAIN,BIIRM_GAIN,BIIF_GAIN] = GPS_GAIN_data(plot_flag)

% Function to extrapolated the data of the gain pattern of the GPS
% constellation (for each block) in an usable form starting from the
% published data avaible.


% PROTOTYPE:
    % [BIIR_GAIN,BIIRM_GAIN,BIIF_GAIN] = GPS_GAIN_data(1);
    
% INPUT:
   % plot_flag: plotting flag, 1:plot, 0: no plot
% OUTPUT:
    % BIIR_GAIN: matrix containting antenna gain values for each off-bore
    % and zimuth slice available. GPS Block IIR SVs.
    % BIIRM_GAIN: matrix containting antenna gain values for each off-bore
    % and zimuth slice available. GPS Block IIR-M SVs.
    % BIIF_GAIN: matrix containting antenna gain values for each off-bore
    % and zimuth slice available. GPS Block IIF SVs.


% VERSIONS:
    % 1/4/2021: First version

    
%% B-IIF Block

% mean gain pattern data for BII-F SV avaible 

BIIF = readmatrix('TX_GAIN\Block-IIF_txgain_gpsace.txt');

% Off-bore angle

theta = (0:1:90);

BIIF_GAIN = zeros(80,360);

% The data avaible is down to 16 deg of off-bore angle. To complete the
% pattern mean azimuthal values are included from other sources

BIIF_GAIN(1,:)  = 13;         BIIF_GAIN(11,:) = 15;
BIIF_GAIN(2,:)  = 13.54;      BIIF_GAIN(12,:) = 14.89;
BIIF_GAIN(3,:)  = 13.98;      BIIF_GAIN(13,:) = 14.59;
BIIF_GAIN(4,:)  = 14.31;      BIIF_GAIN(14,:) = 14.15;
BIIF_GAIN(5,:)  = 14.57;      BIIF_GAIN(15,:) = 13.61;
BIIF_GAIN(6,:)  = 14.75;      BIIF_GAIN(16,:) = 13;
BIIF_GAIN(7,:)  = 14.87;
BIIF_GAIN(8,:)  = 14.95;
BIIF_GAIN(9,:)  = 14.98;
BIIF_GAIN(10,:) = 15;

% Extracting data from the table

for k = 17:91
    
     point = BIIF(:,1) == theta(k);   

     BIIF_GAIN(k,:) = BIIF(point,3);
                            
end

% Adding a column for 2pi azimuth ( equal to 0 azimuth) for interpolation
% purposes

BIIF_GAIN(:,361) = BIIF_GAIN(:,1);

% Plotting some azimuthal slices

if plot_flag == 1
    
    for z = 10:10:359
              
        figure(1)
        subplot(1,2,1)
        hold on 
        grid on

        plot(theta,BIIF_GAIN(:,1),'linewidth',1.5)   

        plot(theta,BIIF_GAIN(:,z),'linewidth',1.5) 
        
    end
    set(gca, 'FontSize', 20) ;  
    title('Block IIF L1 antenna gain pattern') ; 
    xlabel('off-bore angle [Deg]')
    ylabel('Gain [dB]') 
     set(gca, 'FontSize', 16) ; 
     legend('?? = 10 ','?? = 20 ','?? = 30 ','?? = 40 ','?? = 50 ','?? = 60 ','?? = 70 ','?? = 80 ','?? = 90 ','?? = 100 ','?? = 110 ',...
           '?? = 120','?? = 130 ','?? = 140 ','?? = 150 ','?? = 160 ','?? = 170 ','?? = 180 ','?? = 190 ','?? = 200 ','?? = 110 ',...
           '?? = 220 ','?? = 230 ','?? = 240 ','?? = 250 ','?? = 260 ','?? = 270 ','?? = 280 ','?? = 290 ','?? = 300 ','?? = 310 ',...
           '?? = 320 ','?? = 330 ','?? = 340 ','?? = 350 ');

    V = axis;
    V(1) = 0;
    V(2) = 90;
    V(4) = 16;
    axis(V);
    
    subplot(1,2,2)
    hold on 
    grid on
    
    [X_F,Y_F] = meshgrid(0:1:360,0:1:90);
    surf(X_F,Y_F,BIIF_GAIN);
    
    set(gca, 'FontSize', 18) ;  
    title('BIIF L1 antenna gain pattern') ;
    xlabel('azimuth angle [Deg]', 'Interpreter', 'Latex') 
    ylabel('off-bore angle [Deg]', 'Interpreter', 'Latex')
    zlabel('Gain [dB]', 'Interpreter', 'Latex') 

end


%% B-IIR

% data for four BIIR SV avaible (modernized panels)

SVN47 = readtable('TX_GAIN\L1_BIIR_47.xlsx','range','A4:AK94');
SVN59 = readtable('TX_GAIN\L1_BIIR_59.xlsx','range','A4:AK94');
SVN60 = readtable('TX_GAIN\L1_BIIR_60.xlsx','range','A4:AK94');
SVN61 = readtable('TX_GAIN\L1_BIIR_61.xlsx','range','A4:AK94');


% Off-bore angle

theta = SVN47.(1);

BIIR_GAIN = zeros(91,36);

% Gain correction factor

GCF = 1.3; % [dB]

% Mean value of the directivity of the block in each azimuth slice avaible

for k = 2:37
    
     BIIR_GAIN(:,k-1) = (SVN47.(k) + SVN59.(k) + SVN60.(k) + SVN61.(k))/4 - GCF;  
                         
end

% Adding a column for 2pi azimuth ( equal to 0 azimuth) for interpolation
% purposes

BIIR_GAIN(:,37) = BIIR_GAIN(:,1);

 if plot_flag == 1
     
    figure(2)

    subplot(1,2,1)
    hold on 
    grid on
    
    plot(theta,BIIR_GAIN(:,1 : (end-1)),'linewidth',1.5)    
    
    set(gca, 'FontSize', 20) ;  
    title('Block IIR L1 antenna gain pattern') ;
    xlabel('off-bore angle [Deg]')
    ylabel('Gain [dB]') 

    V = axis;
    V(1) = -90;
    V(2) =  90;
    V(4) =  16;
    axis(V);
    
    
    subplot(1,2,2)
    hold on 
    grid on
    
    [X_R,Y_R] = meshgrid(0:10:360,-90:2:90);
    surf(X_R,Y_R,BIIR_GAIN);
    
    set(gca, 'FontSize', 18) ;  
    title('Block IIR L1 antenna gain pattern') ;
    xlabel(' azimuth angle [Deg]', 'Interpreter', 'Latex') 
    ylabel('off-bore angle [Deg]', 'Interpreter', 'Latex')
    zlabel('Gain [dB]', 'Interpreter', 'Latex') 
    

 end

 
%% B-IIR-M

% data for nine B-IIR-M SV avaible

SVN48 = readtable('TX_GAIN\L1_IIRM_48.xlsx','range','A4:AK94');
SVN49 = readtable('TX_GAIN\L1_IIRM_49.xlsx','range','A4:AK94');
SVN50 = readtable('TX_GAIN\L1_IIRM_50.xlsx','range','A4:AK94');
SVN51 = readtable('TX_GAIN\L1_IIRM_51.xlsx','range','A4:AK94');
SVN52 = readtable('TX_GAIN\L1_IIRM_52.xlsx','range','A4:AK94');
SVN53 = readtable('TX_GAIN\L1_IIRM_53.xlsx','range','A4:AK94');
SVN55 = readtable('TX_GAIN\L1_IIRM_55.xlsx','range','A4:AK94');
SVN57 = readtable('TX_GAIN\L1_IIRM_57.xlsx','range','A4:AK94');
SVN58 = readtable('TX_GAIN\L1_IIRM_58.xlsx','range','A4:AK94');


BIIRM_GAIN = zeros(91,36);

% Gain correction factor

GCF = 1.4; % [dB]

% Mean value of the directivity of the block in each azimuth slice avaible


for k = 2:37
    
     BIIRM_GAIN(:,k-1) = (SVN48.(k) + SVN49.(k) + SVN50.(k) + SVN51.(k) + SVN52.(k) + ...
                                 SVN53.(k)  + SVN55.(k) + SVN57.(k) + SVN58.(k))/9 - GCF;  

end

% Adding a column for 2pi azimuth ( equal to 0 azimuth) for interpolation
% purposes

BIIRM_GAIN(:,37) = BIIRM_GAIN(:,1);


if plot_flag == 1

    figure(3)
    
    subplot(1,2,1)
    hold on 
    grid on

    plot(theta,BIIRM_GAIN(:,1 : (end-1)),'linewidth',1.5)    

    set(gca, 'FontSize', 20) ;  
    title('Block IIR-M L1 antenna gain pattern') ;
    xlabel('off-bore angle [Deg]')
    ylabel('Gain [dB]') 

    V = axis;
    V(1) = -90;
    V(2) =  90;
    V(4) =  16;
    axis(V);
    
    
    subplot(1,2,2)
    hold on 
    grid on
    
    surf(X_R,Y_R,BIIRM_GAIN);
    
    set(gca, 'FontSize', 18) ;  
    title('BIIR-M L1 antenna gain pattern') ;
    xlabel('azimuth angle [Deg]', 'Interpreter', 'Latex') 
    ylabel('off-bore angle [Deg]', 'Interpreter', 'Latex')
    zlabel('Gain [dB]', 'Interpreter', 'Latex') 
    
end

end

