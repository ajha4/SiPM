function [meas_sipm_par,fixed_sipm_par,crystal_par,flags] = read_input_array()
%Read all the input parameters of all the SiPMs connected to the
%scintillator. Typically, to a single scintillator, multiple SiPMs are
%connected, so in this code, we read the measured parameters for all the SiPMs.
%A lot of symbols are the same as used and defined in the following paper:
%Jha, A. K., van Dam, H. T., Kupinski, M. A., & Clarkson, E. (2013).
%Simulating Silicon Photomultiplier Response to Scintillation Light.
%IEEE Transactions on Nuclear Science, 30(1), 336?351. http://doi.org/10.1109/TNS.2012.2234135

sipm_no = 16; % No of SiPMs

Td = 16e-9; % Decay time of the scintillator.
Tr = 0.93e-9; % Rise time of the scintillator
Qsc = 70; % Photon yield for LaBR3
f = 1; % Collection efficiency of the SiPM
Ncells = 3600; % The number of cells in the SiPM
tmeas = 500e-9; %The entire measurement time
%We now define the geometry of the system. The SiPMs are placed in the x-y
%plane.
xmax = 1.5e-3; %The maximum coordinate along the x axis
xmin = -1.5e-3; %The mimimum coordinate along the x axis
ymax = 1.5e-3; %Same for y axis.
ymin = -1.5e-3;

sipm_pix_xno = sqrt(Ncells);
sipm_pix_yno = sqrt(Ncells);
wid_mcell_x = (xmax - xmin)/sipm_pix_xno;
wid_mcell_y = (ymax - ymin)/sipm_pix_yno;

% The parameters that will typically be the same for all SiPMs.
fixed_sipm_par = struct('Ncells',3600,'tau_aval',0.1e-9,'tau_phd2',0.1e-12,'tau_phr',0.1e-9,'Rs',15,...
    'wid_mcell_x',wid_mcell_x,'wid_mcell_y',wid_mcell_y,'sipm_pix_xno',60,...
    'sipm_pix_yno',60,'tmeas',tmeas);

% The parameters that describe the scintillation pulse
crystal_par = struct('Td',0,'Qsc',0);
crystal_par.Td = Td; %Decay time of scintillation pulse
crystal_par.Tr = Tr; %Rise time of scintillation pulse
crystal_par.Qsc = Qsc;

%Simulation flags
flags = struct('oct',1,'ap',1,'deadtime',1,'overvoltage',1,'output_voltage',1,'dc',1,'partial_recovery',1,...
    'overvoltage_pde',0,'overvoltage_pde_meas',0,'overvoltage_pde_lin',0,...
    'overvoltage_dc',0,'overvoltage_dc_meas',0,'overvoltage_dc_lin',0,...
    'overvoltage_oct',0,'overvoltage_oct_meas',0,'overvoltage_oct_quad',0,...
    'overvoltage_ap',0,'overvoltage_ap_meas',0,'overvoltage_ap_quad',0);

% The measured parameters of the SiPM at V_{ob} = V_{ob0}. Many a times,
% the values at V_{ob0} are the only ones that are available. Therefore,
% creating a separate structure for them.

% All the error bars correspond to 95 % confidence intervals, which corresponds to 2*sigma
qe = [0.157 0.172]; % The internal microcell QE
Nbar_oct_mean =  [0.140 0.132]; % The mean number of optical cross-talk events when one discharge occurs
Nbar_oct_var = [0.005 0.005]/2; % % The standard deviation of optical cross-talk events when one discharge occurs
Nbar_ap_short_mean = [0.124 0.132]; % The mean number of afterpulse events when one discharge occurs
Nbar_ap_short_var =  [0.005 0.005]/2; %The corresponding standard deviation
Nbar_ap_long_mean = [0 0 ]; % This was not computed in the experimental setup in Herman et al.
Nbar_ap_long_var =  [0 0];

G_mean = [8.3*1e5 7.9*1e5]; % Average gain of the SiPM
G_var = [0.5*1e5 0.5*1e5]/2; %Standard deviation of the gain
rdc_mean = [5.6*1e6 5.1*1e6]; % Rate of dark noise
rdc_var = [0.1e6 0.1e6]/2; %Standard deviation of dark noise
tau_phd1_mean = [15.3e-9 12.8e-9]; % The mean recovery time of the SiPM
tau_phd1_var = [1e-9 1e-9]/2; %Standard dev. of the recovery time
tau_ap_short_mean = [25e-9 26e-9]; % Mean decay time of the afterpulse event
tau_ap_short_var = [8e-9 8e-9]/2; % Standard deviation of the decay time of afterpulse event
tau_ap_long_mean = tau_ap_short_mean;
tau_ap_long_var = [0 0];
Vob_mean = [ 1.26 1.39]; % The measured voltage over breakdown V_{ob0}
Vob_var = [0.05 0.05]/2; %Symbols as defined in the paper
Rq_mean = [145.3e3 140.4e3];
Rq_var = [ 0.5e3 0.5e3]/2;

meas_sipm_par(sipm_no) = struct('qe',0,'Nbar_oct_mean',0, 'Nbar_oct_var',0,...
    'Nbar_ap_short_mean',0,'Nbar_ap_short_var',0,'Nbar_ap_long_mean',0,'Nbar_ap_long_var',0,...
    'G_mean',0,'G_var',0,'rdc_mean',0,'rdc_var',0,...
    'tau_phd1_mean',0, 'tau_phd1_var',0,'tau_ap_short_mean',0, 'tau_ap_short_var',0,...
    'tau_ap_long_mean',0, 'tau_ap_long_var',0,...
    'V_ob_mean',0,'V_ob_var',0,'f',0,'Rq_mean',0,'Rq_var',0,...
    'pde_ratio',0,'ov_ratio_pde',0);
%     'xmax',1.5e-3,'xmin',-1.5e-3,'ymax',1.5e-3,'ymin',-1.5e-3);

flag_default_variation_parameters_overvoltage = 1;% This flag can be set to zero if
% the user wants to input specific variations with respect to overvoltage. 
% For example, if the user has measurements for how the different
% SiPM parameters, such as the afterpulsing rate, dark current etc. vary
% with the overvoltage and the user would like to use these measurements.

% Here we assume that all SiPM measurement parameters are the same, which is typically true.
% Thus, we just put in the measured values for the first SiPM and at the
% completion of the for loop, provide the same values for all other SiPMs.
% However, if that is not the case, simply change the for loop index from 1
% to number of SiPMs.

for sipm_ind = 1:1,
    
    %     input_tmp = input('Please enter the input temperature at which the SiPM is operating: ','s');
    %     res_tmp = input('Do you have any measurements of SiPM parameter variation with temperature? Please type yes or no: ','s');
    %     if(strcmp(res_tmp,'yes'))
    %         %Please add the temperature varying component code here.
    %     end
    
    %     display('Reading the SiPM parameters at V_{ob} = V_{ob0} from file sipm_par_Vob0.dat')
    
    %Some initializations for variations with overvoltage
    
    pde_ratio = 1;
    ov_ratio_pde =1;
    ov_ratio_dc = 1;
    dc_ratio = 1;
    oct_ratio =1;
    ov_ratio_oct=1;
    ap_ratio=1;
    ov_ratio_ap=1;
    
    if(flag_default_variation_parameters_overvoltage)
        flags.overvoltage_pde = 1; 
        flags.overvoltage_pde_lin = 1; %PDE varies with overvoltage linearly
        flags.overvoltage_dc = 1; 
        flags.overvoltage_dc_lin = 1; % Dark current rate varies with overvoltage linearly
        flags.overvoltage_oct = 1;
        flags.overvoltage_oct_quad = 1; % Optical crosstalk varies with overvoltage quadratically
        flags.overvoltage_ap = 1;
        flags.overvoltage_ap_quad = 1; %Afterpulsing varies with overvoltage quadratically
    else
        
        display(' ');
        display('---------------------------');
        display(['           SiPM ' num2str(sipm_ind)]);
        display('---------------------------');
        display('Detector PDE variation with over-voltage');
        display('   ---------------------------');
        res = input('Do you have the measurements for PDE measurements with overvoltage (y/n): ','s');
        if(strcmp(res,'y'))
            flags.overvoltage_pde = 1;
            flags.overvoltage_pde_meas = 1;
            no_meas = input('Enter the number of overvoltage values at which the measurements are available: '); %Enter the value of PDE at V_{ob0} as the last entry
            %         no_meas = 5;
            filename = ['measurements\pde_ov_',num2str(sipm_ind),'.dat'];
            display(['Reading the measurements from file ' filename]);
            fid = fopen(filename,'r');
            if(fid == -1)
                error('File for PDE measurements vs overvoltage not found');
            end
            ov = fscanf(fid,'%f',no_meas); % The overvoltages written  as V_{ob}/V_{obz}.
            pde = fscanf(fid,'%f',no_meas); % The pde written as pde(V_{ob})/V_{obz}.
            ov_ratio_pde = ov/ov(no_meas);
            pde_ratio = pde/pde(no_meas);
            fclose(fid);
        else
            res = input('Would you like to use the usual linear relationship between overvoltage and PDE (y/n):','s');
            if(strcmp(res,'y'))
                flags.overvoltage_pde = 1;
                flags.overvoltage_pde_lin = 1;
            else
                flags.overvoltage_pde = 0;
            end
        end
        
        display(' ');
        display('Dark current variation with over-voltage');
        display('   ---------------------------');
        res = input('Do you have the measurements for dark current rate measurements with overvoltage (y/n): ','s');
        if(strcmp(res,'y'))
            flags.overvoltage_dc = 1;
            flags.overvoltage_dc_meas = 1;
            no_meas = input('Enter the number of overvoltage values at which the measurements are available: '); %Enter the value of PDE at V_{ob0} as the last entry
            filename = ['measurements\dc_ov_',num2str(sipm_ind),'.dat'];
            display(['Reading the measurements from file ' filename]);
            fid = fopen(filename,'r');
            if(fid == -1)
                error('File for dark current measurements vs overvoltage not found');
            end
            ov = fscanf(fid,'%f',no_meas); % The overvoltages written  as V_{ob}/V_{obz}.
            dc = fscanf(fid,'%f',no_meas); % The pde written as pde(V_{ob})/V_{obz}.
            ov_ratio_dc = ov/ov(no_meas);
            dc_ratio = dc/dc(no_meas);
            fclose(fid);
        else
            res = input('Would you like to use the usual linear relationship between overvoltage and dark current (y/n):','s');
            if(strcmp(res,'y'))
                flags.overvoltage_dc = 1;
                flags.overvoltage_dc_lin = 1;
            else
                flags.overvoltage_dc = 0;
            end
        end
        
        display(' ');
        display('Optical Crosstalk probability variation with over-voltage');
        display('   ---------------------------');
        res = input('Do you have the measurements for optical cross talk probability with overvoltage (y/n): ','s');
        if(strcmp(res,'y'))
            flags.overvoltage_oct = 1;
            flags.overvoltage_oct_meas = 1;
            no_meas = input('Enter the number of overvoltage values at which the measurements are available: ');
            filename = ['measurements\pde_oct_',num2str(sipm_ind),'.dat'];
            display(['Reading the measurements from file ' filename]);
            fid = fopen('filename','r');
            ov = fscanf(fid,'%f',no_meas);
            Nbar_oct = fscanf(fid,'%f',no_meas);
            fclose(fid);
            ov_ratio_oct = ov/ov(no_meas);
            oct_ratio = Nbar_oct/Nbar_oct(no_meas);
        else
            res = input('Would you like to use the usual quadratic relationship between overvoltage and after pulsing probability (y/n):','s');
            if(strcmp(res,'y'))
                flags.overvoltage_oct = 1;
                flags.overvoltage_oct_quad = 1;
            else
                flags.overvoltage_oct = 0;
            end
        end
        
        display(' ');
        display('Afterpulsing probability variation with over-voltage');
        display('   ---------------------------');
        res = input('Do you have the measurements for afterpulsing probability with overvoltage (y/n): ','s');
        if(strcmp(res,'y'))
            
            flags.overvoltage_ap = 1;
            flags.overvoltage_ap_meas = 1;
            filename = ['measurements\pde_ap_',num2str(sipm_ind),'.dat'];
            display(['Reading the measurements from file ' filename]);
            fid = fopen('filename','r');
            ov = fscanf(fid,'%f',no_meas);
            Nbar_ap = fscanf(fid,'%f',no_meas);
            fclose(fid);
            ov_ratio_ap = ov/ov(no_meas);
            ap_ratio = Nbar_ap/Nbar_ap(no_meas);
        else
            res = input('Would you like to use the usual quadratic relationship between overvoltage and after pulsing probability (y/n):','s');
            if(strcmp(res,'y'))
                flags.overvoltage_ap = 1;
                flags.overvoltage_ap_quad = 1;
            else
                flags.overvoltage_ap = 0;
            end
        end
        
        disp('');
        disp('');
    end
    
    meas_sipm_par(sipm_ind).qe = qe(sipm_ind);
    meas_sipm_par(sipm_ind).Nbar_oct_mean = Nbar_oct_mean(sipm_ind);
    meas_sipm_par(sipm_ind).Nbar_oct_var = Nbar_oct_var(sipm_ind);
    meas_sipm_par(sipm_ind).Nbar_ap_short_mean = Nbar_ap_short_mean(sipm_ind);
    meas_sipm_par(sipm_ind).Nbar_ap_short_var = Nbar_ap_short_var(sipm_ind);
    meas_sipm_par(sipm_ind).Nbar_ap_long_mean = Nbar_ap_long_mean(sipm_ind);
    meas_sipm_par(sipm_ind).Nbar_ap_long_var = Nbar_ap_long_var(sipm_ind);
    meas_sipm_par(sipm_ind).G_mean = G_mean(sipm_ind);
    meas_sipm_par(sipm_ind).G_var = G_var(sipm_ind);
    meas_sipm_par(sipm_ind).rdc_mean = rdc_mean(sipm_ind);
    meas_sipm_par(sipm_ind).rdc_var = rdc_var(sipm_ind);
    meas_sipm_par(sipm_ind).tau_phd1_mean = tau_phd1_mean(sipm_ind);
    meas_sipm_par(sipm_ind).tau_phd1_var = tau_phd1_var(sipm_ind);
    meas_sipm_par(sipm_ind).tau_ap_short_mean = tau_ap_short_mean(sipm_ind);
    meas_sipm_par(sipm_ind).tau_ap_short_var = tau_ap_short_var(sipm_ind);
    meas_sipm_par(sipm_ind).tau_ap_long_mean = tau_ap_long_mean(sipm_ind);
    meas_sipm_par(sipm_ind).tau_ap_long_var = tau_ap_long_var(sipm_ind);
    meas_sipm_par(sipm_ind).V_ob_mean = Vob_mean(sipm_ind);
    meas_sipm_par(sipm_ind).V_ob_var = Vob_var(sipm_ind);
    meas_sipm_par(sipm_ind).Rq_mean = Rq_mean(sipm_ind);
    meas_sipm_par(sipm_ind).Rq_var = Rq_var(sipm_ind);
    meas_sipm_par(sipm_ind).pde_ratio = pde_ratio;
    meas_sipm_par(sipm_ind).ov_ratio_pde = ov_ratio_pde;
    meas_sipm_par(sipm_ind).dc_ratio = dc_ratio;
    meas_sipm_par(sipm_ind).ov_ratio_dc = ov_ratio_dc;
    meas_sipm_par(sipm_ind).oct_ratio = oct_ratio;
    meas_sipm_par(sipm_ind).ov_ratio_oct = ov_ratio_oct;
    meas_sipm_par(sipm_ind).ap_ratio = ap_ratio;
    meas_sipm_par(sipm_ind).ov_ratio_ap = ov_ratio_ap;
    meas_sipm_par(sipm_ind).f = f;
    
end

% Now let all the SiPMs have the same property that the first SiPM has
for sipm_ind = 2:sipm_no,
    meas_sipm_par(sipm_ind) = meas_sipm_par(1);
end
