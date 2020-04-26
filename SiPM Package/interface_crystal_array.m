%This is the main routine. A detailed description of this software is
% provided in the following manuscript
% Jha AK, van Dam HT, Kupinski MA, Clarkson E. Simulating Silicon
% Photomultiplier Response to Scintillation Light. IEEE Trans Nucl Sci. 2013
% Feb;30(1):336-351. PubMed PMID: 26236040; PubMed Central PMCID: PMC4519993.
% Please also refer to this manuscript in any publications that use this
% code.
% For any questions, please email: abhinavkumar.jha@gmail.com
% Copyright: University of Arizona

clear all;
close all;
no_timesteps = 50;
sipm_no = 16;

[meas_sipm_par(1:16), fixed_sipm_par,crystal_par,flags] = read_input_array();

E = 511; % Energy of the gamma-ray photon in keV.
Nbar_ph = crystal_par.Qsc*E;

det_pde = zeros(2,1)``; %The SiPM photon-detection efficiency

for i = 1:sipm_no,
    det_pde(i) = meas_sipm_par(i).qe*meas_sipm_par(i).f;
end

op_charge = zeros(sipm_no,1);
% matlabpool close
matlab_version_flag = 0; %Set this flag to 1 if this software is being executed using a
%Matlab version where the matlabpool command is
%still active.

if(matlab_version_flag)
    matlabpool close force local
    
    if ~(matlabpool('size'))
        matlabpool open
    end
end

display 'Starting the simulation';

parfor sipm_ind = 1:sipm_no,
    
    Nph_t = get_timestamp_photons(Nbar_ph, crystal_par.Tr, crystal_par.Td,fixed_sipm_par.tmeas);
    %These are the time stamps of the optical photons emitted by the
    % scintillation crystal. These time stamps could also be obtained from other more sophisticated
    %studies, such as a GEANT-based study.
    
    tic;
    [op_charge(sipm_ind)] =...
        sipm_sim_array(no_timesteps,Nph_t,meas_sipm_par(sipm_ind),fixed_sipm_par,crystal_par,flags);
    toc;
    
    display(['SiPM ' num2str(sipm_ind) ' output in Coulumbs is ' num2str(op_charge(sipm_ind)) ]);
end

display 'Ending the simulation'


