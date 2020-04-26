function [ph_fired] = sipm_sim_array(no_timesteps,Nph_t,meas_sipm_par,fixed_sipm_par, crystal_par, flags )
% The code to simulate the SiPM operation.
% The following are the input arguments:
% no_timesteps: The number of time steps into which you would like to
% divide the total measurement time
% Nph_t: The timestamps of the scintillation photons
% meas_sipm_par: The measured parameters of the SiPM. Please see
% read_input_array.m
% fixed_sipm_par: The fixed SiPM parameters. Again see read_input_array.m
% crystal_par: The crsytal decay and rise times. Please see read_input_array.m
% flags: The flags to turn on and off the different SiPM mechanisms.
% Default value is 1. Please see read_input_array.m for more details.

% Read in all the input variables to the local variables in this code
% % First the measured parameters of the SiPM
qe = meas_sipm_par.qe;
Nbar_oct_mean = meas_sipm_par.Nbar_oct_mean ;
Nbar_oct_var = meas_sipm_par.Nbar_oct_var;
Nbar_ap_short_mean = meas_sipm_par.Nbar_ap_short_mean;
Nbar_ap_short_var = meas_sipm_par.Nbar_ap_short_var;
Nbar_ap_long_mean = meas_sipm_par.Nbar_ap_long_mean;
Nbar_ap_long_var = meas_sipm_par.Nbar_ap_long_var;
G_mean = meas_sipm_par.G_mean;
G_var = meas_sipm_par.G_var;
rdc_mean = meas_sipm_par.rdc_mean;
rdc_var = meas_sipm_par.rdc_var;
V_ob_mean = meas_sipm_par.V_ob_mean;
V_ob_var = meas_sipm_par.V_ob_var;
Rq_mean = meas_sipm_par.Rq_mean;
Rq_var = meas_sipm_par.Rq_var;

f = meas_sipm_par.f;

% The overvoltage variations
if(flags.overvoltage_pde && flags.overvoltage_pde_meas)
    pde_ratio = meas_sipm_par.pde_ratio;
    ov_ratio_pde = meas_sipm_par.ov_ratio_pde;
end

if(flags.overvoltage_dc && flags.overvoltage_dc_meas)
    dc_ratio = meas_sipm_par.dc_ratio;
    ov_ratio_dc = meas_sipm_par.ov_ratio_dc;
end

if(flags.overvoltage_oct && flags.overvoltage_oct_meas)
    oct_ratio = meas_sipm_par.oct_ratio;
    ov_ratio_oct = meas_sipm_par.ov_ratio_oct;
end

if(flags.overvoltage_ap && flags.overvoltage_ap_meas)
    ap_ratio = meas_sipm_par.ap_ratio;
    ov_ratio_ap = meas_sipm_par.ov_ratio_ap;
end

% Now the fixed parameters of the SiPM
Ncells = fixed_sipm_par.Ncells;
Rs = fixed_sipm_par.Rs;
sipm_pix_xno = fixed_sipm_par.sipm_pix_xno;
sipm_pix_yno = fixed_sipm_par.sipm_pix_yno;
tmeas = fixed_sipm_par.tmeas;

% Finally read the crystal parameters
Td = crystal_par.Td;
Qsc = crystal_par.Qsc;

oct_delay = 1e-9;
charge_elec = 1.6e-19;

if (flags.output_voltage == 1 && no_timesteps > 1)
    num_time_samples = no_timesteps;
    tmp = linspace(0,1 - exp(-6),num_time_samples);
    Tsamp = -Td * log(1 - tmp);
    if(Tsamp(num_time_samples) > tmeas)
        Tsamp(num_time_samples)= tmeas;
    end
    %    Tsamp = linspace(0,tmeas,num_time_samples);
    delt = Tsamp(2:num_time_samples) - Tsamp(1:num_time_samples-1);
    delt(num_time_samples) = tmeas - Tsamp(num_time_samples);
else
    num_time_samples = 1; %ceil(tmeas/tau_phd1);
    Tsamp = 0;
    delt = tmeas;
end

mcell(sipm_pix_xno,sipm_pix_yno,num_time_samples) = struct('count',0,'t',[],'gain',[],'flag',[]);

mcell_hit(sipm_pix_xno,sipm_pix_yno) = struct('count',0,'t',[],'flag',[]);
phot_pot = struct('count',0,'t',[], 'flag',[],'i',[],'j',[],'gain',[],'analysed_count',0);

Cd = 77e-15;
Cq = 23e-15;
Cg = 28e-12;
Rq = Rq_mean;
Rd = 1e3;
Rqp = Rq/(Ncells-1);
Cqp = Cq*(Ncells-1);
Cdp = Cd*(Ncells-1);

qe_effect = 1;
dc_effect = 1;
oct_effect = 1;
ap_effect = 1;


lost_phot = 0;
%% Running the simulation


mcell_hit(sipm_pix_xno,sipm_pix_yno) = struct('count',0,'t',[],'flag',[]);
mcell(sipm_pix_xno,sipm_pix_yno,num_time_samples) = struct('count',0,'t',[],'gain',[],'flag',[]);
phot_pot = struct('count',0,'t',[], 'flag',[],'i',[],'j',[],'gain',[],'analysed_count',0);

Nbar_oct = Nbar_oct_mean + randn(1)*Nbar_oct_var ; %Many Stochastic here
Nbar_ap_short = Nbar_ap_short_mean + randn(1)*Nbar_ap_short_var;
Nbar_ap_long = Nbar_ap_long_mean + randn(1)*Nbar_ap_long_var;
Nbar_ap =  Nbar_ap_short + Nbar_ap_long;
G = G_mean + randn(1)*G_var;
rdc = rdc_mean + randn(1)*rdc_var;
tau_ap_short = meas_sipm_par.tau_ap_short_mean + randn(1)*meas_sipm_par.tau_ap_short_var ;
tau_ap_long = meas_sipm_par.tau_ap_long_mean + randn(1)*meas_sipm_par.tau_ap_long_var ;
V_ob = V_ob_mean + randn(1)*V_ob_var;
Rq = Rq_mean + randn(1)*Rq_var;

tau1 = charge_elec * G*Rq/V_ob;
tau2 = Rq*Cq;
tau3 = Rs*Ncells*Cd;
tau4 = Rs*Cg;

taumr = 1e-12;
taumd = Cd*(Rd*Rq/(Rd+Rq));
tau_aval = taumd;

taucd1 = ((tau1 + tau3 + tau4 ) + sqrt( (tau1 + tau3 + tau4 )^2 - 4 * ( tau1 * tau4 + tau2 * tau3))) / 2;
taucd2 = ((tau1 + tau3 + tau4 ) - sqrt( (tau1 + tau3 + tau4 )^2 - 4 * ( tau1 * tau4 + tau2 * tau3))) / 2;

a1 = taucd1*(taucd1 - tau2)/((taucd1 - taucd2) * ( taucd1 - taumd) * ( taucd1 - taumr));
a2 = taucd2*(taucd2 - tau2)/((taucd2 - taucd1) * ( taucd2 - taumd) * ( taucd2 - taumr));
a3 = taumd*(taumd - tau2)/((taumd - taucd1) * ( taumd - taucd2) * ( taumd - taumr));
a4 = taumr*(taumr - tau2)/((taumr - taucd1) * ( taumr - taucd2) * ( taumr - taumd));

b1 = tau1/( (tau1 - taumr)*(tau1 - taumd));
b2 = taumd/( (taumd - tau1)*(taumd - taumr));
b3 = taumr/( (taumr - tau1)*(taumr - taumd));
b4 = (taucd1 - tau2)*taucd1/( (taucd1 - taucd2)*(taucd1 - tau1)*(taucd1 - taumd)*(taucd1 - taumr));
b5 = (taucd2 - tau2)*taucd2/( (taucd2 - taucd1)*(taucd2 - tau1)*(taucd2 - taumd)*(taucd2 - taumr));
b6 = (taumd - tau2)*taumd/( (taumd - taucd1)*(taumd - tau1)*(taumd - taucd2)*(taumd - taumr));
b7 = (taumr - tau2)*taumr/( (taumr - taucd1)*(taumr - tau1)*(taumr - taumd)*(taumr - taucd2));
b8 = (tau1 - tau2)*tau1/( (tau1 - taucd1)*(tau1 - taucd2)*(tau1 - taumd)*(tau1 - taumr));


for i = 1:sipm_pix_xno,
    for j = 1:sipm_pix_yno,
        mcell_hit(i,j).count = 0;
        mcell_hit(i,j).t = [];
        for k = 1:num_time_samples,
            mcell(i,j,k).count = 0;
            mcell(i,j,k).t = [];
        end
    end
end

Nph = length(Nph_t);
% Now assume that when the photons travel to the SiPM, the order between them is maintained.

% Of the generated photons, some photons reach the SiPM and some do not.
coll_success = rand(Nph,1);
ind = (coll_success < f);
Nph_sipm = Nph_t(ind);
Nph_sipm_no = length(Nph_sipm);


Nph_sipm = sort(Nph_sipm,'ascend');
phot_rej = 0;


for ph_ind = 1:Nph_sipm_no,
    i = ceil(rand(1)*sipm_pix_xno);
    j = ceil(rand(1)*sipm_pix_yno);
    
    if(rand(1) < qe)
        if (mcell_hit(i,j).count > 0 && flags.deadtime )
            if (Nph_sipm(ph_ind) - mcell_hit(i,j).t(mcell_hit(i,j).count) > tau_aval)
                mcell_hit(i, j).t = [mcell_hit(i, j).t Nph_sipm(ph_ind)];
                mcell_hit(i, j).count = mcell_hit(i, j).count + 1;
                mcell_hit(i,j).flag = [mcell_hit(i, j).flag 1];
                phot_pot.t = [phot_pot.t Nph_sipm(ph_ind)];
                phot_pot.flag = [phot_pot.flag 1];
                phot_pot.i = [phot_pot.i i];
                phot_pot.j = [phot_pot.j j];
                phot_pot.gain = [phot_pot.gain 0];
                phot_pot.count = phot_pot.count + 1;
            else
                phot_rej = phot_rej + 1;
            end
        elseif (~flags.deadtime || mcell_hit(i,j).count == 0)
            mcell_hit(i, j).t = [mcell_hit(i, j).t Nph_sipm(ph_ind)];
            mcell_hit(i, j).count = mcell_hit(i, j).count + 1;
            phot_pot.t = [phot_pot.t Nph_sipm(ph_ind)];
            phot_pot.flag = [phot_pot.flag 1];
            phot_pot.i = [phot_pot.i i];
            phot_pot.j = [phot_pot.j j];
            phot_pot.gain = [phot_pot.gain 0];
            phot_pot.count = phot_pot.count + 1;
            
        end
    end
end


if(flags.dc)
    % The times at which the dark current will trigger avalanches
    for i = 1:sipm_pix_xno,
        for j = 1:sipm_pix_yno,
            prob_darkcurr = rdc*tmeas/Ncells;
            Ndc = poissrnd(prob_darkcurr); %Stochastic
            %Ndc = uint8(prob_darkcurr);
            if(Ndc > 0)
                darkcurr_time = rand(Ndc,1)*tmeas;
                mcell_hit(i,j).t = [mcell_hit(i, j).t darkcurr_time'];
                mcell_hit(i,j).count = mcell_hit(i, j).count + Ndc;
                phot_pot.t = [phot_pot.t darkcurr_time'];
                phot_pot.flag = [phot_pot.flag 2*ones(Ndc,1)'];
                phot_pot.i = [phot_pot.i i*ones(Ndc,1)'];
                phot_pot.j = [phot_pot.j j*ones(Ndc,1)'];
                phot_pot.gain = [phot_pot.gain 0*ones(Ndc,1)'];
                phot_pot.count = phot_pot.count + Ndc;
            end
        end
    end
end


for i = 1:sipm_pix_xno,
    for j = 1:sipm_pix_yno,
        [mcell_hit(i,j).t] = sort(mcell_hit(i,j).t,'ascend');
        %             mcell_hit(i,j).flag = mcell_hit(i,j).flag(order);
    end
end

[phot_pot.t,order] = sort(phot_pot.t, 'ascend');
phot_pot.flag = phot_pot.flag(order);
phot_pot.i = phot_pot.i(order);
phot_pot.j = phot_pot.j(order);
phot_pot.gain = phot_pot.gain(order);

% Now that we know which photons are hitting which microcell, and since
% the various parameters of the SiPM are themselves function of
% over-voltage which in turn is a function of output voltage, therefore, we need to sample our time interval so that
% at the end of different time intervals, we evaluate the  output voltage
% and use that to redetermine the various SiPM parameters.

Vs = 0;
ind_t = 1;
phot_pot_samp_count = zeros(num_time_samples,1);
phot_glob_count = 0;
%     ap_count_2 = 0;
%     oct_count_2 = 0;
I_Rs = 0;

for t_samp = Tsamp,
    
    
    if(flags.output_voltage && ind_t > 1)
        Vs = (Rs) * I_Rs ; %(charge_time_int(ind_t-1)) * G * charge_elec/delt(ind_t-1);
        Vs_ind(ind_t) = Vs;
    end
    
    for i = 1:sipm_pix_xno,
        for j = 1:sipm_pix_yno,
            for k = 1:mcell_hit(i,j).count,
                if(mcell_hit(i,j).t(k) >= t_samp && mcell_hit(i,j).t(k) < t_samp + delt(ind_t)),
                    mcell(i,j,ind_t).t = [mcell(i, j,ind_t).t mcell_hit(i,j).t(k)];
                    mcell(i,j,ind_t).count = mcell(i, j,ind_t).count + 1;
                    mcell(i,j,ind_t).gain = [mcell(i, j,ind_t).gain 0];
                    %                         mcell(i,j,ind_t).flag = [mcell(i, j,ind_t).flag mcell_hit(i,j).flag];
                    phot_pot_samp_count(ind_t) = phot_pot_samp_count(ind_t) + 1;
                end
            end
        end
    end
    
    for i = 1:sipm_pix_xno,
        for j = 1:sipm_pix_yno,
            [mcell(i,j,ind_t).t, order] = sort(mcell(i,j,ind_t).t,'ascend');
        end
    end
    
    gain = 1;
    phot_count = 0;
    
    if(phot_pot_samp_count(ind_t) > 0)
        photon_present = 1;
    else
        photon_present = 0;
    end
    
    while photon_present,
        phot_count = phot_count + 1;
        phot_glob_count = phot_glob_count + 1;
        
        i = phot_pot.i(phot_glob_count);
        j = phot_pot.j(phot_glob_count);
        
        for k = 1:mcell(i,j,ind_t).count,
            if( mcell(i,j,ind_t).t(k) == phot_pot.t(phot_glob_count))
                break;
            end
        end
        
        if(isempty(k))
            display('break here due to empty k');
        end
        if(flags.partial_recovery || flags.over_voltage_effect)
            
            
            I_Rq = 0; % Current through quenching resistance due to all the previous avalanches
            if(ind_t >1)
                for prev_ind_t = 1:ind_t-1
                    for prev_phot_ind = 1:mcell(i,j,prev_ind_t).count
                        prev_trigger_time = mcell(i,j,prev_ind_t).t(prev_phot_ind);
                        tp = mcell(i,j,ind_t).t(k) - prev_trigger_time;
                        I_Rq = I_Rq + (G*charge_elec) * mcell(i,j,prev_ind_t).gain(prev_phot_ind)* ...
                            ((b1 * exp(-tp/tau1) + b2*exp(-tp/taumd) + b3*exp(-tp/taumr)) + ...
                            (tau3/Ncells)*(b4*exp(-tp/taucd1) + b5*exp(-tp/taucd2) + b6*exp(-tp/taumd) + b7*exp(-tp/taumr)+ b8*exp(-tp/tau1)));
                    end
                end
            end
            if( k >1)
                for prev_phot_ind = 1:k-1
                    prev_trigger_time = mcell(i,j,ind_t).t(prev_phot_ind);
                    tp = mcell(i,j,ind_t).t(k) - prev_trigger_time;
                    I_Rq = I_Rq + (G*charge_elec) * mcell(i,j,ind_t).gain(prev_phot_ind)* ...
                        ((b1 * exp(-tp/tau1) + b2*exp(-tp/taumd) + b3*exp(-tp/taumr)) - ...
                        (tau3/Ncells)*(b4*exp(-tp/taucd1) + b5*exp(-tp/taucd2) + b6*exp(-tp/taumd) + b7*exp(-tp/taumr) + b8*exp(-tp/tau1)));
                    
                end
            end
            
            
            gain = (1 - Vs/V_ob - I_Rq*Rq/V_ob);
            if(gain < 0)
                gain = 0;
                lost_phot = lost_phot+1 ;
            end
            if(isnan(gain))
                display('gain is nan');
            end
            if(gain > 1)
                gain = 0;
            end
        else
            gain = 1;
        end
        
        if(isnan(gain))
            display('gain is nan');
        end
        
        
        if(flags.overvoltage)
            if(flags.overvoltage_pde_lin)
                qe_effect = gain;
            else if(flags.overvoltage_pde_meas)
                    qe_effect = interp1(ov_ratio_pde,pde_ratio,gain,'splines');
                else
                    qe_effect = 1;
                end
            end
            if(flags.dc)
                if(flags.overvoltage_dc && flags.overvoltage_dc_lin)
                    dc_effect = gain;
                else if(flags.overvoltage_dc && flags.overvoltage_dc_meas)
                        dc_effect = interp1(ov_ratio_dc, dc_ratio,gain,'splines');
                    else
                        dc_effect = 1;
                    end
                end
            end
            if(flags.oct)
                if(flags.overvoltage_oct && flags.overvoltage_oct_quad)
                    oct_effect = gain^2;
                else if(flags.overvoltage_oct && flags.overvoltage_ap_meas)
                        oct_effect = interp1(ov_ratio_oct, oct_ratio,gain,'splines');
                    else
                        oct_effect = 1;
                    end
                end
            end
            if(flags.ap)
                if(flags.overvoltage_ap && flags.overvoltage_ap_quad)
                    ap_effect = gain^2;
                else if(flags.overvoltage_ap && flags.overvoltage_ap_meas)
                        ap_effect = interp1(ov_ratio_ap, ap_ratio,gain,'splines');
                    else
                        ap_effect = 1;
                    end
                end
            end
        end
        sel_rand = rand(1);
        
        if( ( sel_rand < qe_effect && (phot_pot.flag(phot_glob_count) == 1))  || ...
                (rand(1) < dc_effect && (phot_pot.flag(phot_glob_count) == 2))  ||...
                phot_pot.flag(phot_glob_count) == 3 || phot_pot.flag(phot_glob_count) == 4)
            % Check if the photon or dark current or
            % afterpulses ( due to previous photon triggered avalanches) or OCT triggers an avalanche. If so, only
            % then is there scope for OCT, afterpulse etc.
            % Afterpulses do not require the gain parameter
            % to be evaluated for a different overvoltage. We have
            % already considered that when evaluating prob(N_ap).
            if(flags.partial_recovery)
                phot_pot.gain(phot_glob_count) = gain;
            end
            mcell(i,j,ind_t).gain(k) = gain;
            
            
            Nbar_ap_d = Nbar_ap*ap_effect;
            Nbar_oct_d = Nbar_oct*oct_effect;
            
            
            % The times at which the after pulses occur. The afterpulses are
            % inhomogeneous Poisson point processes with intensity function being
            % an exponential with decay time tap, and the mean of the Poisson being
            % N_ap. Therefore, this event will occur only
            % after the avalanche that triggers it, and will be
            % following it in the list.
            
            if(flags.ap)
                Nap = poissrnd(Nbar_ap_d); % Stochastic
                %Nap = uint8(Nbar_ap_d);
                if(Nap > 0)
                    ap_time = get_ap_time(tau_ap_short, tau_ap_long, Nbar_ap_short, Nbar_ap_long,Nap); %exprnd(tau_ap,Nap,1);
                    for kp = 1:Nap,
                        ap_time_ind = 0;
                        time_ap = ap_time(kp) + mcell(i,j,ind_t).t(k);
                        for tmp_ind_t = ind_t+1:num_time_samples,
                            if(time_ap < Tsamp(tmp_ind_t))
                                ap_time_ind = tmp_ind_t-1;
                                break;
                            end
                        end
                        if(time_ap > Tsamp(num_time_samples) && time_ap < tmeas)
                            ap_time_ind = num_time_samples;
                        end
                        if(ap_time_ind)
                            mcell(i,j,ap_time_ind).t = [mcell(i,j,ap_time_ind).t ap_time(kp) + mcell(i,j,ind_t).t(k)];
                            mcell(i,j,ap_time_ind).count = mcell(i,j,ap_time_ind).count + 1;
                            mcell(i,j,ap_time_ind).gain = [mcell(i,j,ap_time_ind).gain 0];
                            %                                 mcell(i,j,ap_time_ind).flag = [mcell(i,j,ap_time_ind).flag 3];
                            %ap_count_2 = ap_count_2 + 1;
                            % Re-sort to get the afterpulses at the correct place in
                            % the list.
                            [mcell(i,j,ap_time_ind).t, order] = sort(mcell(i,j,ap_time_ind).t,'ascend');
                            phot_pot.t = [phot_pot.t time_ap];
                            phot_pot.flag = [phot_pot.flag 3];
                            phot_pot.i = [phot_pot.i i];
                            phot_pot.j = [phot_pot.j j];
                            phot_pot.gain = [phot_pot.gain 0];
                            phot_pot.count = phot_pot.count + 1;
                            
                            phot_pot_samp_count(ap_time_ind) = phot_pot_samp_count(ap_time_ind) + 1;
                            
                            
                        end
                    end
                end
            end
            
            % The OCT events occur within 1 ns in the neighboring 4
            % cells ( with uniform probability).
            
            if(flags.oct)
                Noct = poissrnd(Nbar_oct_d); %Stochastic
                %Noct = uint8(Nbar_oct_d);
                for mp = 1:Noct,
                    nind = ceil(rand(1)*4);
                    if(nind == 1)
                        ip = i;
                        jp = j + 1;
                    else if(nind == 2)
                            ip = i + 1;
                            jp = j;
                        else if(nind == 3)
                                ip = i - 1;
                                jp = j;
                            else if ( nind == 4)
                                    ip = i;
                                    jp = j -1;
                                end
                            end
                        end
                    end
                    time_oct = mcell(i, j,ind_t).t(k) + ( mp/10)*oct_delay;
                    if(ip > 0 && ip <= sipm_pix_xno && jp > 0 && jp <= sipm_pix_yno && time_oct < tmeas)
                        ind_ttmp = ind_t;
                        
                        for tmp_ind_t = ind_t+1:num_time_samples,
                            if( time_oct < Tsamp(tmp_ind_t))
                                ind_ttmp = tmp_ind_t-1;
                                break;
                            end
                        end
                        if(time_oct > Tsamp(num_time_samples) && time_oct < tmeas)
                            ind_ttmp = num_time_samples;
                        end
                        
                        mcell(ip,jp,ind_ttmp ).t = [mcell(ip,jp,ind_ttmp).t , mcell(i, j,ind_t).t(k) + ( mp/10)*oct_delay];
                        mcell(ip,jp,ind_ttmp ).count = mcell(ip,jp,ind_ttmp ).count + 1;
                        mcell(ip,jp,ind_ttmp).gain = [mcell(ip,jp,ind_ttmp).gain 0];
                        [mcell(ip,jp,ind_ttmp).t, order] = sort(mcell(ip,jp,ind_ttmp).t,'ascend');
                        %                                 mcell(ip,
                        %                                 jp,ind_t).flag = mcell(ip,jp,ind_t).flag(order);
                        %                             mcell(ip,jp,ind_ttmp ).flag = [mcell(ip,jp,ind_ttmp ).flag 4];
                        %                                 kp = mcell(ip,jp,ind_t).count;
                        
                        phot_pot.t = [phot_pot.t time_oct];
                        phot_pot.flag = [phot_pot.flag 4];
                        phot_pot.i = [phot_pot.i ip];
                        phot_pot.j = [phot_pot.j jp];
                        phot_pot.count = phot_pot.count + 1;
                        phot_pot.gain = [phot_pot.gain 0];
                        phot_pot_samp_count(ind_ttmp) = phot_pot_samp_count(ind_ttmp) + 1;              
                    end
                end
            end
        end
        
        [phot_pot.t,order] = sort(phot_pot.t, 'ascend');
        phot_pot.flag = phot_pot.flag(order);
        phot_pot.i = phot_pot.i(order);
        phot_pot.j = phot_pot.j(order);
        phot_pot.gain = phot_pot.gain(order);
        
        
        
        if(phot_count == phot_pot_samp_count(ind_t)) % All the photons have been processed. Using a while loop due to afterpulsing and crosstalk, which can
            %create new photons within the same time interval as well
            photon_present = 0;
            phot_pot.analysed_count = phot_pot.analysed_count + phot_pot_samp_count(ind_t);
        end
    end
    
    I_Rs = 0;
    if(flags.output_voltage)
        if(ind_t < num_time_samples)
            for k = 1:phot_pot.analysed_count,
                tp = t_samp + delt(ind_t) - phot_pot.t(k);
                I_Rs = I_Rs + G*(charge_elec) * phot_pot.gain(k) * ...
                    ( a1*exp(-tp/taucd1) + a2* exp(-tp/taucd2) +  a3 *exp(-tp/taumd) + a4 * exp(-tp/taumr));
            end
        end
    end
    %     charge_time_int(ind_t) = I_Rs*delt(ind_t)/(G*charge_elec);
    
    ind_t = ind_t+1;
end

% Neq_fired = sum(charge_time_int);
ph_fired = sum( phot_pot.gain)*G*charge_elec;

