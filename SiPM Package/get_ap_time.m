
function [Nap_t] = get_ap_time(tau_ap_short,tau_ap_long,Nap_short,Nap_long,Nap)
%This routine provides the time stamps for the afterpulse photons. 

    t = linspace(0,6*max(tau_ap_long,tau_ap_short), 1000);
    u = 1 - Nap_short/(Nap_short + Nap_long)*exp(-t/tau_ap_short) + Nap_long/(Nap_short + Nap_long)* exp(-t/tau_ap_long);
    
    Nap_t = interp1(u,t,rand(Nap,1),'splines');
    
    