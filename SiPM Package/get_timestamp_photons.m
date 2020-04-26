
function [Nph_t] = get_timestamp_photons(Nbar_ph, Tr, Td,tmeas)
%This is a simple routine to evaluate the timestamp of the optical photons.
% The timestamps can also be obtained more rigorously using advanced
% software such as GEANT4. 

    Nph = poissrnd(Nbar_ph,1); % Number of photons generated.
    % Times at which the photons are generated.
    t = linspace(0,tmeas,1000);
    u = 1 - Td/(Td - Tr)*exp(-t/Td) + Tr/(Td - Tr)* exp(-t/Tr);
    
    %Nph_t = -Td* log((1-rand(Nph,1))*(1 -Tr/Td));
    Nph_t = interp1(u,t,rand(Nph,1),'splines');
    
    