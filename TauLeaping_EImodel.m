function [spike_times,spike_ids,network_state,Irc,time_points] = TauLeaping_EImodel(W,response_fn,beta,alpha,I,t_min,t_max,init_state,args)
% This is an implementation of Tau Leaping Algorithm for the analysis
% of neuronal avalanches in Wilson-Cowan based stochastic
% population dynamics. Tau Leaping belongs to the family of
% stochastic simulation algoritmhs that are deriven from the
% Master Equation of an ODE. Tau Leaping, just like SSA,
% provides a compromise between continuous-deterministic solutions
% obtainable by direct solutions of ODEs and deterministic-stochastic
% solutions obtainable by exact/approximate solutions of the
% Master Equation. It is possible to interpret stochastic simulators
% also as an extension of systems with Langevin Dynamics. For more
% information, refer to the 2006 paper by Gillespie.
%
% Inputs:
% - W [N,N]: is the weight matrix,  W(i,j) is synaptic weight
% to the ith neuron FROM the jth
% - response_fn [1.N]: handle for response function (sigmoid, 
% hyperbolic tan, etc.)
% - input [N,T] is the net temporal input
% - alpha [1,N]: is the rate at which active neurons decay to being
% quiescent
% - beta [1,N]: is the height of the response function, 
% i.e. the rate at
% which saturated-input quiescent neurons become active
% - init_state [1,N]: is the initial state vector
%
% Method Specific Inputs:
% - mode {0,1,2}: is the operation mode of Tau Leaping, 
% defines how spike times are spread around each leap interval.
%
% Outputs:
% - spike_times [1,E(a)]: timing of each neuronal spike according 
% to the chosen mode
% - spike_ids [1,E(a)]: index of spiking neurons for each leap interval
% - Irc [N,T]: reconstruction of the temporal neuronal inputs 
% from the input at each leap interval
% - tc [1,T]: reconstruction of input timesteps
%
% References:
% Citation: Benayoun M, Cowan JD, van Drongelen W, Wallace E (2010) 
% Avalanches in a Stochastic Model of Spiking Neurons. 
% PLoS Comput Biol 6(7): e1000846. 
% https://doi.org/10.1371/journal.pcbi.1000846
% 
% Citation Yang Cao, Daniel T. Gillespie, Linda R. Petzold (2006); 
% Efficient step size selection for the tau-leaping simulation method. 
% J. Chem. Phys. 28 January 2006; 124 (4): 044109. 
% https://doi.org/10.1063/1.2159468

val = values(args);
mode = val(1);

N = size(W,1);
Tau = 0.5;

time_points = t_min:Tau:t_max;
active = init_state(1,:)';
network_states = zeros(N, length(time_points));
event_count = zeros(N, length(time_points));
Irc = [];
spike_times = [];
spike_ids = [];

network_states(:,1) = active;
event_count(:,1) = zeros(N,1);
leap_no = 0;

current = W*active + I(:,1);
propensity = beta .* ~active .* feval(response_fn, current) + ...
    alpha .* active;

for t=t_min:Tau:t_max
    num_events = poissrnd(propensity * Tau);
    
    for i = 1:N
        if num_events(i) > 0 
            if mod(num_events(i),2)==1
                if active(i) == 1
                    active(i) = 0;
                    current = current - W(:,i);
                else
                    active(i) = 1;
                    current = current + W(:,i);
                end
            else
                continue;
            end
        end
    end
    leap_no = leap_no + 1;

    Ipast = I(:,floor(t+1));
    Icurr = I(:,floor(t+Tau+1));

    current = current + (Icurr - Ipast);
    propensity = beta .* ~active .* feval(response_fn, current) + ...
        alpha .* active;

    network_states(:,leap_no) = active;
    event_count(:,leap_no) = num_events;
    Irc = [Irc Icurr];
end

for i=1:size(network_states,2)
    if mode == 0 % DOES NOT WORK - NEED TO CHANGE BIN STRATEGY
        add = repmat(time_points(i),1,sum(network_states(:,i)));
    elseif mode == 1
        add = linspace(time_points(i), time_points(i)+Tau, sum(network_states(:,i)));
    elseif mode == 2
        add = linspace(time_points(i), time_points(i)+Tau, sum(network_states(:,i))) + rand(1, sum(network_states(:,i)))*0.1*Tau;
        if ~isempty(spike_times) && (spike_times(1)<=0)
            spike_times(1) = 0;
        end
    elseif mode == 3
        U = time_points(i)+Tau;
        L = time_points(i);
        add = rand(1, sum(network_states(:,i)))*(U-L) + L;
    end
    spike_times = [spike_times add];
    spike_ids = [spike_ids find(network_states(:,i)')];
end
if mode == 3
    [spike_times,key] = sort(spike_times);
    spike_ids = spike_ids(key);
end
network_state = [network_states(:,end)'; ~network_states(:,end)'];

end