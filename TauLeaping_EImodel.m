function [spike_times,spike_ids,Irecons] = TauLeaping_EImodel(W,response_fn,beta,alpha,I,t_min,t_max,init_state,mode)
% This is an implementation of Tau Leaping Algorithm for the analysis
% of neuronal avalanches in Wilson-Cowan based stochastic
% population dynamics. Tau Leaping belongs to the family of
% stochastic simulation algoritmhs that are deriven from 
% Master Equation of an ODE. Tau Leaping, just like SSA,
% provides a compromise between continuous-deterministic solutions
% obtainable by direct solutions of ODEs and deterministic-stochastic
% solutions obtainable by exact/approximate solutions of the
% Master Equation. It is possible to interpret stochastic simulators
% also as an extension of systems with Langevin Dynamics. For more
% information, refer to the 2006 paper of Gillespie.
%
% Inputs:
% -W: is the weight matrix, n_neurons*n_neurons: W(i,j) is synaptic
% weight TO the ith neuron FROM the jth
% -response_fn: handle for response function (sigmoid, hyperbolic tan, etc.)
% -input is the net input, 1*n_neurons
% -alpha: is the rate at which active neurons decay to being
% quiescent, 1*n_neurons
% -beta: is the height of the response function, i.e. the rate at
% which saturated-input quiescent neurons become active, 1*n_neurons
% -init_state: is the initial state vector, 2*n_neurons
% -mode: is the operation mode of Tau Leaping, defines how spike times are
% spread around each leap interval.
%
% Outputs:
% -spike_times: timing of each neuronal spike according to the chosen mode
% -spike_ids: index of spiking neurons for each leap interval
% -Irecons: reconstruction of the temporal neuronal inputs from the input
% at each leap interval
%
% References:
% 
% Citation: Benayoun M, Cowan JD, van Drongelen W, Wallace E (2010) 
% Avalanches in a Stochastic Model of Spiking Neurons. 
% PLoS Comput Biol 6(7): e1000846. 
% https://doi.org/10.1371/journal.pcbi.1000846
% 
% Citation Yang Cao, Daniel T. Gillespie, Linda R. Petzold (2006); 
% Efficient step size selection for the tau-leaping simulation method. 
% J. Chem. Phys. 28 January 2006; 124 (4): 044109. 
% https://doi.org/10.1063/1.2159468

N = size(W,1);
Tau = 0.1;

time_points = t_min:Tau:t_max;
active = init_state(1,:)';
network_states = zeros(N, length(time_points));
event_count = zeros(N, length(time_points));

network_states(:,1) = active;
event_count(:,1) = zeros(N,1);
count = 1;

current = W*active + I(:,1);
propensity = beta .* ~active .* feval(response_fn, current) ...
    + alpha .* active;

for t=t_min:Tau:t_max
    num_events = poissrnd(propensity * Tau);
    
    for i = 1:N
        if num_events(i) > 0 
            if (num_events(i) // 2)
                if active(i) == 1
                    current += W(:,i);
                    active(i) = 0;
                else
                    current -= W(:,i);
                    active(i) = 1;
                end
            else
                continue;
            end
        else
            print('No event occurred for neuron %d between time %.2f and %.2f\n', i, t, t+Tau');
        end
    end
    count += 1;

    Ipast = I(:,floor(t+1));
    Icurr = I(:,floor(t+Tau+1));

    current = W*active + Icurr - Ipast;

    propensity = beta .* ~active .* feval(response_fn, current) ...
        + alpha .* active;
    network_states(:,count) = active;
    event_count(:,count) = num_events;
end

if mode == 0
    print('Accumulating spike times at the end of each time leap.\n')
    time_points = time_points
elseif mode == 1
    print('Uniformly distributing spike times throughout each time leap.\n')
    trecons = [];
    for i=1:size(event_count,2)
        trecons = [trecons linspace(time_points(i), time_points(i)+Tau, sum(event_count(:,i)))];
    end
    time_points = trecons
elseif mode == 2
    print('Uniformly distributing spike times with Gaussian Noise throughout each time leap.\n')
    trecons = [];
    for i=1:size(event_count,2)
        trecons = [trecons linspace(time_points(i), time_points(i)+Tau, sum(event_count(:,i))...
            + rand(1, sum(event_count(:,i))*0.1*Tau))]
        if trecons[1]<=0
            trecons[1] = 0;
        end
    end
    time_points = trecons
end

spike_times = [];
for i=1:len(time_points)
    spike_times = [spike_times repmat(time_points(i), 1, sum(event_count(:,i)))];
end
spike_ids = network_states;
end