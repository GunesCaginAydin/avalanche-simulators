function [spike_times,spike_ids,network_state,Irc,times] = Gillespie_EImodel(W,response_fn,beta,alpha,I,t_min,t_max,init_state,varargin)

% Simulates a 2-state Wilson-Cowan model with the Gillespie algorithm and
% arbitrary weight matrix W, inputs I, and transition rates (alpha and 
% tranfert function), between times t_min and t_max.
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
%--------------------------------------------------------------------------

N = size(W,1); 

factor=10;
expected_events=N*(t_max-t_min)*factor;
times = zeros(1, expected_events);
updates = zeros(1, expected_events);
new_states = zeros(1, expected_events);
Irc = [];

event_no = 0;
curr_time = t_min;

active = init_state(1,:)';

currents = W*active + I(:,1);
trans = beta .* ~active .* feval(response_fn,currents) + ...
    alpha .* active;

while (curr_time < t_max)
    cum_trans=cumsum(trans);
    total_trans = cum_trans(end);
    if(cum_trans(N)==0)
        fprintf('Total propensity is effectively 0 - terminating simulation')
        return
    end

    dt = -log(rand)/total_trans;
    test_variable = total_trans*rand;
    i_update = find(cum_trans >= test_variable,1,'first');
    active(i_update) = double(~active(i_update));
    
    if active(i_update) == 1
        currents = currents + W(:,i_update);
    else
        currents = currents - W(:,i_update);
    end
    Ipast = I(:,floor(curr_time+1));
    curr_time = curr_time + dt;
    Icurr = I(:,floor(curr_time+1));
    
    currents = currents +  (Icurr - Ipast);
    trans = beta .* ~active .*feval(response_fn,currents) + ...
        alpha.*active;

    event_no = event_no+1;
    times(event_no) = curr_time;
    updates(event_no) = i_update;
    new_states(:, event_no) = active(i_update);
    Irc = [Irc Icurr];
end

event_no = event_no - 1;
times=times(1:event_no);
updates=updates(1:event_no);
new_states=new_states(1:event_no);
network_state = [active'; ~active'];

spikes = find(new_states);

spike_times = times(spikes);
spike_ids = updates(spikes);

return