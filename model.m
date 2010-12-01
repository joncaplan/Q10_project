% This model is based on Soto-Trevino et all 1995 and adds Q10's.

function model()
    % Model parameters %
    dt = 0.05;              % [ms] Time step
    sim_length = 10 * 1000; % [ms] Simulation length
    neuron_num      = 2;    % Number of neurons
    channels        = 9;    % Number of ion channel types.
    compartment_num = 2;    % Number of compartments per neuron. % Compartment 1 is soma. Compartment 2 is spike initiation zone.
    Ca_steady_state = 1;

    % Set up model starting conditions
    starting_voltage = -70; % [mV] 
    for neuron = 1:neuron_num
        for compartment = 1:compartment_num 
            neurons(neuron).compartments(compartment).voltage = starting_voltage;
        end
    end
    
    for neuron = 1:neuron_num
        neurons(neuron).Ca = Ca_steady_state;
    end

    %%%%%%%%%%%%%%%%%
    %   Main Loop   %
    %%%%%%%%%%%%%%%%%
    for sim_time = dt:dt:sim_length
        for neuron = 1:neuron_num
            for compartment = 1:compartment_num

                % Find channel activation and inactivations
                for channel = 1:channels
                    channel_state(channel) = get_channel_state(neuron, channel);
                end

            end
        end    
    end

%%%%%%%%%%%%%%
% Functions  %
%%%%%%%%%%%%%%

function [m, h] = get_channel_state(neuron, channel, compartment)
    V = neurons(neuron).voltage(compartment); % [mV]
    if neuron == 1, neuron_type = 1; else neuron_type = 2; end

    switch channel
        case 1 % I_Na
            m_inf = sigmoid(-1, 24.7, 5.29, V);
            h_inf = sigmoid(+1, 48.9, 5.18, V);
            tau_m = tau_sigmoid(1.32, -1.26, -1,  120, 25, V);
            tau_h = tau_sigmoid(0,     0.67, -1, 62.9, 10, V)  * tau_sigmoid(1.5, 1, +1, 34.9, 3.6, V);

        case 2 % I_CaT
            m_inf = sigmoid(-1, 25, 7.2, V);
            h_inf = sigmoid(+1, 36, 7  , V);
            tau_m = tau_sigmoid(55, -49.5, -1, 58, 17, V);
            switch neuron_type
                case 1
                    tau_h = tau_sigmoid(87.5,  -75, -1, 50, 16.9, V);
                case 2
                    tau_h = tau_sigmoid(350,  -300, -1, 50, 16.9, V);
                otherwise 
                    STOP
            end
            
        case 3 % I_CaS
            m_inf = sigmoid(-1, 22.5, 8.5, V);
            h_inf = 1;
            tau_m = tau_sigmoid(16, -13.1, -1, 25.1, 26.4, V);
            tau_h = 1;
            
        case 4 % I_NaP
            m_inf = sigmoid(-1, 26.8, 8.2, V);
            h_inf = sigmoid( 1, 48.5, 4.8, V);
            tau_m = tau_sigmoid(19.8, -10.7, -1, 26.5,  8.6,  V);
            tau_h = tau_sigmoid(666,  -379,  -1, 33.6,  11.7, V);
            
        case 5 % I_h
            m_inf = sigmoid( 1, 70, 6, V);
            h_inf = 1;
            tau_m = tau_sigmoid(272, +1499, -1, 42.2, 8.73, V);
            tau_h = 1;
            
        case 6 % I_K
            m_inf = sigmoid(-1, 14.2, 11.8, V);
            h_inf = 1;
            tau_m = tau_sigmoid( 7.2, -6.4, -1, 28.3, 19.2, V);
            tau_h = 1;
            
        case 7 % I_KCa
            Ca = neurons(neuron).Ca;
            switch neuron_type
                case 1
                    m_inf = Ca/(Ca + 30) * sigmoid(-1, 51, 4);
                case 2
                    m_inf = Ca/(Ca + 30) * sigmoid(-1, 51, 8);
                otherwise 
                    STOP
            end        
            h_inf = 1;
            tau_m = tau_sigmoid(90.3, -75.09, -1, 46, 22.7, V);
            tau_h = 1;
            
        case 8 % I_A 
            m_inf = sigmoid( -1, 27,   8.7,  V);
            h_inf = sigmoid( +1, 56.9, 4.9,  V);
            tau_m = tau_sigmoid( 11.6, -10.4, -1, 32.9, 15.2,  V);
            tau_h = tau_sigmoid( 38.6, -29.2, -1, 38.9, 26.5,  V);    
            
        case 9 % I_proc
            m_inf = sigmoid( -1, 12, 3.05,   V);
            h_inf = 1;
            tau_m = 0.5;
            tau_h = 1;
            
        otherwise
            STOP
    end
    % Integrate to find m and h at next time step. 
    m = m_inf * tau_m; %%%%%%%%% Placeholder code %%%%%%%%%%
    h = h_inf * tau_h;

% Sigmoid functions for calculating channel dynamics
function [answer] = sigmoid(a,b,c, V)
    answer = 1/(1+exp(a*(V+b)/c));

function[tau] = tau_sigmoid(a, b, c, d, e, V)
        tau = a + b(1+exp(c*(V+d)/e));
