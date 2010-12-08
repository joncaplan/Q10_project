% This model is based on Soto-Trevino et all 2005 and adds Q10's.

% TODO
%   Add Q10s
%   Save history
%   Display history
%   Add RK4 integration
%   Check signs on all currents - IN PROGRESS
%   Check unit prefixes in all calculations.

% Bugs
%   [Ca++] can go negative. 

function model()
    tic
    AB = 1; % Enumerate cell types.
    PD = 2;
    Soma = 1; % Enumerate compartments.
    AIS  = 2; % (Axon initial segment.)
    
    % Model parameters %
    dt                  = 0.05*10^-3;  % [s] Time step
    sim_length          = 0.0742;      % [s] Simulation length
    num_neurons         = 2;           % Number of neurons
    num_channels        = 10;          % Number of ion channel types.
    num_compartments    = 2;           % Number of compartments per neuron. Compartment 1 is soma. Compartment 2 is spike initiation zone.
%     g_gap               = 0.75*10^-6; % [S] Conductance between the soma of the two neurons.
%     neurons(AB).g_axial = 0.3 *10^-6; % [S] Conductance between the two compartments of the AB neuron.
%     neurons(PD).g_axial = 1.05*10^-6; % [S] Conductance between the two compartments of the PD neuron.

    % Testing with no connections between compartments.
    g_gap               = 0;     % [S] Conductance between the soma of the two neurons.
    neurons(AB).g_axial = 0;     % [S] Conductance between the two compartments of the AB neuron.
    neurons(PD).g_axial = 0;     % [S] Conductance between the two compartments of the PD neuron.    
    
    neurons(AB).compartments(AIS).C  =  1.5*10^-9; % [F] Capacitance of each compartment.
    neurons(AB).compartments(Soma).C =  9.0*10^-9; % [F]
    neurons(PD).compartments(AIS).C  =  6.0*10^-9; % [F]
    neurons(PD).compartments(Soma).C = 12.0*10^-9; % [F]
    
    % Constants for Ca++ dynamics
    R = 8.31447215;  % [J/mol/K] Ideal gas constant
    T = 273.15 + 18; % [Kelvin]  Temperature
    z = +2;          % Charge of Ca++.
    F = 96485.3399;  % [C/mol] Faraday's constant
    Ca_out             = 13000*10^-6; % [M]
    Ca_steady_state    = 0.5*10^-6;   % [M]
    neurons(AB).tau_Ca = 0.303;       % [s]
    neurons(PD).tau_Ca = 0.300;       % [s]
    neurons(AB).F_Ca   = 0.418*10^3;  % [M/A]
    neurons(PD).F_Ca   = 0.515*10^3;  % [M/A]

    % Reversal potentials.
    E_Na   =  50; % [mV] 
    E_K    = -80;
    E_Ca   =   0; % (Placeholder. Nernst calculation later.)
    E_Nap  =  50;
    E_H    = -20;
    E_KCa  = -80;
    E_A    = -80;
    E_proc =   0;
    E_all = [E_Na E_K E_Ca E_Ca E_Nap E_H E_K E_KCa E_A E_proc]*10^-3; % [V] Store reversal potentials in vector which matches the channel vector. (Ca values are placeholders.)

    % Maximal conductances
    %      Channels are:
    %      (AIS)  NA   K     (soma) CaT   CaS nap   H      K       KCa     A      proc   
    g_max_all = [ 300  52.5         55.2  9   2.7   0.054  1890    6000    200    570; ...       % AB % [S] Conductance values for all (non-leak) channels.
                 1100  150          22.5  60  4.38  0.219  1576.8  251.85  39.42  0   ]*10^-6;   % PD
             
    g_max_all = [0 0    0 0 0 0 0 0 0 0; 0 0    0 0 0 0 0 0 0 0 ];

    % Leak conductances
    %              Soma  AIS
	g_leaks   = [  0.045 0.0018;   ...     % AB % [S] Leak conductances [AB-soma, AB-AIS; PD-soma PD-AIS] (Note: Compartment order different than for g_max_all.)
                   0.105 0.00081 ]*10^-6;  % PD
               
    
    
    % Starting conditions
    starting_voltage = -70*10^-3; % [V] 
    for neuron = 1:num_neurons
        for compartment = 1:num_compartments 
            neurons(neuron).compartments(compartment).voltage = starting_voltage; 
        end
        neurons(neuron).Ca = Ca_steady_state; % Ca++ concentration in soma. (Axon lacks Ca++ gated or permiable channels.)
        for channel = 1:num_channels
            neurons(neuron).channels(channel).m = 0; % Begin with all channels closed.
            neurons(neuron).channels(channel).h = 0;
        end
    end
        

    % History variables
    V_hist = zeros(num_neurons, num_compartments, ceil(sim_length/dt));
    time_step = 1;
    
    %%%%%%%%%%%%%%%%%
    %   Main Loop   %
    %%%%%%%%%%%%%%%%%
    for sim_time = dt:dt:sim_length

        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find channel activation and inactivations  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for neuron = 1:num_neurons
            Ca = neurons(neuron).Ca;
            for channel = 1:num_channels
                if channel <= 2, compartment = 2; else compartment = 1; end % First two channels are in axon initial segment.
                V = neurons(neuron).compartments(compartment).voltage; % [V]
                [m h] = get_channel_state(neuron, channel, V, Ca, neurons(neuron).channels(channel).m, neurons(neuron).channels(channel).h, dt);
                neurons(neuron).channels(channel).m = m;
                neurons(neuron).channels(channel).h = h;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate channel currents  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for neuron = 1:num_neurons
            % Find Ca reversal potential via Nernst equation.
            Ca_in = neurons(neuron).Ca;
            E_Ca = R*T/(z*F)*log(Ca_out/Ca_in);
            Ca_in
            E_Ca
            E_all(3) = E_Ca; % Set reversal values for Ca++. 
            E_all(4) = E_Ca;
            
            for channel = 1:num_channels
                m = neurons(neuron).channels(channel).m;
                h = neurons(neuron).channels(channel).h;
                if (channel<=2), V = neurons(neuron).compartments(2).voltage; else V = neurons(neuron).compartments(1).voltage; end % Get voltage from correct compartment.
                g_max = g_max_all(neuron, channel); 
                switch channel % Lookup the channel activation and inactivation exponents.
                    case 1      % I_Na
                        a = 3; b = 1; 
                    case 2      % I_K (AIS)
                        a = 4; b = 0;
                    case 3      % I_CaT
                        a = 3; b = 1;
                    case 4      % I_CaS
                        a = 3; b = 0;
                    case 5      % I_nap
                        a = 3; b = 1;
                    case 6      % I_h
                        a = 1; b = 0;
                    case 7      % I_K (soma)
                        a = 4; b = 0;
                    case 8      % I_KCa
                        a = 4; b = 0;
                    case 9      % I_A
                        if (neuron == AB), a = 3; elseif (neuron == PD), a = 4; else STOP; end
                               b = 1;
                    case 10     % I_proc
                        a = 1; b = 0;
                end
                
                g_channel = g_max * m^a * h^b;
                if (g_channel < 0)
                    disp ('Negative conductance. Stopping.')
                    g_max
                    m
                    a
                    h
                    b
                    STOP
                end
                
                neurons(neuron).channels(channel).I = g_channel * (V - E_all(channel));
                if (neurons(neuron).channels(channel).I ~= abs(neurons(neuron).channels(channel).I) && neurons(neuron).channels(channel).I ~= -abs(neurons(neuron).channels(channel).I))
                    disp ('Current has complex component or is NaN. Stopping.')
                    neurons(neuron).channels(channel).I
                    neuron
                    channel
                    V
                    g_channel
                    E_all(channel)
                    STOP
                end
                if (g_channel ~= abs(g_channel) && -g_channel ~= abs(g_channel)) % Testing for imaginary component.
                    neuron
                    channel
                    V
                    g_channel
                    E_all(channel)
                    STOP
                end
            end
            for compartment = 1:2 % Calculate leak currents
                I_leak(neuron, compartment) = -g_leaks(neuron, compartment) * neurons(neuron).compartments(compartment).voltage;
                leak = I_leak(neuron, compartment);
                if (leak == Inf || leak == -Inf)
                    disp('+/- Inf leak current. Stopping.')
                    leak
                    STOP
                end
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find new voltage of each compartment %      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        % Find coupling current
        dV_neurons = neurons(PD).compartments(Soma).voltage - neurons(AB).compartments(Soma).voltage;
        coupling_current =  dV_neurons * g_gap; % NOTE: *Positive* current if PD is at higher potential. (Current flows into cell AB.)
        
        % Find currents within each neuron
        for neuron = 1:num_neurons
            % Sum the currents for each compartment
            dV_axial = neurons(neuron).compartments(Soma).voltage - neurons(AIS).compartments(Soma).voltage; % NOTE: Positive for current flowing into AIS from soma.
            axial_current   = dV_axial * neurons(neuron).g_axial;
            AIS_current     = 0;
            somatic_current = 0;
            for channel = 1:2
                AIS_current     = AIS_current     + neurons(neuron).channels(channel).I
            end
            AIS_current = AIS_current + I_leak(neuron, AIS)

            for channel = 3:num_channels
                %channel
                %neurons(neuron).channels(channel).I
                somatic_current = somatic_current + neurons(neuron).channels(channel).I;
            end
            somatic_current = somatic_current + I_leak(neuron, Soma)
                        
            coupling_current_direction = ((neuron==1)*2 -1); % +1 for cell 1. -1 for cell 2.
            total_soma_current = somatic_current - axial_current + coupling_current*coupling_current_direction;  % FIX: Verify signs on currents. 
            total_AIS_current  = AIS_current     + axial_current;
            
            dV_soma = total_soma_current / neurons(neuron).compartments(Soma).C * dt;
            dV_axon = total_AIS_current  / neurons(neuron).compartments(AIS ).C * dt;
            
            neurons(neuron).compartments(Soma).voltage = neurons(neuron).compartments(Soma).voltage + dV_soma;
            neurons(neuron).compartments(AIS ).voltage = neurons(neuron).compartments(AIS ).voltage + dV_axon;
            
            V_hist(neuron, Soma, time_step) = neurons(neuron).compartments(Soma).voltage;
            V_hist(neuron, AIS,  time_step) = neurons(neuron).compartments(AIS ).voltage;
            
            soma_voltage = neurons(neuron).compartments(Soma).voltage
            AIS_voltage  = neurons(neuron).compartments(AIS ).voltage
            
            if (soma_voltage == Inf || AIS_voltage == Inf)
                disp ('Infinite voltage. Stopping.')
                STOP
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%
        % Calculate [Ca++]  %
        %%%%%%%%%%%%%%%%%%%%%

        for neuron = 1:num_neurons
            I_CaT = neurons(neuron).channels(3).I;
            I_CaS = neurons(neuron).channels(4).I;            
            I_Ca  = I_CaT + I_CaS; 
            Ca = neurons(neuron).Ca;
            F_Ca = neurons(neuron).F_Ca;
            tau_Ca = neurons(neuron).tau_Ca;
            integration = 'ExpEuler';
            switch integration
                
                case 'RK4'
                    k1 = ( -F_Ca * I_Ca -  Ca              + Ca_steady_state) / tau_Ca;
                    k2 = ( -F_Ca * I_Ca - (Ca + 0.5*dt*k1) + Ca_steady_state) / tau_Ca;
                    k3 = ( -F_Ca * I_Ca - (Ca + 0.5*dt*k2) + Ca_steady_state) / tau_Ca;
                    k4 = ( -F_Ca * I_Ca - (Ca +     dt*k3) + Ca_steady_state) / tau_Ca;
                    dCa = (1/6)*dt*(k1 + 2*k2 + 2*k3 + k4);
                    neurons(neuron).Ca = Ca + dCa;

                case 'ForwardEuler'
                    dCa = (-F_Ca * I_Ca - Ca + Ca_steady_state) / tau_Ca * dt;
                    neurons(neuron).Ca = neurons(neuron).Ca + dCa;
                
                case 'ExpEuler'
                    neurons(neuron).Ca = Ca_steady_state + (Ca - Ca_steady_state)*exp(-dt/tau_Ca);  %inf + (curr - inf)*exp(-dt/tau); 
                
            end
            if (neurons(neuron).Ca < 0)
                disp('Negative [Ca++] value. Stopping')
                Ca 
                dCa
                I_Ca

                F_Ca
                tau_Ca
                disp('.')
                STOP
            end
        end
        if sim_time/100 == floor(sim_time/100),
            sim_time
        end
        time_step = time_step+1
    end % End main loop
    toc
    
    
    figure
    hold on
        plot(squeeze(V_hist(AB,Soma,:)),'r')
        plot(squeeze(V_hist(AB,AIS, :)),'g')
        plot(squeeze(V_hist(PD,Soma,:)),'b')
        plot(squeeze(V_hist(PD,AIS, :)),'c')
        legend('AB Soma','AB AIS','PD Soma','PD AIS')
    hold off
    
end % End function


%%%%%%%%%%%%%%
% Functions  %
%%%%%%%%%%%%%%

function [m, h] = get_channel_state(neuron, channel, V, Ca, old_m, old_h, dt)
    if neuron == 1, neuron_type = 1; else neuron_type = 2; end

    switch channel
        case 1 % I_Na (Axon initial segment)
            m_inf = sigmoid(-1, 24.7, 5.29, V);
            h_inf = sigmoid(+1, 48.9, 5.18, V);
            tau_m = tau_sigmoid(1.32, -1.26, -1,  120, 25, V);
            tau_h = tau_sigmoid(0,     0.67, -1, 62.9, 10, V)  * tau_sigmoid(1.5, 1, +1, 34.9, 3.6, V);

        case 3 % I_CaT
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
            
        case 4 % I_CaS
            m_inf = sigmoid(-1, 22.5, 8.5, V);
            h_inf = 1;
            tau_m = tau_sigmoid(16, -13.1, -1, 25.1, 26.4, V);
            tau_h = 1;
            
        case 5 % I_NaP
            m_inf = sigmoid(-1, 26.8, 8.2, V);
            h_inf = sigmoid( 1, 48.5, 4.8, V);
            tau_m = tau_sigmoid(19.8, -10.7, -1, 26.5,  8.6,  V);
            tau_h = tau_sigmoid(666,  -379,  -1, 33.6,  11.7, V);
            
        case 6 % I_h
            m_inf = sigmoid( 1, 70, 6, V);
            h_inf = 1;
            tau_m = tau_sigmoid(272, +1499, -1, 42.2, 8.73, V);
            tau_h = 1;
            
        case {2, 7} % I_K {axon initial segment, soma}
            m_inf = sigmoid(-1, 14.2, 11.8, V);
            h_inf = 1;
            tau_m = tau_sigmoid( 7.2, -6.4, -1, 28.3, 19.2, V);
            tau_h = 1;
            
        case 8 % I_KCa
            switch neuron_type
                case 1
                    m_inf = Ca/(Ca + 30) * sigmoid(-1, 51, 4, V);
                case 2
                    m_inf = Ca/(Ca + 30) * sigmoid(-1, 51, 8, V);
                otherwise 
                    STOP
            end        
            h_inf = 1;
            tau_m = tau_sigmoid(90.3, -75.09, -1, 46, 22.7, V);
            tau_h = 1;
            
        case 9 % I_A 
            m_inf = sigmoid( -1, 27,   8.7,  V);
            h_inf = sigmoid( +1, 56.9, 4.9,  V);
            tau_m = tau_sigmoid( 11.6, -10.4, -1, 32.9, 15.2,  V);
            tau_h = tau_sigmoid( 38.6, -29.2, -1, 38.9, 26.5,  V);    
            
        case 10 % I_proc
            m_inf = sigmoid( -1, 12, 3.05,   V);
            h_inf = 1;
            tau_m = 0.5;
            tau_h = 1;
            
        otherwise
            STOP
    end
    % Integrate to find m and h at next time step.
    m = m_inf + (old_m - m_inf)*exp(-dt/tau_m);       %inf + (curr - inf)*exp(-dt/tau); 
    h = h_inf + (old_h - h_inf)*exp(-dt/tau_h);       %inf + (curr - inf)*exp(-dt/tau); 
    if (false)
        dm = (m_inf - old_m)/tau_m * dt;  % FIX: Forward Euler to expEuler.
        dh = (h_inf - old_h)/tau_h * dt;
        m = old_m + dm;
        h = old_h + dh;
    end
    if (m < 0 || m > 1 ), disp ('m out of range. Stopping.')
        m
        old_m
        m_inf
        tau_m
        STOP
    end
    if (h < 0 || h > 1 ), disp ('h out of range. Stopping.')
        h
        STOP
    end
    
end

% Sigmoid functions for calculating channel dynamics
function [answer] = sigmoid(a,b,c, V)
    answer = 1/(1+exp(a*(V+b)/c));
end

function[tau] = tau_sigmoid(a, b, c, d, e, V)
        tau = a + b/(1+exp(c*(V+d)/e));
end