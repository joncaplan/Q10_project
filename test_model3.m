% This model is based on Soto-Trevino et all 2005 and adds Q10's.
% This branch created to test effect of different temperatures and
% timesteps on the soma spiking. (12/30/2010)

% TODO
%   Add Q10s
%   Add RK4 integration


function test_model3()
    tic
    AB = 1; % Enumerate cell types.
    PD = 2;
    Soma = 1; % Enumerate compartments.
    AIS  = 2; % (Axon initial segment.)
    
    % Model parameters %
    dt                  = 0.05*10^-3;  % [s] Time step
%    I_inj               = [0 0];       % [Amp] Current injection [AB_soma PD_soma]
    sim_length          = 2;         % [s] Simulation length
    num_neurons         = 2;           % Number of neurons
    num_channels        = 10;          % Number of ion channel types.
    num_compartments    = 2;           % Number of compartments per neuron. Compartment 1 is soma. Compartment 2 is spike initiation zone.
%     g_gap               = 0.75*10^-6; % [S] Conductance between the soma of the two neurons.
%     neurons(AB).g_axial = 0.3 *10^-6; % [S] Conductance between the two compartments of the AB neuron.
%     neurons(PD).g_axial = 1.05*10^-6; % [S] Conductance between the two compartments of the PD neuron.

    % Testing with no connections between compartments.
%     g_gap               = 0;     % [S] Conductance between the soma of the two neurons.
     neurons(AB).g_axial = 0;     % [S] Conductance between the two compartments of the AB neuron.
     neurons(PD).g_axial = 0;     % [S] Conductance between the two compartments of the PD neuron.    
    
    neurons(AB).compartments(AIS).C  =  1.5*10^-9; % [F] Capacitance of each compartment.
    neurons(AB).compartments(Soma).C =  9.0*10^-9; % [F]
    neurons(PD).compartments(AIS).C  =  6.0*10^-9; % [F]
    neurons(PD).compartments(Soma).C = 12.0*10^-9; % [F]
    
    % Constants for Ca++ dynamics
    R = 8.31447215;  % [J/mol/K] Ideal gas constant
    T = 273.15 + 10; % [Kelvin]  Temperature (Soto-trevino's animals are tested at 18C. Buchholtz et al 1992 use 10C.)
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

             
%    g_max_all = [0 0         0 0 0 0 0 0 0 0; 0 0          0 0 0 0 0 0 0 0 ];
    %g_max_all = [300 52.5    0 0 0 0 0 0 0 0; 1100 150     0 0 0 0 0 0 0 0 ]*10^-6; % Just turn on AIS conductances.
    
    % Just turn on soma conductances.
%     g_max_all = [ 0  0          55.2  9   2.7   0.054  1890    6000    200    570; ...       % AB % [S] Conductance values for all (non-leak) channels.
%                   0  0          22.5  60  4.38  0.219  1576.8  251.85  39.42  0   ]*10^-6;   % PD    

         % Channel  Na  K(AIS)  CaT  CaS  nap  h   K(soma) KCa  A   proc
    a_exponents = [ 3   4       3    3    3    1   4       4    3   1     ; ... % AB
                    3   4       3    3    3    1   4       4    4   1    ];    % PD
    b_exponents = [ 1   0       1    0    1    0   0       0    1   0    ];    % Both AB and PD    
    
    % Leak conductances
    %              Soma  AIS
	g_leaks   = [  0.045 0.0018;   ...     % AB % [S] Leak conductances [AB-soma, AB-AIS; PD-soma PD-AIS] (Note: Compartment order different than for g_max_all.)
                   0.105 0.00081 ]*10^-6;  % PD
               
    % Leak reversal potentials
    %              Soma  AIS
    E_leak     = [ -50   -60;  ...         % AB [V]
                   -55   -55  ] * 10^-3;   % PD
    
    % Starting conditions
    starting_voltage = +20*10^-3; % [V] 
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
    %I_hist = zeros(num_neurons, num_channels+1,     ceil(sim_length/dt)); % +1 is for leak channel.
    global Na K_AIS CaT CaS nap H K_soma KCa A proc Leak_AIS Leak_soma
    Na=1; K_AIS=2; CaT=3; CaS=4; nap=5; H=6; K_soma=7; KCa=8; A=9; proc=10;Leak_AIS=11;Leak_soma=12; % Enumeration of all the currents.
    time_step = 1;
    
    %%%%%%%%%%%%%%%%%
    %   Main Loop   %
    %%%%%%%%%%%%%%%%%
    tic
    for sim_time = dt:dt:sim_length

        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find channel activation and inactivations  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for neuron = 1:num_neurons
            Ca = neurons(neuron).Ca;
            for channel = 1:num_channels
                if channel <= 2, compartment = AIS; else compartment = Soma; end % First two channels are in axon initial segment.
                V = neurons(neuron).compartments(compartment).voltage; % [V]
                [m h] = get_channel_state(neuron, channel, V, Ca, neurons(neuron).channels(channel).m, neurons(neuron).channels(channel).h, dt);
                neurons(neuron).channels(channel).m = m;
                neurons(neuron).channels(channel).h = h;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate channel conductances  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Find channel conductances
        for neuron = 1:num_neurons
            % Find Ca reversal potential via Nernst equation.
            Ca_in = neurons(neuron).Ca;
            E_Ca = R*T/(z*F)*log(Ca_out/Ca_in);
            E_all(3) = E_Ca; % Set reversal values for Ca++. 
            E_all(4) = E_Ca;
            
            for channel = 1:num_channels
                m = neurons(neuron).channels(channel).m;
                h = neurons(neuron).channels(channel).h;
                g_max = g_max_all(neuron, channel); 
                
                a = a_exponents(neuron, channel);
                b = b_exponents(channel);
                
                g_channels(neuron, channel) = g_max * m^a * h^b;
                g_channel=g_channels(neuron, channel);
                if ((g_channel ~= abs(g_channel) && -g_channel ~= abs(g_channel)) || (g_channel < 0) ) % Testing for imaginary or negative component.
                    disp ('Imaginary or negative conductance. Stopping.'); STOP
                end                
            end
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find new voltage of each compartment %      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        % Find coupling current
%         dV_neurons = neurons(PD).compartments(Soma).voltage - neurons(AB).compartments(Soma).voltage;
%         coupling_current =  dV_neurons * g_gap; % NOTE: *Positive* current if PD is at higher potential. (Current flows into cell AB.)
        
        % Find currents within each neuron
        for neuron = 1:num_neurons
            % Sum the currents for each compartment
%             dV_axial = neurons(neuron).compartments(Soma).voltage - neurons(neuron).compartments(AIS).voltage; % NOTE: Positive for current flowing into AIS from soma.
%             axial_current   = dV_axial * neurons(neuron).g_axial;
%             if axial_current ~= 0,
%                 disp ('Non zero axial current')
%                 axial_current
%                 disp ('Stopping.')
%                 STOP; 
%             end %%%%% TESTING uncoupled compartments should have no axial current.
%             coupling_current_direction = ((neuron==1)*2 -1); % +1 for cell 1. -1 for cell 2.
    
            I_weighted_sum    = zeros(1, num_compartments); % Sum of single steady state channel currents gives weighted sum of their contributions at steady state.
            current_sum       = zeros(1, num_compartments); % Sum of all currents     (active, leak, axial, gap and applied).
            conductance_sum   = zeros(1, num_compartments); % Sum of all conductances (active, leak, axial and gap)
            
            for compartment = 1:num_compartments
                V = neurons(neuron).compartments(compartment).voltage;
                if compartment == AIS, channels = 1:2; else channels = 3:num_channels; end
                for channel = channels % Iterate through the set of channels.

                    channel_current = g_channels(neuron, channel)* (E_all(channel) -V);
%                     if (channel ==3 || channel == 4) % Record Ca++ channel currents (CaT and CaS) for later use in finding new [Ca++].
%                         neurons(neuron).channels(channel).I = channel_current;
%                     end
                    I_hist(neuron, channel, time_step) =  channel_current;
                    current_sum(compartment) = current_sum(compartment) + channel_current; 
                    I_weight     = g_channels(neuron, channel)*  E_all(channel); % Single channel current weight.
                    I_weighted_sum(compartment)  = I_weighted_sum(compartment)  + I_weight;
                    conductance_sum(compartment) = conductance_sum(compartment) + g_channels(neuron, channel);
                end

                I_leak_weight  =  g_leaks(neuron, compartment)*(E_leak(neuron, compartment) ) ;
                leak_current   =  g_leaks(neuron, compartment)*(E_leak(neuron, compartment) -V) ;
                I_hist(neuron, Leak_soma, time_step) = leak_current;
                if compartment == Soma
                    conductance_sum(Soma) = conductance_sum(Soma) + g_leaks(neuron, Soma); % TESTING %%%%%%%%%%%%+ neurons(neuron).g_axial + g_gap;
                    I_weighted_sum(Soma)  = I_weighted_sum(Soma)  + I_leak_weight; % TESTING %%%%%%%%%%%% + axial_current + coupling_current*coupling_current_direction + I_inj(compartment);
%                     current_sum(Soma)     = current_sum(Soma)     + leak_current  + axial_current + coupling_current*coupling_current_direction + I_inj(compartment);
                else
                    if compartment ~= AIS, STOP; end
                    conductance_sum(AIS)  = conductance_sum(AIS)  + g_leaks(neuron, AIS); %%%%%%%%% TESTING + neurons(neuron).g_axial;
                    I_weighted_sum(AIS)   = I_weighted_sum(AIS)   + I_leak_weight; % TESTING %%%%%%%%%%%% + axial_current;
%                     current_sum(AIS)      = current_sum(AIS)      + leak_current  + axial_current;                
                end
                
            end
 
            % Now calculate new voltage from currents.
            for compartment = 1:num_compartments
                V = neurons(neuron).compartments(compartment).voltage;
                V_inf  =  I_weighted_sum(compartment)/conductance_sum(compartment);   % Weighted sum of reversal potentials.
                integration_method = 'exp_Euler';
                switch integration_method
                    case 'exp_Euler'
                        tau_V  =  neurons(neuron).compartments(compartment).C/(conductance_sum(compartment));   % Membrane time constant.
                        neurons(neuron).compartments(compartment).voltage      =  V_inf + (V-V_inf)*exp(-dt/tau_V);
%                     case 'forward_Euler'
%                         dV = current_sum(compartment)/neurons(neuron).compartments(compartment).C*dt;
%                         neurons(neuron).compartments(compartment).voltage = V + dV;
                    otherwise
                        STOP
                end
            end

            % Record history.
            V_hist(neuron, Soma, time_step) = neurons(neuron).compartments(Soma).voltage;
            V_hist(neuron, AIS,  time_step) = neurons(neuron).compartments(AIS ).voltage;
        end
        
        %%%%%%%%%%%%%%%%%%%%%
        % Calculate [Ca++]  %
        %%%%%%%%%%%%%%%%%%%%%
        
        for neuron = 1:num_neurons       
            I_CaT = I_hist(neuron, 3, time_step); % WAS: neurons(neuron).channels(3).I;
            I_CaS = I_hist(neuron, 4, time_step); % WAS: neurons(neuron).channels(4).I; 
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
                otherwise
                    STOP
            end
            if (neurons(neuron).Ca < 0)
                disp('Negative [Ca++] value. Stopping'); STOP
            end
        end
        
        %%% Update and display time step.
        if sim_time/100 == floor(sim_time/100),
            disp(num2str(sim_time))
        end
        time_step = time_step+1;
        if (time_step/100 == floor(time_step/100))
            disp(num2str(time_step))
        end
    end % End main loop
    toc
    
    %%%%%%%%%%%%%%%%%
    % Draw Figures  %
    %%%%%%%%%%%%%%%%%
    
    figure
    hold on
        plot( (dt:dt:sim_length)*1000, squeeze(V_hist(AB,Soma,:))*1000, 'r')
        plot( (dt:dt:sim_length)*1000, squeeze(V_hist(AB,AIS, :))*1000, 'g')
        plot( (dt:dt:sim_length)*1000, squeeze(V_hist(PD,Soma,:))*1000, 'b')
        plot( (dt:dt:sim_length)*1000, squeeze(V_hist(PD,AIS, :))*1000, 'c')
        title(strcat('Membrane potential (T=', num2str(T-273.15),'C dt=',num2str(dt*1000),' ms)' ))
        legend('AB Soma','AB AIS','PD Soma','PD AIS')
        xlabel('Time (ms)')
        ylabel('Membrane potential (mV)')
    hold off
    
    for neuron = 1:num_neurons
    figure
    hold on
        plot( (dt:dt:sim_length)*1000, squeeze(I_hist(neuron,       CaT,:))*1000, 'r')
        plot( (dt:dt:sim_length)*1000, squeeze(I_hist(neuron,       CaS,:))*1000, 'g')
        plot( (dt:dt:sim_length)*1000, squeeze(I_hist(neuron,       nap,:))*1000, 'b')
        plot( (dt:dt:sim_length)*1000, squeeze(I_hist(neuron,         H,:))*1000, 'c')
        plot( (dt:dt:sim_length)*1000, squeeze(I_hist(neuron,    K_soma,:))*1000, 'y')
        plot( (dt:dt:sim_length)*1000, squeeze(I_hist(neuron,       KCa,:))*1000, 'm')
        plot( (dt:dt:sim_length)*1000, squeeze(I_hist(neuron,         A,:))*1000, 'k')
        if neuron==AB, plot( (dt:dt:sim_length)*1000, squeeze(I_hist(neuron,      proc,:))*1000, 'r.'), end % No proctolin current in PD.
        plot( (dt:dt:sim_length)*1000, squeeze(I_hist(neuron, Leak_soma,:))*1000, 'g.')
        if neuron== AB
            legend('Ca_T','Ca_S','nap','H','K soma','K_C_a','A','proctolin', 'Leak soma')
        else
            legend('Ca_T','Ca_S','nap','H','K soma','K_C_a','A',             'Leak soma')
        end
        
        if neuron==1
            title(strcat('AB soma currents (T=', num2str(T-273.15),'C dt=',num2str(dt*1000),' ms)' ))
        else
            title(strcat('PD soma currents (T=', num2str(T-273.15),'C dt=',num2str(dt*1000),' ms)' ))
        end
    hold off
    end
    
end % End function


%%%%%%%%%%%%%%
% Functions  %
%%%%%%%%%%%%%%

function [m, h] = get_channel_state(neuron, channel, V, Ca, old_m, old_h, dt)
    if neuron == 1, neuron_type = 1; else neuron_type = 2; end
    
    V  = V  * 1000; % Convert voltage from V to mV for use in functions table. (Elsewhere voltage is always in Volts.)
    Ca = Ca * 10^6; % Convert intracellular [Ca++] from Mol to μM. (Elsewhere [Ca++] is always in Mol.)

    % Channel dynamics functions table.
    Na=1; K_AIS=2; CaT=3; CaS=4; nap=5; H=6; K_soma=7; KCa=8; A=9; proc=10; % Enumeration of all the currents.

    switch channel
        case Na % I_Na (Axon initial segment)
            m_inf = sigmoid(-1, 24.7, 5.29, V);
            h_inf = sigmoid(+1, 48.9, 5.18, V);
            tau_m = tau_sigmoid(1.32, -1.26, -1,  120, 25, V);
            tau_h = tau_sigmoid(0,     0.67, -1, 62.9, 10, V)  * tau_sigmoid(1.5, 1, +1, 34.9, 3.6, V);

        case CaT % I_CaT
            m_inf = sigmoid(-1, 25, 7.2, V);
            h_inf = sigmoid(+1, 36, 7  , V);
            tau_m = tau_sigmoid(55, -49.5, -1, 58, 17, V);
            switch neuron_type
                case 1 % AB
                    tau_h = tau_sigmoid(87.5,  -75, -1, 50, 16.9, V);
                case 2 % PD
                    tau_h = tau_sigmoid(350,  -300, -1, 50, 16.9, V);
                otherwise 
                    STOP
            end
            
        case CaS % I_CaS
            m_inf = sigmoid(-1, 22.5, 8.5, V);
            tau_m = tau_sigmoid(16, -13.1, -1, 25.1, 26.4, V);
            
        case nap % I_NaP
            m_inf = sigmoid(-1, 26.8, 8.2, V);
            h_inf = sigmoid( 1, 48.5, 4.8, V);
            tau_m = tau_sigmoid(19.8, -10.7, -1, 26.5,  8.6,  V);
            tau_h = tau_sigmoid(666,  -379,  -1, 33.6,  11.7, V);
            
        case H % I_h
            m_inf = sigmoid( 1, 70, 6, V);
            tau_m = tau_sigmoid(272, +1499, -1, 42.2, 8.73, V);
            
        case {K_AIS, K_soma} % I_K {axon initial segment, soma}
            m_inf = sigmoid(-1, 14.2, 11.8, V);
            tau_m = tau_sigmoid( 7.2, -6.4, -1, 28.3, 19.2, V);
            
        case KCa % I_KCa
            switch neuron_type
                case 1
                    m_inf = Ca/(Ca + 30) * sigmoid(-1, 51, 4, V);
                case 2
                    m_inf = Ca/(Ca + 30) * sigmoid(-1, 51, 8, V);
                otherwise 
                    STOP
            end        
            tau_m = tau_sigmoid(90.3, -75.09, -1, 46, 22.7, V);
            
        case A % I_A 
            m_inf = sigmoid( -1, 27,   8.7,  V);
            h_inf = sigmoid( +1, 56.9, 4.9,  V);
            tau_m = tau_sigmoid( 11.6, -10.4, -1, 32.9, 15.2,  V);
            tau_h = tau_sigmoid( 38.6, -29.2, -1, 38.9, 26.5,  V);    
            
        case proc % I_proc
            m_inf = sigmoid( -1, 12, 3.05,   V);
            tau_m = 0.5*10^-3;
            
        otherwise
            STOP
                        
    end
    
    
    % Integrate to find m and h at next time step.
    m = m_inf + (old_m - m_inf)*exp(-dt/tau_m);       %inf + (curr - inf)*exp(-dt/tau); % Exp Euler.
    switch channel
        case {K_AIS, CaS, H, K_soma, KCa, proc} % Channels which do not inactivate.
            h = 1;
        case {Na, CaT, nap, A} % Channels which inactivate.
            h = h_inf + (old_h - h_inf)*exp(-dt/tau_h);       %inf + (curr - inf)*exp(-dt/tau); 
        otherwise
            STOP
    end
    if (false)
        dm = (m_inf - old_m)/tau_m * dt;  % OLD CODE: Forward Euler.
        dh = (h_inf - old_h)/tau_h * dt;
        m = old_m + dm;
        h = old_h + dh;
    end
    if (m < 0 || m > 1 ), disp ('m out of range. Stopping.'); STOP; end
    if (h < 0 || h > 1 ), disp ('h out of range. Stopping.'); STOP; end
    
end

% Sigmoid functions for calculating channel dynamics
function [answer] = sigmoid(a,b,c, V)
    answer = 1/(1+exp(a*(V+b)/c));
end

function[tau] = tau_sigmoid(a, b, c, d, e, V)
        tau = (a + b/(1+exp(c*(V+d)/e)))*10^-3;
end