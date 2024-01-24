classdef FCGLE_simulator
    %FCGLE_simulator Solver of the forced complex Ginzburg-Landau equation
    % 15/12/2022: Corrected G = -cD; Delta = cNL;
    % (Note: Nature --> phys. convention; PRL --> eng. convention)
    % Noise (normally not uniformly distributed) on BOTH real and imaginary parts (important for homoclones)
    % Removed feedback term (et) appearing for long iterations
    
    % Marco Piccardo - Harvard University
    % Script to solve the FCGLE for an active cavity adapted from a LLE solver
    % of a passive cavity. The original LLE solver was based on the formulation
    % in M. Tlidi et al., Chaos 27, 114312 (2017):
    % --> https://github.com/tony-ko/LLEF
    % The formalism for the FCGLE is based on L. Columbo, M. Piccardo et al.,
    % PRL 126, 173903 (2021).
    % The equation is solved by the split-step Fourier method for the linear
    % part, while the nonlinear part is by 2nd order Runge-Kutta method.
    %
    % The LLE equation with delayed feedback is:
    %
    % D[E(t,xi),t] = i*b*D[E(t,xi),xi,xi]
    %                -(1+i*th)*E(t,xi)+i*abs(E(t,xi))^2*E(t,xi)+E_in(t)
    %                +et*exp(i*ph)*E(t-tau,xi)
    %
    % et can be zero (regular LLE), constant (LLE w/feedback) or a function.
    % xi is fast time inside the cavity (subject to a periodic boundary
    % condition), t is slow evolution time.
    % The linear part is dispersion, losses, detuning, and pump:
    %
    % L = i*b*D[E(t,xi),xi,xi]-(1+i*th)*E(t,xi)+E_in(t)
    %
    % The nonlinear part is Kerr frequency shift and feedback:
    %
    % N = i*abs(E(t,xi))^2*E(t,xi) + et*exp(i*ph)*E(t-tau,xi)
    %
    % The FCGLE equation is:
    %
    % D[E(t,xi),t] = (1+i*G)*D[E(t,xi),xi,xi]
    %                +gamma*(1-i*Theta)*E(t,xi)-(1-i*Delta)*abs(E(t,xi))^2*E(t,xi)+E_in(t)
    %
    % G is related to diffusion and dispersion; gamma is a sign determining
    % above (+1) or below (-1) threshold operation; Theta depends on detuning,
    % Delta, and pumping; Delta depends on LEF and Kerr coefficient.
    % All parameters are defined in the PRL with the following renaming of the
    % variables: E --> F; t --> tau; E_in --> F_I; xi --> eta
    % Note that the field E, and the space-time coordinates xi and t are scaled.

    properties

        % Runtime, save, load and noise

        saveoutput = 0;
        loadoutput = 0;         % if 1, the initial condition is set by the last recorded state
        addnoise = 1;           % add noise to input state
        noiselevel = 0.005;
        spacetimescale = 0;     % plot output space and time coordinates in physical units
        solvertype = 1;         % 1 = FCGLE; 2 = CGLE; 3 = LLE

        run_simulation = false; % Toggle switch for running/pausing simulation

        % Space and time grids
            
        nF          % number of modes in Fourier space for the fast intracavity time
        dt          % slow time step duration
        nRK         % RK substeps
        
        % integration ranges
        
        S           % fast time
        endtime     % slow time
        
        % the further parameters are calculated upon object creation
        
        M           % number of steps to take in slow time 
        dtRK        % RK substep
        h           % fast time step
        n           % spectral mode indices   
        x           % fast time grid points
        k           % fast time wavenumbers
        
        
        % Default parameters for the equation
        
        % FCGLE parameters from L. Columbo, M. Piccardo et al., PRL 126, 173903 (2021)
        
        alpha = 2;      % LEF
        beta = 0;       % Kerr
        xi = 1;         % GVD (positive -> anomalous)
        
        G               % G = alpha + xi;
        b               % b = (1+1i*G)/1i;
         
        Theta = 4.7;
        gamma = 1;
        th              % th = (gamma*(-1+1i*Theta) - 1)/1i;
        
        et = 0;         % eliminates the delayed feedback term
        ph = pi;        % when et = 0, it becomes unrelevant
        tau = 100;      % when et = 0, it becomes unrelevant
        
        Y = 6.5;        % injected power
        E_in            % E_in = sqrt(Y);
        
        Delta           % Delta = alpha + beta;
        
        % CGLE parameters from M. Piccardo et al., Nature 582, 360 (2020)
        
        cD = -1.75;     % CGLE dispersive coefficient
        cNL = 0.7;      % CGLE nonlinear coefficient
        
        
        % Parameters required for rescaling space-time coordinates
        
        n_g = 3.3;          % group index
        tau_d = 60*1e-15;   % [s] dephasing time
        tau_p = 50*1e-12;   % [s] damping time of the cavity field
        r = 0.01;           % fractional pumping (>0 above threshold; <0 below threshold)
        c_light = 3e8;      % [m/s] light speed
        
        c_g
        d_r                 % [m^2]
        d_i                 % [m^2]

        % Current solution

        re_u
        im_u
        u
        U

        u_mean
        u_spatiotemporal
        u_spectrum
        u_spectrum_freq

        % plotting utilities
        UI_axis_cavity_intensity
        UI_axis_mean_intensity
        UI_axis_spacetime
        UI_axis_spectrum

        N_memory = 2000;
        colors

        

    end

    methods
        function obj = FCGLE_simulator()
            %FCGLE_simulator Construct an instance of this class
            %   Detailed explanation goes here
            obj.nF = 128;
            obj.dt = 0.05;
            obj.nRK = 16;       
            obj.S = 100;
            obj.endtime = 0.5;
            
            obj.M = round(obj.endtime/obj.dt);
            obj.dtRK = obj.dt/obj.nRK;
            obj.h = obj.S/obj.nF;
            obj.n = (-obj.nF/2:1:obj.nF/2-1)';
            obj.x = obj.n*obj.h;
            obj.k = 2*pi*obj.n/obj.S;
            
            switch obj.solvertype
                case 1
                    % FCGLE parameters from L. Columbo, M. Piccardo et al., PRL 126, 173903 (2021)
                    obj.alpha = 2; %LEF
                    obj.beta = 0; %Kerr
                    obj.xi = 1;
                    
                    obj.G = obj.alpha + obj.xi;
                    obj.b = (1+1i*obj.G)/1i;
                    
                    obj.Theta = 4.7;
                    obj.gamma = 1;
                    obj.th = (obj.gamma*(-1+1i*obj.Theta) - 1)/1i;
                    
                    obj.et = 0;  % eliminates the delayed feedback term
                    obj.ph = pi; % when et = 0, it becomes unrelevant
                    obj.tau = 1; % when et = 0, it becomes unrelevant
                    
                    obj.Y = 6; %injected power
                    obj.E_in = sqrt(obj.Y);
                    
                    obj.Delta = obj.alpha + obj.beta;
                case 2
                    % CGLE parameters
            
                    obj.cD = -1.5; % phase turbulence
                    obj.cNL = 1; % phase turbulence
%                     obj.cD = -3.5; % defect turbulence
%                     obj.cNL = 1; % defect turbulence
%                     obj.cD = -4.5; % defect turbulence
%                     obj.cNL = 1; % defect turbulence
%                     obj.cD = -1; % homoclones
%                     obj.cNL = 1.1; % homoclones
%                     obj.cD = 0.5; % Nozaki-Bekki holes
%                     obj.cNL = 2; % Nozaki-Bekki holes
%                     obj.cD = -1.1; % PRL
%                     obj.cNL = 1.1; % PRL
            
                    obj.G = -obj.cD;
                    obj.Delta = obj.cNL;
            
                    obj.b = (1+1i*obj.G)/1i;
                    
                    obj.Theta = 0;
                    obj.gamma = 1;
                    obj.th = (obj.gamma*(-1+1i*obj.Theta) - 1)/1i;
                    
                    obj.et = 0; % eliminates the delayed feedback term
                    obj.ph = pi; % when et = 0, it becomes unrelevant
                    obj.tau = 1; % when et = 0, it becomes unrelevant
                    
                    obj.Y = 0; %injected power
                    obj.E_in = sqrt(obj.Y);
                    
                    obj.beta = 0;
                    obj.alpha = obj.Delta - obj.beta;
                    obj.xi = obj.G - obj.alpha;
            end
        
            % Parameters required for rescaling space-time coordinates
        
            if obj.spacetimescale
                obj.n_g = 3.3;                                      % group index
                obj.tau_d = 60*1e-15;                               % [s] dephasing time
                obj.tau_p = 50*1e-12;                               % [s] damping time of the cavity field
                obj.r = 0.01;                                       % fractional pumping (>0 above threshold; <0 below threshold)
                obj.c_light = 3e8;                                  % [m/s] light speed
                
                obj.c_g = obj.c_light/obj.n_g;
                obj.d_r = (obj.c_g*obj.tau_d)^2/(1 + obj.alpha^2);  % [m^2]
                obj.d_i = obj.d_r*(obj.alpha + obj.xi);             % [m^2]
            end
            
%             if obj.loadoutput
%                 obj.re_u = importdata('re_u0.csv');
%                 obj.im_u = importdata('im_u0.csv');
%             else
%                 % obj.re_u = sech(obj.x/4);
%                 obj.re_u = ones(size(obj.x));
%                 obj.im_u = zeros(size(obj.x));
%             end
% 
%             obj.u = obj.re_u + 1i*obj.im_u;
% 
%             if obj.addnoise
%                 obj.u = obj.u + (randn(size(obj.x)) + 1i*randn(size(obj.x)))*obj.noiselevel;
%             end

        end
        
        function obj = initialize_solver(obj)
            %FCGLE_simulator Construct an instance of this class
            %   Detailed explanation goes here
            obj.nF = 2^10;
            obj.dt = 0.05;
            obj.nRK = 16;       
            obj.S = 100;
            obj.endtime = 0.5;
%             obj.endtime = 0.05;
            
            obj.M = round(obj.endtime/obj.dt);
            obj.dtRK = obj.dt/obj.nRK;
            obj.h = obj.S/obj.nF;
            obj.n = (-obj.nF/2:1:obj.nF/2-1)';
            obj.x = obj.n*obj.h;
            obj.k = 2*pi*obj.n/obj.S;
            
            switch obj.solvertype
                case 1 % FCGLE
                    % FCGLE parameters from L. Columbo, M. Piccardo et al., PRL 126, 173903 (2021)
                    obj.alpha = 2;              % LEF
                    obj.beta = 0;               % Kerr
                    obj.xi = 1;                 % GVD
                    
                    obj.G = obj.alpha + obj.xi;
                    obj.b = (1+1i*obj.G)/1i;
                    
                    obj.Theta = 4.7;            % detuning
                    obj.gamma = 1;             % 1 – active resonator above threshold, -1 – passive resonator
                    obj.th = (obj.gamma*(-1+1i*obj.Theta) - 1)/1i;
                    
                    obj.et = 0;                 % eliminates the delayed feedback term
                    obj.ph = pi;                % when et = 0, it becomes unrelevant
                    obj.tau = 1;                % when et = 0, it becomes unrelevant
                    
                    obj.Y = 6;                  %injected power
                    obj.E_in = sqrt(obj.Y);
                    
                    obj.Delta = obj.alpha + obj.beta;
                case 2 % CGLE
                    %CGLE parameters
            
            %         cD = -1.5; %phase turbulence
            %         cNL = 1; %phase turbulence
            %         cD = -3.5; %defect turbulence
            %         cNL = 1; %defect turbulence
            %         cD = -4.5; %defect turbulence
            %         cNL = 1; %defect turbulence
            %         cD = -1; %homoclons
            %         cNL = 1.1; %homoclons
                    obj.cD = 0.5; %Nozaki-Bekki holes
                    obj.cNL = 2; %Nozaki-Bekki holes
            
                    obj.G = obj.cD; %PRL with implicit sign "i" correction due to phys. vs. eng. notation of field
                    obj.Delta = -obj.cNL; %PRL with implicit sign "i" correction due to phys. vs. eng. notation of field
            
                    obj.b = (1+1i*obj.G)/1i;
                    
                    obj.Theta = 0;
                    obj.gamma = 1;
                    obj.th = (obj.gamma*(-1+1i*obj.Theta) - 1)/1i;
                    
                    obj.et = 0; % eliminates the delayed feedback term
                    obj.ph = pi; % when et = 0, it becomes unrelevant
                    obj.tau = 1; % when et = 0, it becomes unrelevant
                    
                    obj.Y = 0; %injected power
                    obj.E_in = sqrt(obj.Y);
                    
                    obj.beta = 0;
                    obj.alpha = obj.Delta - obj.beta;
                    obj.xi = obj.G - obj.alpha;

                case 3 % LLE
                    obj.b = 0;
                    obj.th = -1;

            end
        
            % Parameters required for rescaling space-time coordinates
        
            if obj.spacetimescale
                obj.n_g = 3.3;                                      % group index
                obj.tau_d = 60*1e-15;                               % [s] dephasing time
                obj.tau_p = 50*1e-12;                               % [s] damping time of the cavity field
                obj.r = 0.01;                                       % fractional pumping (>0 above threshold; <0 below threshold)
                obj.c_light = 3e8;                                  % [m/s] light speed
                
                obj.c_g = obj.c_light/obj.n_g;
                obj.d_r = (obj.c_g*obj.tau_d)^2/(1 + obj.alpha^2);  % [m^2]
                obj.d_i = obj.d_r*(obj.alpha + obj.xi);             % [m^2]
            end
            
            if obj.loadoutput
                obj.re_u = importdata('re_u0_512.csv');
                obj.im_u = importdata('im_u0_512.csv');
            else
                % obj.re_u = sech(obj.x/4);
                obj.re_u = ones(size(obj.x));
                obj.im_u = zeros(size(obj.x));
            end

            obj.u = obj.re_u + 1i*obj.im_u;
            obj.U = zeros(obj.M,obj.nF);

            if obj.addnoise
                obj.u = obj.u + (randn(size(obj.x)) + 1i*randn(size(obj.x)))*obj.noiselevel;
            end

            obj.u_mean = zeros(obj.N_memory,1);
            obj.u_spatiotemporal = zeros(obj.N_memory,length(obj.x));
%             addpath('utils/cbrewer')
%             cmap = cbrewer('div', 'RdBu', 512);  
%             cmap = cbrewer('seq', 'YlOrRd', 512); 
            load('lajolla.mat');
            cmap = lajolla;
            cmap(cmap<0) = 0;
            
            obj.colors = flipud(cmap);
        end
        
        function obj = advance_simulation(obj)
            %run_simulation Numerical integration of the equation by split-step Fourier method
%             if obj.run_simulation
%                 disp(obj.solvertype)
                % the matrix with all calc results
                obj.U = zeros(obj.M,obj.nF);
                
                % the number of previous iteration corresponding to the delay term
                nhist = ceil(obj.tau/obj.dt);
                
                fE_in=repmat(obj.E_in,obj.nF,1);
                fE_in=fftshift(fft(fE_in));
%                 u_curr = obj.re_u + 1i*obj.im_u;
                u_curr = obj.u; 
                
                hist_len = nhist*obj.nRK;
                
                init_hist = repmat({obj.u},hist_len,1);
                
                % our initial history for RK method
                qY1 = init_hist;
                qY1_head = 1;
                qY1_tail = hist_len;
                
                qY2 = init_hist;
                qY2_head = 1;
                qY2_tail = hist_len;
                
                clear init_hist
%                 tic
                % the loop in slow time
%                 disp(obj.Y)
                for m = 1:1:obj.M
%                     if mod(m,1000) == 0
%                         disp(m/obj.M*obj.endtime)
%                     end
                    % linear part regain only pump and losses terms
                    const_term = fE_in;
                    prop = -1i*obj.b*obj.k.^2 - (1 + 1i*obj.th);
%                     prop = -1i*obj.b*obj.k.^2  +0.005 * obj.b * obj.k.^3 - (1 + 1i*obj.th);
                    exp_prop = exp(prop*obj.dt/2);
                    % this can seem incoherent, but this is how you introduce a constant
                    % term in SSFM
                    const_add = const_term.*(exp_prop - 1)./prop;
                    
                    % we propagate the solution half-step in Fourier domain
                    v = fftshift(fft(u_curr));
                    v = const_add + exp_prop.*v;
                    v = ifft(fftshift(v));                    
                    % gradual feedback increase
                %     et = 0.02*fix(m/round(3000/dt)); %It's crucial to comment out this or there will be feedback after some time!
                    obj.et = 0;
                    
                    % we propagate the nonlinear part by RK substeps
                    % we use queue structure to keep intermediate steps
                    
                    
                    for l=1:1:obj.nRK
                        %%boresome operators over queues
                        %%q 1
                %         Y1T = qY1{qY1_head};
                %         qY1_head = qY1_head + 1;
                %         if qY1_head > hist_len
                %             qY1_head = 1;
                %         end
                        %%q 2
                %         Y2T = qY2{qY2_head};
                %         qY2_head = qY2_head + 1;
                %         if qY2_head > hist_len
                %             qY2_head = 1;
                %         end
                        
                        %%yay, some RK calculations!
                        Y1T = 0;
                        Y2 = v + 0.5.*obj.dtRK.*FCGLEwF_LLN(v,Y1T,obj.et,obj.ph,obj.Delta);
                        
                        %%the same boredom for the queues tails
                %         qY1_tail = mod(qY1_tail, hist_len) + 1;
                %         qY1{qY1_tail} = v;
                %  
                %         qY2_tail = mod(qY2_tail, hist_len) + 1;
                %         qY2{qY2_tail} = Y2;
                        Y2T = 0;
                        f2 = FCGLEwF_LLN(Y2,Y2T,obj.et,obj.ph,obj.Delta);
                        v = v +obj.dtRK.*f2;
                    end
                    
                    % then propagate it half-step further
                    v = fftshift(fft(v));
                    v = const_add + exp_prop.*v;
                    u_curr = ifft(fftshift(v));
                      
                    % add to the field storage
                    obj.U(m,:) = u_curr'+ ((randn(size(obj.x)) + 1i*randn(size(obj.x)))*obj.noiselevel)';
%                     pause(0.01)
%                     plot(obj.UI_axis_cavity_intensity,abs(U(m,:)).^2)
                end
                pause(0.01)

                obj = obj.update_plots();
                
%                 obj.u = U(end,:)';
%                 obj.re_u = real(U(end,:))';
%                 obj.im_u = imag(U(end,:))';
% 
%                 plot(obj.UI_axis_cavity_intensity,obj.x,abs(U(end,:)).^2,'Color',[0 0 0],'Linewidth',1.5)
%                 obj.UI_axis_cavity_intensity.YLim = [0 5];
% 
%                 obj.u_mean = circshift(obj.u_mean,-1);
%                 obj.u_mean(end) = mean(abs(U(end,:)).^2);
%                 plot(obj.UI_axis_mean_intensity,obj.u_mean,'Color',[0 0 0],'Linewidth',1.5)
%                 obj.UI_axis_mean_intensity.YLim = [0 5];
% 
%                 obj.u_spatiotemporal = circshift(obj.u_spatiotemporal,-1,1);
%                 obj.u_spatiotemporal(end,:) = U(end,:);
% 
%                 imagesc(obj.UI_axis_spacetime,obj.x,linspace(1,obj.N_memory,obj.N_memory),abs(obj.u_spatiotemporal).^2)
% %                 obj.UI_axis_spacetime.XLim = [obj.x(1) obj.x(end)];
%                 obj.UI_axis_spacetime.CLim = [0 5];
%                 obj.UI_axis_spacetime.Colormap = obj.colors;
%                 
%                 N_rt = 6;
%                 L = N_rt*length(obj.x);
%                 u_for_fft = reshape(U',1,size(U,1).*size(U,2));
% %                 u_for_fft = reshape(obj.u_spatiotemporal',[1 size(obj.u_spatiotemporal,1)*size(obj.u_spatiotemporal,2)]);
% %                 u_for_fft = u_for_fft.*blackman(1,length(u_for_fft));
% %                 u_for_fft = abs(u_for_fft);
%                 T = abs(obj.x(2)-obj.x(1));
%                 Fs = 1/T;
%                 f = Fs*(-L/2:L/2)/L*(abs(obj.x(end)-obj.x(1)));
%                 spectrum = fftshift(fft(u_for_fft(end-L:end)));
%                 obj.u_spectrum_freq = f;
%                 obj.u_spectrum = spectrum;
%                 
% %                 spectrum = circhift(spectrum,floor(length(spectrum)/2));
%                 plot(obj.UI_axis_spectrum,f,10*log10(abs(spectrum)),'Color',[0 0 0],'Linewidth',1.0)
%                 obj.UI_axis_spectrum.XLim = [-100 100];
%                 obj.UI_axis_spectrum.YLim = [-10 60];

                % clear qY1 qY2
                % clear qY1_tail qY2_tail
                % clear qY1_head qY2_head
                % clear f1 f2 exp_prop prop const_add const_term
                % clear Y1 Y1T Y2 Y2T u0
%                 toc
                
                %% Export the last roundtrip to use it further
%                 if obj.saveoutput
%                     u0 = U(end,:)';
%                     dlmwrite('re_u0.csv',real(u0),'delimiter',',','precision',5);
%                     dlmwrite('im_u0.csv',imag(u0),'delimiter',',','precision',5);
%                 end
%             end
        end
    
        function obj = update_params(obj)
            
            obj.E_in = sqrt(obj.Y);
            obj.th = (obj.gamma*(-1+1i*obj.Theta) - 1)/1i;
            obj.G = obj.alpha + obj.xi;
            obj.b = (1+1i*obj.G)/1i;
            obj.Delta = obj.alpha + obj.beta;
        end

        function obj = update_plots(obj)
                obj.u = obj.U(end,:)';
                obj.re_u = real(obj.U(end,:))';
                obj.im_u = imag(obj.U(end,:))';

                plot(obj.UI_axis_cavity_intensity,obj.x,abs(obj.U(end,:)).^2,'Color',[0 0 0],'Linewidth',1.5)
                obj.UI_axis_cavity_intensity.YLim = [0 5];

                obj.u_mean = circshift(obj.u_mean,-1);
                obj.u_mean(end) = mean(abs(obj.U(end,:)).^2);
                plot(obj.UI_axis_mean_intensity,obj.u_mean,'Color',[0 0 0],'Linewidth',1.5)
                obj.UI_axis_mean_intensity.YLim = [0 5];

                obj.u_spatiotemporal = circshift(obj.u_spatiotemporal,-1,1);
                obj.u_spatiotemporal(end,:) = obj.U(end,:);

                imagesc(obj.UI_axis_spacetime,obj.x,linspace(1,obj.N_memory,obj.N_memory),abs(obj.u_spatiotemporal).^2)
%                 obj.UI_axis_spacetime.XLim = [obj.x(1) obj.x(end)];
                obj.UI_axis_spacetime.YLim = [0 obj.N_memory];
                obj.UI_axis_spacetime.YTick = [0 obj.N_memory];
                obj.UI_axis_spacetime.YTickLabel = {'0' num2str(obj.N_memory)};
                obj.UI_axis_spacetime.CLim = [0 5];
                obj.UI_axis_spacetime.Colormap = obj.colors;
                
                N_rt = 6;
%                 N_rt = 1;
                L = N_rt*length(obj.x);
                u_for_fft = reshape(obj.U',1,size(obj.U,1).*size(obj.U,2));
%                 u_for_fft = reshape(obj.u_spatiotemporal',[1 size(obj.u_spatiotemporal,1)*size(obj.u_spatiotemporal,2)]);
%                 u_for_fft = u_for_fft.*blackman(1,length(u_for_fft));
%                 u_for_fft = abs(u_for_fft);
                T = abs(obj.x(2)-obj.x(1));
                Fs = 1/T;
                f = Fs*(-L/2:L/2)/L*(abs(obj.x(end)-obj.x(1)));
                spectrum = fftshift(fft(u_for_fft(end-L:end)));
                obj.u_spectrum_freq = f;
                obj.u_spectrum = spectrum;
                
%                 spectrum = circhift(spectrum,floor(length(spectrum)/2));
                plot(obj.UI_axis_spectrum,obj.u_spectrum_freq,10*log10(abs(obj.u_spectrum)),'Color',[0 0 0],'Linewidth',1.0)
                obj.UI_axis_spectrum.XLim = [-100 100];
                obj.UI_axis_spectrum.YLim = [-10 60];
        end
    
        function obj = changeMemorySize(obj,N)
            obj.N_memory = N;
            obj.u_mean = zeros(obj.N_memory,1);
            obj.u_spatiotemporal = zeros(obj.N_memory,length(obj.x));
        end
        
        function dy = LLwF_LLN(y,yT,et,ph,Delta)

            %LLwF_LLN Lugiato-Lefever nonlinear RHS (FCGLE)
            dy = -(1 - 1i*Delta).*abs(y).^2.*y + et.*exp(-1i.*ph).*yT;

        end

        
    end
end
