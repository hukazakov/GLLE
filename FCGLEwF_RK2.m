%% Forced complex Ginzburg-Landau equation solver
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
%

clear all
close all

%% here we define the space and time grids
% number of modes in Fourier space for the fast intracavity time
nF = 512;
% slow time step duration
dt = 0.05;
% RK substeps
nRK = 16;

% integration ranges
% fast time
S = 50;
% slow time
endtime = 500;

% the further parameters are calculated
% number of steps to take in slow time
M = round(endtime/dt);
% RK substep
dtRK = dt/nRK;

% fast time step
h = S/nF;
% spectral mode indices
n = (-nF/2:1:nF/2-1)';
% fast time grid points
x = n*h;
% fast time wavenumbers
k = 2*pi*n/S;

%% Save, load and noise

saveoutput = 1;
loadoutput = 1; % if 1, the initial condition is set by the last recorded state
addnoise = 1; %add noise to input state
noiselevel = 0.01;
spacetimescale = 0; %plot output space and time coordinates in physical units

%% Parameters for the equation

%FCGLE parameters from L. Columbo, M. Piccardo et al., PRL 126, 173903 (2021)
alpha = 2; %LEF
beta = 0; %Kerr
xi = 1;

G = alpha + xi;
b = (1+1i*G)/1i;

Theta = 4.7;
gamma = 1;
th = (gamma*(-1+1i*Theta) - 1)/1i;

et = 0; % eliminates the delayed feedback term
ph = pi; % when et = 0, it becomes unrelevant
tau = 100; % when et = 0, it becomes unrelevant

Y = 6.5; %injected power
E_in = sqrt(Y);

Delta = alpha + beta;

%CGLE parameters from M. Piccardo et al., Nature 582, 360 (2020)
% cD = -1.75;
% G = cD;
% b = (1+1i*G)/1i;
% 
% Theta = 0;
% gamma = 1;
% th = (gamma*(-1+1i*Theta) - 1)/1i;
% 
% et = 0; % eliminates the delayed feedback term
% ph = pi; % when et = 0, it becomes unrelevant
% tau = 100; % when et = 0, it becomes unrelevant
% 
% Y = 0; %injected power
% E_in = sqrt(Y);
% 
% cNL = 0.7;
% Delta = -cNL;

%% Parameters required for rescaling space-time coordinates

n_g = 3.3; % group index
tau_d = 60*1e-15; %[s] dephasing time
tau_p = 50*1e-12; %[s] damping time of the cavity field
r = 0.1; % fractional pumping (>0 above threshold; <0 below threshold)
c_light = 3e8; %[m/s] light speed

c_g = c_light/n_g;
d_r = (c_g*tau_d)^2/(1 + alpha^2); %[m^2]
d_i = d_r*(alpha + xi); %[m^2]


%% Numerical integration of the equation by split-step Fourier method

if loadoutput
    % or we can import the precomputed cavity soliton solution
    re_u0 = importdata('re_u0_512.csv');
    im_u0 = importdata('im_u0_512.csv');
    u0 = re_u0 + 1i*im_u0;
    if addnoise
        u0 = u0 + rand(size(x))*noiselevel;
    end
    clear re_u0 im_u0
else
    % initial condition u(t;t<=0,xi), and also this is the history function
    % it can be just some sech pulse
    % u0 = sech(x/4);
    u0 = ones(size(x));
    if addnoise
        u0 = u0 + rand(size(x))*noiselevel;
    end
end



% the matrix with all calc results
U = zeros(M,nF);

% the number of previous iteration corresponding to the delay term
nhist = ceil(tau/dt);

fE_in=repmat(E_in,nF,1);
fE_in=fftshift(fft(fE_in));
u = u0;

hist_len = nhist*nRK;

init_hist = repmat({u0},hist_len,1);

% our initial history for RK method
qY1 = init_hist;
qY1_head = 1;
qY1_tail = hist_len;

qY2 = init_hist;
qY2_head = 1;
qY2_tail = hist_len;

clear init_hist
tic
% the loop in slow time
for m = 1:1:M
    if mod(m,1000) == 0
        disp(m)
    end
    % linear part regain only pump and losses terms
    const_term = fE_in;
    prop = -1i*b*k.^2 - (1 + 1i*th);
    exp_prop = exp(prop*dt/2);
    % this can seem incoherent, but this is how you introduce a constant
    % term in SSFM
    const_add = const_term.*(exp_prop - 1)./prop;
    
    % we propagate the solution half-step in Fourier domain
    v = fftshift(fft(u));
    v = const_add + exp_prop.*v;
    v = ifft(fftshift(v));
    
    % gradual feedback increase
    et = 0.02*fix(m/round(3000/dt));
    
    % we propagate the nonlinear part by RK substeps
    % we use queue structure to keep intermediate steps
    
    
    for l=1:1:nRK
        % boresome operators over queues
        % q 1
        Y1T = qY1{qY1_head};
        qY1_head = qY1_head + 1;
        if qY1_head > hist_len
            qY1_head = 1;
        end
        
        % q 2
        Y2T = qY2{qY2_head};
        qY2_head = qY2_head + 1;
        if qY2_head > hist_len
            qY2_head = 1;
        end
        
        % yay, some RK calculations!
        Y2 = v + 0.5.*dtRK.*FCGLEwF_LLN(v,Y1T,et,ph,Delta);
        
        % the same boredom for the queues tails
        qY1_tail = mod(qY1_tail, hist_len) + 1;
        qY1{qY1_tail} = v;
 
        qY2_tail = mod(qY2_tail, hist_len) + 1;
        qY2{qY2_tail} = Y2;
        
        f2 = FCGLEwF_LLN(Y2,Y2T,et,ph,Delta);
        v = v +dtRK.*f2;
    end
    
    % then propagate it half-step further
    v = fftshift(fft(v));
    v = const_add + exp_prop.*v;
    u = ifft(fftshift(v));
      
    % add to the field storage
    U(m,:) = u';
end

clear qY1 qY2
clear qY1_tail qY2_tail
clear qY1_head qY2_head
clear f1 f2 exp_prop prop const_add const_term
clear Y1 Y1T Y2 Y2T u0
toc

%% Export the last roundtrip to use it further
if saveoutput
    u0 = U(end,:)';
    dlmwrite('re_u0_512.csv',real(u0),'delimiter',',','precision',5);
    dlmwrite('im_u0_512.csv',imag(u0),'delimiter',',','precision',5);
end

%% Stacked plot
figure
xscale = 3000;
yscale = 1;
h = waterfall(abs(U(1:xscale:end,1:yscale:end)).^2);
colormap([1 1 1]);
color = repmat([0 0 0],length(h.FaceVertexCData),1);
set(h, 'FaceVertexCData', color)
set(h, 'FaceColor', 'flat')
set(gca,'Color','k')
set(h, 'EdgeColor', [1 1 1])
view([-30 60])
axis tight
zlabel('|E|^2')
if spacetimescale
    S_scaled = S/sqrt(r/d_r)*1e6; %[um]
    endtime_scaled = endtime*tau_p/abs(r)*1e9; %[ns]
    set(gca, 'xtick', 1:nF/(10*yscale):nF)
    set(gca, 'xticklabels', round(-S_scaled/2:S_scaled/10:S_scaled/2,0))
    set(gca, 'ytick', 1:M/(10*xscale):M/xscale)
    set(gca, 'yticklabels', round(0:endtime_scaled/10:endtime_scaled,0))
    xlabel('x (\mum)')
    ylabel('t (ns)')
else
    set(gca, 'xtick', 1:nF/(10*yscale):nF)
    set(gca, 'xticklabels', -S/2:S/10:S/2)
    set(gca, 'ytick', 1:M/(10*xscale):M/xscale)
    set(gca, 'yticklabels', 0:endtime/10:endtime)
    xlabel('x (a.u.)')
    ylabel('t (a.u.)')
end

clear h xscale yscale color

%% Projected surface plot
figure
xscale = 5;
yscale = 1;
imagesc(abs(U(1:xscale:end,1:yscale:end)).^2)
h = colorbar;
ylabel(h, '|E|^2')
axis tight
if spacetimescale
    S_scaled = S/sqrt(r/d_r)*1e6; %[um]
    endtime_scaled = endtime*tau_p/abs(r)*1e9; %[ns]
    set(gca, 'xtick', 1:nF/(10*yscale):nF)
    set(gca, 'xticklabels', round(-S_scaled/2:S_scaled/10:S_scaled/2,0))
    set(gca, 'ytick', 1:M/(10*xscale):M/xscale)
    set(gca, 'yticklabels', round(0:endtime_scaled/10:endtime_scaled,0))
    xlabel('x (\mum)')
    ylabel('t (ns)')
else
    set(gca, 'xtick', 1:nF/(10*yscale):nF)
    set(gca, 'xticklabels', -S/2:S/10:S/2)
    set(gca, 'ytick', 1:M/(10*xscale):M/xscale)
    set(gca, 'yticklabels', 0:endtime/10:endtime)
    xlabel('x (a.u.)')
    ylabel('t (a.u.)')
end

clear xscale yscale h

