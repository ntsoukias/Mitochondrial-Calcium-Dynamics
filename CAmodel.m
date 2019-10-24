function [t, c, e, m, u, NonStateVar] = CAmodel(varargin)

dbstop if error

defaultClampJIP3R = 'false';
defaultClampJIP3R_u = 'false';
defaultRemoveMt = 'false';
defaultRemoveMd = 'false';
defaultTspan = [0, 600];
defaultStimOnset = 100;

p = inputParser;
addRequired(p, 'Pset', @isstruct)
addParameter(p, 'ClampJIP3R', defaultClampJIP3R)
addParameter(p, 'ClampJIP3R_u', defaultClampJIP3R_u)
addParameter(p, 'RemoveMt', defaultRemoveMt)
addParameter(p, 'RemoveMd', defaultRemoveMd)
addParameter(p, 'Tspan', defaultTspan)
addParameter(p, 'StimOnset', defaultStimOnset)

parse(p, varargin{:})

valueClampJIP3R = p.Results.ClampJIP3R;
valueClampJIP3R_u = p.Results.ClampJIP3R_u;
flagRemoveMt = p.Results.RemoveMt;
flagRemoveMd = p.Results.RemoveMd;

P = p.Results.Pset;
tspan = p.Results.Tspan;
stimonset = p.Results.StimOnset;

volCt = P.volCt;      % [pL] volume cytosol

% volMd = P.volMd;      % [pL] volume microdomain

distance = P.distance; % [nm] mean intermembrane distance

dist = distance*1e-9; % [m]
r = 0.58*1e-6; % [m]
SA = 2*0.2*4*pi*r^2; % [m^2]
nMtObject = 200;
volMd = SA*nMtObject*dist; % [m^3]
volMd = volMd * 1e15; % [pL]

volER = P.volER;      % [pL] volume ER
volMt = P.volMt;      % [pL] volume mitocondria
Vip3r = P.Vip3r;      % [1/s] max flux of IP3R
Vserca = P.Vserca;    % [uM/s] max flux of SERCA pump
kserca = P.kserca;    % [uM] activation constant for SERCA pump
ip3 = P.ip3;          % [uM] IP3 in the cytosol
a2 = P.a2;            % [uM^-1*s^-1] IP3R binding rate ca inhibition sites
d1 = P.d1;            % [uM] IP3R dissociation constant for IP3 sites
d2 = P.d2;            % [uM] IP3R dissociation constant for ca inhibition sites
d3 = P.d3;            % [uM] IP3R dissociation constant for IP3 sites
d5 = P.d5;            % [uM] IP3R dissociation constant for ca activation sites
Vmcu = P.Vmcu;        % [uM/s] max rate of ca uptake by MCU
kmcu = P.kmcu;        % [uM] half-max rate of ca pumping from c to m
Vncx = P.Vncx;        % [uM/s] max rate of ca release through NCX
kncx = P.kncx;        % [uM] activation constant for NCX
kna = P.kna;          % [mM] Na activation constant for MCU
N = P.N;              % [mM] Na in cytosol
N_u = P.N_u;          % [mM] Na in microdomain
leak_e_u = P.leak_e_u; % [1/s] leak constant from ER to Md
leak_e_c = P.leak_e_c; % [1/s] leak constant from ER to Ct
leak_u_c = P.leak_u_c; % [1/s] leak constant from Md to Ct
cI = P.cI;            % fraction of IP3R facing microdomain
cS = P.cS;            % fraction of SERCA facing microdomain
cM = P.cM;            % fraction of MCU facing microdomain
cN = P.cN;            % fraction of mNCX facing microdomain
bt_c = P.bt_c;      % [uM] total buffer concentration in cytosol
K_c = P.K_c;        % buffer rate constant ratio (Qi 2015)
bt_e = P.bt_e;      % [uM] total buffer concentration in ER
K_e = P.K_e;        % buffer rate constant ratio (Qi 2015)
bt_m = P.bt_m;      % [uM] total buffer concentration in mitocondria
K_m = P.K_m;        % buffer rate constant ratio (Qi 2015)
bt_u = P.bt_u;      % [uM] total buffer concentration in micro-domain
K_u = P.K_u;        % buffer rate constant ratio (Qi 2015)

%% Initial Conditions
cInit = 0.1079; % [uM]
eInit = 247.4026; % [uM]
mInit = 0.1845; % [uM]
uInit = 0.2384; % [uM]

hInit = 0.7;
h_uInit = 0.7;

initStateVar = [cInit; eInit; mInit; uInit; hInit; h_uInit];
X0 = initStateVar;
FSoptions = optimset('Display','off');
X0 = fsolve(@(X) equations(0,X), X0, FSoptions);

%% Solve
[t,X] = ode45(@equations, tspan, X0);

%% Assign solutions to state variables
k = 1;
c = X(:,k);     k = k + 1;
e = X(:,k);     k = k + 1;
m = X(:,k);     k = k + 1;
u = X(:,k);     k = k + 1;
h = X(:,k);     k = k + 1;
h_u = X(:,k);   k = k + 1;

%% Compute non-state
[~, ns] = equations(t, X.');

k = 1;
NonStateVar.Jip3r = ns(k,:);     k = k + 1;
NonStateVar.Jip3r_u = ns(k,:);   k = k + 1;
NonStateVar.Jserca = ns(k,:);    k = k + 1;
NonStateVar.Jserca_u = ns(k,:);  k = k + 1;
NonStateVar.Jncx = ns(k,:);      k = k + 1;
NonStateVar.Jncx_u = ns(k,:);    k = k + 1;
NonStateVar.Jmcu = ns(k,:);      k = k + 1;
NonStateVar.Jmcu_u = ns(k,:);    k = k + 1;
NonStateVar.Jleak_u_c = ns(k,:); k = k + 1;
NonStateVar.Jleak_e_u = ns(k,:); k = k + 1;
NonStateVar.Jleak_e_c = ns(k,:); k = k + 1;
NonStateVar.theta_c = ns(k,:);   k = k + 1;
NonStateVar.theta_e = ns(k,:);   k = k + 1;
NonStateVar.theta_m = ns(k,:);   k = k + 1;
NonStateVar.theta_u = ns(k,:);   k = k + 1;
NonStateVar.Poip3r = ns(k,:);    k = k + 1;
NonStateVar.Poip3r_u = ns(k,:);

%% Setup differential equations to be solved using ODE solver

    function [dXdt, nonState] = equations(t,X)
        
        if t <= stimonset
            ip3 = 0;
        else
            ip3 = P.ip3;
        end
        
        % State Variables
        k = 1;
        c = X(k, :);     k = k + 1;
        e = X(k, :);     k = k + 1;
        m = X(k, :);     k = k + 1;
        u = X(k, :);     k = k + 1;
        h = X(k, :);     k = k + 1;
        h_u = X(k, :);   k = k + 1;

        %% Compute Non-State Variables
        % IP3R
        sact = h.*((ip3/(ip3 + d1)) .* c./(c + d5));
        Poip3r = sact.^4 + 4*sact.^3.*(1 - sact);
        
        % IP3R_u
        sact_u = h.*((ip3/(ip3 + d1)) .* u./(u + d5));
        Poip3r_u = sact_u.^4 + 4*sact_u.^3.*(1 - sact_u);

        if t >= stimonset
            if isnumeric(valueClampJIP3R) && t >= stimonset
                Jip3r = valueClampJIP3R;
            else
                Jip3r = (1 - cI)*(Vip3r*Poip3r).*(e - c);
            end
            
            if isnumeric(valueClampJIP3R_u) && t >= stimonset
                Jip3r_u = valueClampJIP3R_u;
            else
                Jip3r_u = cI*(Vip3r*Poip3r_u).*(e - u);
            end
        else
            Jip3r = (1 - cI)*(Vip3r*Poip3r).*(e - c);
            Jip3r_u = cI*(Vip3r*Poip3r_u).*(e - u);
        end
        
        % SERCA
        Jserca = (1 - cS)*Vserca*c.^2./(kserca^2 + c.^2); 
        % SERCA_u
        Jserca_u = cS*Vserca*u.^2./(kserca^2 + u.^2);
        % mNCX
        Jncx = (1 - cN)*Vncx*(N^3/(kna^3 + N^3)).*(m./(kncx + m));    
        % mNCX_u
        Jncx_u = cN*Vncx*(N_u^3/(kna^3 + N_u^3)).*(m./(kncx + m));
        % MCU
        Jmcu = (1 - cM)*Vmcu*(c.^2./(kmcu^2 + c.^2));       
        % MCU_u
        Jmcu_u = cM*Vmcu*(u.^2./(kmcu^2 + u.^2));
        % leaks
        Jleak_u_c = leak_u_c.*(u - c);
        Jleak_e_u = leak_e_u.*(e - u);
        Jleak_e_c = leak_e_c.*(e - c);        
        % h
        ah = a2*d2*(ip3 + d1)/(ip3 + d3);
        bh = a2*c;        
        % h_u
        bh_u = a2*u;        
        % buffering
        theta_c = bt_c*K_c./((K_c + c).^2); % buffer factor
        theta_e = bt_e*K_e./((K_e + e).^2); % buffer factor
        theta_m = bt_m*K_m./((K_m + m).^2); % buffer factor
        theta_u = bt_u*K_u./((K_u + u).^2); % buffer factor   
        
        if strcmp(flagRemoveMt, 'true')
            Jmcu = zeros(1, length(theta_u));
            Jmcu_u = zeros(1, length(theta_u));
            Jncx = zeros(1, length(theta_u));
            Jncx_u = zeros(1, length(theta_u));
        end
        
        if strcmp(flagRemoveMd, 'true')
            Jleak_u_c = zeros(1, length(theta_u));
            Jleak_e_u = zeros(1, length(theta_u));
            Jip3r_u = zeros(1, length(theta_u));
            Jserca_u = zeros(1, length(theta_u));
            Jmcu_u = zeros(1, length(theta_u));
            Jncx_u = zeros(1, length(theta_u));
            h_u = ones(1, length(theta_u));
            bh_u = zeros(1, length(theta_u));
        end
        
        %% Compute State Variables
        dcdt = (Jip3r + Jleak_u_c + Jleak_e_c + Jncx...
            - Jserca - Jmcu)./(1 + theta_c);
        
        dedt = (volCt/volER*(Jserca + Jserca_u ...
            - Jip3r - Jip3r_u - Jleak_e_u - Jleak_e_c))./(1 + theta_e);
        
        dmdt = (volCt/volMt*(Jmcu + Jmcu_u - Jncx...
            - Jncx_u))./(1 + theta_m);
        
        dudt = (volCt/volMd*(Jip3r_u + Jncx_u + Jleak_e_u ...
            - Jserca_u - Jmcu_u - Jleak_u_c))./(1 + theta_u);
        
        dhdt = ah.*(1 - h) - bh.*h;
        
        dh_udt = ah.*(1 - h_u) - bh_u.*h_u;

        %% Assign equations to function output
        dXdt = zeros(size(X));
        
        % State Variables
        k = 1;
        dXdt(k, :) = dcdt;     k = k + 1;
        dXdt(k, :) = dedt;     k = k + 1;
        dXdt(k, :) = dmdt;     k = k + 1;
        dXdt(k, :) = dudt;     k = k + 1;
        dXdt(k, :) = dhdt;     k = k + 1;
        dXdt(k, :) = dh_udt;
         
        nonState = [Jip3r; Jip3r_u; Jserca; Jserca_u; Jncx; Jncx_u; ...
            Jmcu; Jmcu_u; Jleak_u_c; Jleak_e_u; Jleak_e_c; theta_c; ...
            theta_e; theta_m; theta_u; Poip3r; Poip3r_u];

    end

end

