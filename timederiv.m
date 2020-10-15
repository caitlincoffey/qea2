function dTdt=timederiv(t,T,L_absorber, L_wall)
% Equation: Cabsorber * dTroom/dt = q*Awindow - (Troom - Tout)/Rtotal

% temperature, areas, thickness, solar flux
T_out = -3; % degrees celsius, outside
A_wall = 99.6; % m^2, surface area of all walls
A_window = 2.6*5; %m^2, surface area of window
% decrease surface area of the window! 
A_floor = 5.1*5*2; %m^2 we are taking top and bottom surfaces of insulator!
%A_floor = 5.1*5*2 + L_absorber*5.1*2 + L_absorber*5*2;
q = -361*cos(pi*t/(12*3600)) + 224*cos(pi*t/(6*3600)) + 210; % q in W/m^2, t in seconds

% extra terms for conductivity: k/(L*A), convection: 1/(h*A)
k_wall = .04; % thermal conductivity walls, in W/m-K
h_window = .7; % in W/m-K, equivalent thermal conductivity of glass window
% note: we hardcode h_inside and h_outside into our equations :)

c = 800; %J/kg-K specific heat
density = 3000; % kg/m^3
volume = L_absorber*A_floor; %m^3
Ca = density*c*volume; % in J/K, heat capacity of the absorber

% convection terms for R window and R wall (in series for Rtotal):
convection_floor_innerair = 1/(15*A_floor);

% convection terms for R window:
convection_innerair_insidesurface = 1/(15*A_window); % inside air -> inside surface
convection_insidesurface_outsidesurface = 1/(h_window*A_window);
convection_outsidesurface_outerair = 1/(30*A_window); 

% convection and conduction terms for R wall:
convection_indoor_wall = 1/(15*A_wall);
conduction_wall = k_wall/(L_wall*A_wall); 
convection_outdoor_wall = 1/(30*A_wall);

% Resitances  
R_window = convection_innerair_insidesurface + convection_insidesurface_outsidesurface + convection_outsidesurface_outerair;
% no conduction because there's so many layers! the vacuum is so small
% (fluid involved with contact with surfaces). So we are using a convection
% term
R_walls = convection_indoor_wall + convection_outdoor_wall + conduction_wall;

Rtotal = convection_floor_innerair + 1/(1/R_walls + 1/R_window); % we know Rtotal is too small, making dTdt too big
% R_walls gets smaller with larger wall thickness values, Rtotal gets
% smaller. Means that dTdt gets bigger?
dTdt= (q*A_window)/Ca - ((T - T_out)/Rtotal)/Ca;
end
