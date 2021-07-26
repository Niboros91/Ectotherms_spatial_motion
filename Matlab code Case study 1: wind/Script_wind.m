% ODE model for Drosophila suzukii
% Time varying (based on temperature)

%Initialization
clear, clc % To clear the workspace and the command function

%% Tunning parameters
% Wind
k_factor = 1/100; %
saturation = 5; %Wind at which there is no movement
%% PARAMETERS functions

%Growth rate (Drosophila suzukii)
 a = 1.2*(10^(-4));
 T_l = 3;
 T_m = 30;
 m = 6;
%  

%Mortality rate (Drosophila suzukii)
a1 = -5.4E-06;
b1 = 0.0005194;
c1 = -0.0116827;
d1 = 2.16E-05;
e1 = 1.3146586;
% 

% Birth rate (Drosophila suzukii)
alpha = 659.06;
gamma = 88.53;
lambda = 52.32;
delta = 6.06;
tau = 22.87;

% Sex ratio
s_r = 0.5;

% Mating ratio
r_remate = 0;
r_mate = 1;

% Additional parameters
N_stages = 8; % Eggs/3 larva stages/ pupa /male/ unmated female/ mated female
N_trees = 16; %Number of tress

%% DATA

%Load additional data. 

%Temperature
load('data_pantheon_2020'); %Temperature is Temp_avg
temperature(277:288) = 15; %Temperature from 277 to 288 missing (malfunctioning of weather station)
temperature = temperature(92:end); %Selecting data fromt the 1st of April
wind_direction = wind_direction(92:end);
wind_speed = wind_speed(92:end);


%Small temperature variation for each area
max_t = 1.5;
min_t = -1.5;
t_var = min_t + (max_t-min_t) .* rand(length(temperature'),N_trees);
Temp_areas = repmat(temperature,N_trees,1);
Temp_areas = Temp_areas + t_var';

check = find(Temp_areas >30);
Temp_areas(check) = 29.8;

% Wind direction is degree angles with respect to the north
wind_x = wind_speed .* cosd(wind_direction);
wind_y = wind_speed .* sind(wind_direction);
w=[wind_x', wind_y'];

% Divide by the a correcting factor
w(find(abs(w)>saturation)) = 0; %If wind more than the saturation limit 0.
w(isnan(w))=0; %If there is no measurement from the weather station. Wind 0
w = w*k_factor;

%% Compute rates
%Compute the rates based on temperature
for i =1: N_trees
    G_R(i) = growth_rate(Temp_areas(i,1),a,T_l,T_m,m); %Calling the growth rate function
    B_R(i) = birth_rate_suzuki(Temp_areas(i,1),alpha,gamma,lambda,delta,tau);%Calling the birth rate function
    M_R(i) = mortality_rate(Temp_areas(i,1),a1,b1,c1,d1,e1); %Calling the mortality rate function  
    S_R(i) = s_r;
    R_remate(i) = r_remate;
    R_mate(i) = r_mate;
    
end

%Initialize stages
Pest_stages(1:N_stages,N_trees) = stage; % We create an array of the stage class
Pest_stages = Initialize_stages_ode(B_R,M_R,G_R,S_R,R_mate,R_remate,Pest_stages,N_trees,N_stages); %We initialize the parameters associated to each stage class

%Spatial configuration
[Adj,Adj_w,Link]=Adjency_matrix(4,4);


%% DYNAMIC MODEL

%Initial values for the variables
x = zeros(N_stages*N_trees,1); %Number of total states
A_cont = compute_A_multi_16(N_trees,Pest_stages,N_stages,Adj_w,w(1,:));

%A_cont =compute_A_continous(Pest_stages); %We build the continuous time A matrix
% Discretization of the system based on the zero order holder
sysc = ss(A_cont,[],eye(N_stages*N_trees),[]);
sysd = c2d(sysc,1,'zoh');
A_dis = sysd.A; %Discretized matrix A

%% SIMULATIONS
Simulation_time =length(temperature); %Simulation lenght based on the temperature array introduced

%Plants initially infected
inif = [7,11,12];

%Initial conditions for the infested parcels
ind = [((inif(1)-1)*(N_stages)+8),((inif(2)-1)*N_stages+8),((inif(3)-1)*N_stages+8)];

x(ind) = 1000000; % adult mated females

x_hist=x; %We start the array x_hist to keep a historic of the evolution of the states

for t=1:(Simulation_time-1) %Loop for the simulation
    
    x = A_dis*x; %Compute the state at the next time step
    x_hist = [x_hist,x]; %Store the states
    
    %We recompute A based on the new temperature (for each area)
    for i=1:N_trees
        G_R(i) = growth_rate(Temp_areas(i,t+1),a,T_l,T_m,m); %Calling the growth rate function
        B_R(i) = birth_rate_suzuki(Temp_areas(i,t+1),alpha,gamma,lambda,delta,tau);%Calling the birth rate function
        M_R(i) = mortality_rate(Temp_areas(i,t+1),a1,b1,c1,d1,e1); %Calling the mortality rate function
        S_R(i) = s_r;
        R_remate(i) = r_remate;
        R_mate(i) = r_mate;   
    end
    
    Pest_stages = Initialize_stages_ode(B_R,M_R,G_R,S_R,R_mate,R_remate,Pest_stages,N_trees,N_stages); %We initialize the parameters associated to each stage class

    %Update the matrix A
    A_cont = compute_A_multi_16(N_trees,Pest_stages,N_stages, Adj_w ,w(t+1,:));
    sysc = ss(A_cont,[],eye(N_stages*N_trees),[]);
    sysd = c2d(sysc,1,'zoh');
    A_dis = sysd.A;
    
end


%% PLOTS
t1 = datetime(2018,4,1,12,0,0); %Simulations start on the 12 of April
t=t1+days(0:Simulation_time-1);

% Plot of every stage development

figure
plot(t,x_hist(6,:),t,x_hist(6+(N_stages),:),t,x_hist(6+(N_stages*2),:),t,x_hist(6+(N_stages*3),:),t,x_hist(6+(N_stages*4),:),...
t,x_hist(6+(N_stages*5),:),t,x_hist(6+(N_stages*6),:),t,x_hist(6+(N_stages*7),:),t,x_hist(6+(N_stages*8),:),...
t,x_hist(6+(N_stages*9),:),t,x_hist(6+(N_stages*10),:),t,x_hist(6+(N_stages*11),:),t,x_hist(6+(N_stages*12),:)...
,t,x_hist(6+(N_stages*13),:),t,x_hist(6+(N_stages*14),:),t,x_hist(6+(N_stages*15),:),'LineWidth',2);
legend("Area 1","Area 2","Area 3","Area 4","Area 5","Area 6","Area 7","Area 8","Area 9",...
"Area 10","Area 11","Area 12","Area 13","Area 14","Area 15","Area 16");
xlabel('time (days)');
ylabel('Number insects');





