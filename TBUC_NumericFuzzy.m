clc;close all;clear all;
%% Loading the Datasets generated from 14 initial states
table=[];
datFiles = dir("*.dat"); 
for ND=1:length(datFiles)    
     Filename=datFiles(ND).name;
     tables{ND} = load(Filename);
     table=[table;(tables{ND})];
end
%% Step 1:- Divide the Input and Output spaces into fuzzy regions

% Generating Membership Functions for Input x
x=0:0.1:20 ; %Universe of discourse of Input x
MF_x(:,1)=trapmf(x,[0 0 1.5 6]);
MF_x(:,2)=trimf(x,[3 6 9]);
MF_x(:,3)=trimf(x,[8 9.5 11]);
MF_x(:,4)=trimf(x,[9 12 15]);
MF_x(:,5)=trapmf(x,[12 18.5 20 20]);
subplot(312);
plot(x,MF_x,'Linewidth',1.5);
xlabel('x');ylabel('\mu(x)');
title('Memebrship Function for Input x');
legend('S2','S1','CE','B1','B2');
ylim([0 1]);

% Generating Membership Functions for Input phi
phi = -115:0.1:295; % Universe of discourse of input phi.
MF_phi(:,1) = trimf(phi,[-115 -60 -9]);
MF_phi(:,2) = trimf(phi,[-19 26.5 40]);
MF_phi(:,3) = trimf(phi,[14 50 81]);
MF_phi(:,4) = trimf(phi,[76 88 97]);
MF_phi(:,5) = trimf(phi,[88 138.9 157]);
MF_phi(:,6) = trimf(phi,[130 174.32 230]);
MF_phi(:,7) = trimf(phi,[188 237.9 295]);
subplot(311);
plot(phi,MF_phi,'Linewidth',1.8);
xlabel('\phi');ylabel('\mu(\phi)');
title('Memebrship Functions for Input \phi');
legend('S3','S2','S1','CE','B1','B2','B3');
ylim([0 1]);

% Generating Membership Functions for output theta
theta=-40:0.1:40 ;%universe of discourse of output theta
MF_theta(:,1)=trimf(theta,[-40 -30.2 -20]);
MF_theta(:,2)=trimf(theta,[-33.2 -22.1 -8]);
MF_theta(:,3)=trimf(theta,[-13.8 -9 0.5]);
MF_theta(:,4)=trimf(theta,[-4.7 0.2 5]);
MF_theta(:,5)=trimf(theta,[0 8.2 16]);
MF_theta(:,6)=trimf(theta,[7.4 18.9 34]);
MF_theta(:,7)=trimf(theta,[20.5 39 40]);
subplot(313);
plot(theta,MF_theta,'Linewidth',2.4);
xlabel('\theta');ylabel('\mu(\theta)');
title('Memebrship Functions for Ouput \theta');
legend('S3','S2','S1','CE','B1','B2','B3');
ylim([0 1]);

%% Step 2: Generating Fuzzy Rules from Given Data Pairs

% Determining the degree for input x in different regions
mu_x_S2 = trapmf(table(:,1),[0 0 1.5 7]);
mu_x_S1 = trimf(table(:,1),[4 7 10]);
mu_x_CE = trimf(table(:,1),[9 10 11]);
mu_x_B1 = trimf(table(:,1),[10 13 16]);
mu_x_B2 = trapmf(table(:,1),[13 18.5 20 20]);
mu_x = [mu_x_S2 mu_x_S1 mu_x_CE mu_x_B1 mu_x_B2];

% Determining the degree for input phi in different regions
mu_phi_S3 = trimf(table(:,2),[-115 -65 -15]);
mu_phi_S2 = trimf(table(:,2),[-45 0 45]);
mu_phi_S1 = trimf(table(:,2),[15 52.5 90]);
mu_phi_CE = trimf(table(:,2),[80 90 100]);
mu_phi_B1 = trimf(table(:,2),[90 127.5 165]);
mu_phi_B2 = trimf(table(:,2),[135 180 225]);
mu_phi_B3 = trimf(table(:,2),[195 245 295]);
mu_phi = [mu_phi_S3 mu_phi_S2 mu_phi_S1 mu_phi_CE mu_phi_B1 mu_phi_B2 ...
    mu_phi_B3];
 
% Determining the degree for input theta in different regions
mu_theta_S3 = trimf(table(:,3),[-40 -40 -20]);
mu_theta_S2 = trimf(table(:,3),[-33 -20 -7]);
mu_theta_S1 = trimf(table(:,3),[-14 -7 0]);
mu_theta_CE = trimf(table(:,3),[-4 0 4]);
mu_theta_B1 = trimf(table(:,3),[0 7 14]);
mu_theta_B2 = trimf(table(:,3),[7 20 33]);
mu_theta_B3 = trimf(table(:,3),[20 40 40]);
mu_theta = [mu_theta_S3 mu_theta_S2 mu_theta_S1 mu_theta_CE mu_theta_B1 ...
    mu_theta_B2 mu_theta_B3];
 
% Calculates average degree of inputs and outputs and their corresponding
% row index
[avg_degree_x,avg_index_x] = max(mu_x,[],2);
[avg_degree_phi,avg_index_phi] = max(mu_phi,[],2);
[avg_degree_theta,avg_index_theta] = max(mu_theta,[],2);
% Defining Rules
Rules = [avg_index_x avg_index_phi avg_index_theta];
Rules_degree = avg_degree_x.*avg_degree_phi.*avg_degree_theta;

%% Step 3: Assign a Degree to Each Rule

R_temp=[Rules Rules_degree];
n=length(R_temp);
for i=1:n
    for j=1:n
        if((R_temp(i,1)==R_temp(j,1))&&(R_temp(i,2)==R_temp(j,2)))
            if(R_temp(i,4)>=R_temp(j,4))
            R_temp(j,3)=R_temp(i,3);
            R_temp(j,4)=R_temp(i,4);
            else
               R_temp(i,3)=R_temp(j,3);
               R_temp(i,4)=R_temp(j,4); 
            end
        end
    end
end
final_rules=unique(R_temp,'rows','sorted');%It removes redundent rows 
 % and sort the Data
%% Step 4 :-Create a Combined Fuzzy Rule base
% This step is followed connectors are used.
 
%% Step 5: Determine a Mapping Based on the Combined Fuzzy Rule Base.

figure;
axis([0 20 -10 100]); % Defining x-axis and y-axis limits
iter = 55; % number of iterations for trajectory plotting.
y_t = [-40 -20 -7 0 7 20 40]; % the points where centroid is lying .
b=4; %Length of the truck 
 %% Testing for different cases of Inputs
    x_input= [10 3 13]; % Input x (by user)
    phi_input= [220 -30 30]; % Input phi (by user)
for j=1:3
    y_sample = 2; % initial y position
    y_final = zeros(1,1);
    x_final = zeros(1,1);
    x_final(1,1) = x_input(j);
    y_final(1,1) = y_sample;
    Input = strcat('(',num2str(x_input(j)),',',num2str(phi_input(j)),')');
    for i = 1:iter % no. of iterations for trajectory tracking.

        % value of x.
        mu_xS1_test = trapmf(x_input(j),[0 0 1.5 7]);
        mu_xS2_test = trimf(x_input(j),[4 7 10]);
        mu_xCE_test = trimf(x_input(j),[9 10 11]);
        mu_xB1_test = trimf(x_input(j),[10 13 16]);
        mu_xB2_test = trapmf(x_input(j),[13 18.5 20 20]);

        % value of phi.
        mu_phiS3_test = trimf(phi_input(j),[-115 -65 -15]);
        mu_phiS2_test = trimf(phi_input(j),[-45 0 45]);
        mu_phiS1_test = trimf(phi_input(j),[15 52.5 90]);
        mu_phiCE_test = trimf(phi_input(j),[80 90 100]);
        mu_phiB1_test = trimf(phi_input(j),[90 127.5 165]);
        mu_phiB2_test = trimf(phi_input(j),[135 180 225]);
        mu_phiB3_test = trimf(phi_input(j),[195 245 295]);

        mu_x_test = [mu_xS1_test mu_xS2_test mu_xCE_test mu_xB1_test ...
            mu_xB2_test];
        mu_p_test = [mu_phiS3_test mu_phiS2_test mu_phiS1_test ...
            mu_phiCE_test mu_phiB1_test mu_phiB2_test mu_phiB3_test];
        mo = mu_x_test(final_rules(:,1)).*mu_p_test(final_rules(:,2)); 
        % product operation to determine the degree of output control.
        y_bar = y_t(final_rules(:,3));
        theta_output = sum(mo.*y_bar)/sum(mo); % Centroid defuzzification
        % Value of theta from the rule base.
        % Approximate kinematics of truck backer upper control 
        % for calculating the next states.
        x_input(j) = x_input(j) + cosd(phi_input(j) + theta_output) ... 
            + sind(theta_output)*sind(phi_input(j));
        x_final(i+1,1) = x_input(j);
        y_sample = y_sample + sind(phi_input(j) + theta_output) ...
          -sind(theta_output)*cosd(phi_input(j)); % displacement on y axis.
        y_final(i+1,1) = y_sample;
        phi_input(j) = phi_input(j) - asind(2*sind(theta_output)/b);

    end
    % Plot the trajectory.
    plot(x_final,y_final,'.-','MarkerSize',11.0);
    text(x_final(1,1),y_final(1,1)+1,Input);
    hold on;
end 
xlabel('x (in meters)');
ylabel('y (in meters)');
ylim([0 65]);
title('Truck trajectories using numerical-fuzzy controller');
legend('Trajectory for (10,220)','Trajectory for (3,-30)','Trajectory for (13,30)');