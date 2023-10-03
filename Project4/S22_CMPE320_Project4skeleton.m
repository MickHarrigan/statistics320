% S22 CMPE320 Project 4 BASK simulation
%
close all
clear

disp('CMPE320 Spring 2022 Project 4:  BASK');

% 2.1
p0=[0.01:.01:0.99]; % range of p0 given in 2.1
gamma_dB = 10; % given in 2.1
A=1; % given in 2.1
sigma2 = 10^(-gamma_dB/10); % convert dB to power, then to sigma^2, see 2.1
tauMAP = (sigma2 * log((1 - p0) / p0)) / (2 * A); % put the functional expression for your tauMAP derived in 2.1
figure; % new figure
% p0 and tauMAP, per instructions in 2.1
%
% Professional labels, etc.
plot(p0, tauMAP);


%2.2

figure; %new figure
p0=0.8; % per 2.2
gamma_dB=10; % per 2.2
sigma2=  % convert gamma_dB to sigma^2 as before
Ntrials = ; % set a large number of trials (500,000?)
A = 1; % per 2.1
B = ; % create random variable B based on p0, as in Project 2
M = ; % "modulate" by converting B per equation equation (1.1) in writeup
N = ; % generate noise from N(0,sigma2) as in Project 2
R=M+N; % per signal model for Project 4
Rhist = ; % generate appropriately normalized histogram for R
r = [-2:.01:2]; % fine grid for plotting fRr
fRr = ; % expression for the fRr you derived.  This is very similar to Project 2, but explicitly includes p0

% plot the smooth fRr over the histogram, as in previous projects
% make the plot professional
    

% Prepare for simulations of 2.3.3 and 2.4.3
% This does the simulations in a loop on values of p0

%Set Parameters

gamma_dB=[0:0.5:10 10.25:.25:14]; % A^2/sigma^2 in dB
gamma = 10.^(gamma_dB/10);  % as power ratio

A=1; % assume unity amplitude, per project description

Nbits =;% lots and lots of bits, 5,000,000 if your computer can handle it
p0=[0.5, 0.8]; % values of p0 for ML and MAP

thresholds=zeros(length(p0),length(gamma)); %thresholds vary with sigma2 and p0
results = zeros(length(p0),length(gamma));
pBT = zeros(length(p0),length(gamma));

for kP0=1:length(p0)
    
    
% START THE SIMULATION

    b= ; % create the bits with appropriate probability 0 with p0, 1 with 1-p0
    m= ; % convert to +/-A (A=1) per equation 1.1, remeber 0->+A, 1->-A
    
%Loop on SNR
    for kSNR=1:length(gamma)
    
    % One long continuous string of coded bits and add noise
    % This step both modulates and adds AWGN
           sigma2 = ;  % use the energy ratio in gamma, not gamma_dB
           sigma =  ;  % compute sigma from variance, sigma2
           n = ; % compute noise values from N(0,sigma2);
           r = ; % create the received signal 
           
           %Compute the thresholds as a function of sigma2, A, and p0, per
           %your derivation
           thresholds(kP0,kSNR) = ;
           
           bkhat = (r <= thresholds(kP0,kSNR)); % 1 if less than threshold, 0 if greater
           errors = mod(bkhat - b,2); % 1 = error, 0 = no error
           results(kP0,kSNR) = sum(errors)/Nbits; % pBX for this SNR and this p

           
           pBT_given0 = ;  % per your derivation
           pBT_given1 = ;  % per your derivation
           pBT(kP0,kSNR) = % Law of Total Probability;
         
    end; %loop on SNR
    
    %Plot the results for this P0 and all SNRs
    
        figure; % set new figure
        h=semilogy(gamma_dB,pBT(kP0,:),'k-',gamma_dB,results(kP0,:),'ro');
        set(h,'LineWidth',1.5);
        
        % Make your plots professional with labels, scales, grids, legends
        % etc.
        
 
    end;% loop on P0
    
    
    % Section 2.5
    
    % Combined plot
    figure; % set new figure
    % Plot both of the theoretical curves on the same axes
        h=semilogy(gamma_dB,pBT(1,:),'k-',gamma_dB,pBT(2,:),'r--'); % plot both theoretical on same scale
        set(h,'LineWidth',1.5);
  
    % Make the plot professional

    % Now plot the ratio of the theoretical errors, with p0=0.5 in
    % denominator
    figure; % another new figure
    % compute the ratio of the MAP probability of error (pBT(2,:)) and the
    % ML probability of error (pBT(1,:)) and plot with gamma_dB on the
    % x-axis
        
        
    % Make the plot professional

        