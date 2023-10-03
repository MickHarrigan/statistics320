%% S22  CMPE320 Proj 2 Skeleton
close all;
clear;

PrA = 0.5;  % per the project
Ntrials = 500000; % make this as large as you can for your machine.  
%  From Project 1, more trials give results closer to the pdf.
A_minusA = (rand(1,Ntrials)<=PrA); % 1  = A, 0 = -A;
A_minusA = 2*(A_minusA-0.5);% convert to  +/-A;

Avalue = 2; % per assignment
sigma2 = 9/16; %per assignment;

N =  sqrt(sigma2)*randn(1,Ntrials); % zero mean variance = sigma2
R =  Avalue*A_minusA+N;  % R = (+/-A)+N;

tenSigma = sqrt(sigma2)*10;
dr=0.05;
rEdge=[-tenSigma-Avalue-dr/2:dr:tenSigma+Avalue]; % force bin center to zero

% Figure (1) is the scatterplot
% plots each output value of R for each Ntrials value
figure(1)
x = [1:Ntrials];
y = R;
plot(x,y,'b.'); %create the scatterplot use an appropriate x, an appropriate y, 
% the 'b.' will plot individual points in blue.

% Figure(2) is the histogram
% Now create the histogram, normalized to pdf, as in Project 1.
figure(2)
spdfR = histogram(R, 'BinEdges', rEdge, 'Normalization', 'pdf');

[Vr,Nbinr,r]=unpackHistogram(spdfR);  %I've provided a helper function to assist with histogram management

% Vr is values of the histogram bins
% Nbinr is number of bins
% r is the bin centers

edges = rEdge;
rGivenA = exp(-(edges-Avalue).^2/(2*sigma2))/sqrt(2*pi*sigma2);
rGivenNegA = exp(-(edges-(-Avalue)).^2/(2*sigma2))/sqrt(2*pi*sigma2);
fRr = rGivenA * 0.5 + rGivenNegA * 0.5; % put the equation for your fR(r) here

hold on;
plot(edges, fRr, 'LineWidth', 2); % plot your fRr
hold off;

% Make the plot look professional
xlabel('Voltage');
ylabel('Probability Density');
grid on;
legend('Random Variable R', 'Theoretical Value of R');
title('Probability Density of R');



figure(3); %Scatterplot for 2.2

%  Method 1
%  Notice the trick here.  (R>=0) will be 1 when true and 0 when false.
%  Multiplying point by point using .* will set all negative values to zero
%  and leave all postive values unchanged, thus creating the S for 2.1

S = (R>=0).*R; %  only accept R>=0;
x = [1:Ntrials];
y = S;
plot(x,y,'b.'); % scatterplot

% figure(5); % extra scatterplot
% plot(R, S, 'b.');

ds = dr;
sEdge = [-tenSigma:ds:tenSigma];

figure(4);

% Plot the CDF from the histogram and Theoretical CDF
% this can effectively be ignored in terms of the project document
CDF_S = histogram(S, 'BinEdges', sEdge, 'Normalization', 'cdf');
FSs = (sEdge>=0) .* ((1 - QQ((sEdge - Avalue)/sqrt(sigma2))) * 0.5 + 0.5); % Hint:  Use the QQ helper function to express the integral in terms of Q(s)

% differentiate the FSs for the pdf from the CDF function
SpdfDiff = diff(FSs);

hold on

plot(sEdge, FSs, 'LineWidth', 2); % plot FSs
hold off

%Make the plot look professional
grid on;
xlabel('Value of S');
ylabel('Cumulative Probability (CDF)');
title('CDF of S Given Perfect Diode Detector');
legend('S Normalized as CDF', 'Actual CDF Function');


figure(5);

% you may use subplots or not, as you desire.  If not, then you'll need new figures
spdfS = histogram(S, 'BinEdges', sEdge, 'Normalization', 'pdf'); %generate normalized histogram as in Project 1

[Values,Ns,s]=unpackHistogram(spdfS); %Use the helper function

i0 = min(find(s>=0)); % locate s nearest to zero
fSs = Values.*(s>0); % use the trick again 

% set the value to that of the middle value of R, which is the bin at zero
% times the bin width of that bin
PrS_is_0 = Vr(length(Vr)/2 + 1) * spdfR.BinWidth; % you will have some value


fSs(i0)=PrS_is_0/ds; % create an effective Dirac Delta function at zero
hold on
plot(s,fSs,'r','LineWidth',2); % plot in red on top of histogram
% this plots the pdf function as generated from the CDF above
plot(s+0.025, SpdfDiff/ds, 'LineWidth', 2); % concatenate a zero on the end as the last value is lost
hold off

% make the plot look professional
grid on;
xlabel('Value of S');
ylabel('Probability Density');
legend('Real Probability Density', 'Effective Probability Density');
   
% Would a rescaled version make things easier to see?


meanS =  mean(S); % sample mean from S array, not the histogram
varS = var(S);  % sample variance from S array, not the histogram
meanR = mean(R); % sample mean from R array, not the histogram
varR=var(R);  % sample variance from the R array.

% Print the results (example only, do what you want)
disp('--------');
disp('Section 2.1');
disp(['For the input, with A = +/-',num2str(Avalue),' variance ',num2str(sigma2),' and ',int2str(Ntrials),' trials,']);
disp(['the mean of R is ',num2str(meanR),' with variance ',num2str(varR),' =  ',num2str(Avalue^2),'  +  ',num2str(sigma2)]);
disp(['For Method 1 (ideal diode), the mean of S is ',num2str(meanS),' and the variance is ',num2str(varS)])
disp(['For Method 1 (ideal diode), diode(',num2str(meanR),') = ',num2str(meanR*(meanR>=0))]);

%================
% 2.2 uses abs as the function, but the same signal model.
%      We'll retain N and R, and just replace S

%New figures as necessary

S2 = abs(R); % the second method is absolute value.

figure(6);
% Plot the scatter plot
x = [1:Ntrials];
y = S2;
plot(x,y,'b.');

% New figure
figure(7);
s2Edge = sEdge;
spdfS2 = histogram(S2, 'BinEdges', s2Edge, 'Normalization', 'pdf');% generate the normalized histogram

[Values2,Ns2,s2]=unpackHistogram(spdfS2); % unpack for ease of use

% i0 = min(find(s>=0)); % locate s nearest to zero
% fSs = Values.*(s>0); % use the trick again 
% 
% % set the value to that of the middle value of R, which is the bin at zero
% % times the bin width of that bin
% PrS_is_0 = Vr(length(Vr)/2 + 1) * spdfR.BinWidth; % you will have some value

% fSs(i0)=PrS_is_0/ds; % create an effective Dirac Delta function at zero

fS2s = Values2; % create code for your described fSs
ds = dr;
hold on
% plot(s,fSs,'r','LineWidth',2); % plot in red on top of histogram
plot(s2, fS2s, 'r', 'LineWidth', 2); % your fSs
hold off

% Make your plot professional
grid on;
xlabel('');
ylabel('');
title('');
legend('');

% Compute and plot the CDF and print resultsmodifying lines 74-96 as necessary for this section
figure(8);

CDF_S2 = histogram(S2, 'BinEdges', s2Edge, 'Normalization', 'cdf');
FS2s = (1-QQ((s2Edge-Avalue) / sqrt(sigma2)));

hold on;
plot(s2Edge, FS2s, 'LineWidth',2);
hold off;

meanS2 =  mean(S2); % sample mean from S array, not the histogram
varS2 = var(S2);  % sample variance from S array, not the histogram

% Print the results (example only, do what you want)
disp('--------');
disp('Section 2.2');
disp(['For the input, with A = +/-',num2str(Avalue),' variance ',num2str(sigma2),' and ',int2str(Ntrials),' trials,']);
disp(['the mean of R is ',num2str(meanR),' with variance ',num2str(varR),' =  ',num2str(Avalue^2),'  +  ',num2str(sigma2)]);
disp(['For Method 2 (absolute value detector), the mean of S is ',num2str(meanS2),' and the variance is ',num2str(varS2)]);
disp(['For Method 2 (absolute value detector), diode(',num2str(meanR),') = ',num2str(meanR*(meanR>=0))]);

%================
% 2.3 uses S = R.^2 as the function, but the same signal model.
%      We'll retain N and R, and just replace S

% New plots
figure(9); % this is the scatterplot
S3 = R.^2;

% scatterplot
% this shows the bounds needed for sxEdge
x = [1:Ntrials];
y = S3;
plot(x,y,'b.');

s3Edge = [0:ds:30];
% you may use subplots or not, as you desire.  If not, then you'll need new figures
% this should be the 10th figure
figure(10);
spdfS3 = histogram(S3, 'BinEdges', s3Edge, 'Normalization', 'pdf'); %generate normalized histogram as in Project 1

[Values3,Ns3,s3]=unpackHistogram(spdfS3); %Use the helper function


% get the cdf from the histogram
figure(11);
CDF_S3 = histogram(S3, 'BinEdges', s3Edge, 'Normalization', 'cdf');

% get the histogram from the Q function (FS3s)
FS3s = (1 - QQ( (sqrt(s3Edge) - Avalue) / sqrt(sigma2) ));

% Do it all again, using Lines 54 to 96 as a model

hold on;
plot(s3Edge, FS3s, 'LineWidth',2);
hold off;

meanS3 =  mean(S3); % sample mean from S array, not the histogram
varS3 = var(S3);  % sample variance from S array, not the histogram

disp('--------');
disp('Section 2.3');
disp(['For the input, with A = +/-',num2str(Avalue),' variance ',num2str(sigma2),' and ',int2str(Ntrials),' trials,']);
disp(['the mean of R is ',num2str(meanR),' with variance ',num2str(varR),' =  ',num2str(Avalue^2),'  +  ',num2str(sigma2)]);
disp(['For Method 3 (squeare law detector), the mean of S is ',num2str(meanS3),' and the variance is ',num2str(varS3)]);
disp(['For Method 3 (square law detector), diode(',num2str(meanR),') = ',num2str(meanR*(meanR>=0))]);


% Print output table for use in report
%  This table provides the means and variances for the various options all
%  in one place.


% Jensens holds true due to the gER values being less than or equal to ES for all values
% func(E(x)) <= E(func(x)) <- this is the correct Jensens which is shown in the table through the collected data
disp('Output Table');
method =  [1:3]';
ES = [meanS,meanS2,meanS3]';
gER = [meanR*(meanR>=0), abs(meanR), meanR^2]'; 
table  =  [method ES gER];
%sprintf('%10.5f',table)
disp(table);