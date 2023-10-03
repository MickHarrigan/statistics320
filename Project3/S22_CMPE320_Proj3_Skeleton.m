%S21 CMPE320 Project3 Solns
% Modified for S21, EFCL 4/15/2021

close all;
clear;  

Ntrials = 100000;
kplot=0; 

% This next line is what MATLAB calls an anonymous function, which is a
% function that we use only in one script. the function call is
% (x,mu,sig2), and the returned value is a number.  sig2 = sigma^2

fgauss = @(x,mu,sig2) exp(-0.5*((x-mu).^2/sig2))/sqrt(2*pi*sig2);


%% Problem 2.1  Uniformly distributed
disp(['---------------']);
disp(['Section 2.1']);
Nsum = [2,6,12]; % checked  with S21 assignment
for k=1:length(Nsum)
    
    xd = rand(Nsum(k),Ntrials);  % generate [Nsum(k) by Ntrials] array of random values
    xs = sum(xd); % "sum" adds down the MATLAB columns, thereby giving us the sum of Nsum(k) values
    xmin = 0;
    xmax = Nsum(k); % largest value of the sum
    mu = Nsum(k)*0.5; % 0.5 is the E[X]
    sig2 = Nsum(k)*(1/12); % 1/12 is var[X]
    m = mean(xs); % sample mean
    S = var(xs); % sample variance
    disp(['For ',int2str(Ntrials),' independent trials of the sum of ',int2str(Nsum(k)),' iid rv from U(0,1)']);
    disp(['   the theoretical mean is ',num2str(mu),' and the sample mean is ',num2str(m)]);
    disp(['   the theoretical variance is ',num2str(sig2),' and the sample variance is ',num2str(S)]);
    
    dx=0.1; % fine grain dx for plotting
    x = [0:dx:Nsum(k)+1]; % fine grain for fY(y)
    
    % create new figure...
    figure();
    %...and then a new scaled histogram using the values of xs
    scaledHist = histogram(xs, 'BinEdges', x, 'Normalization', 'pdf');
    
    % And unpack the data using unpackHistogram from Project 2
    [vals, numBins, binCenters] = unpackHistogram(scaledHist);

    % And plot the Gaussian pdf on top of the histogram with labels and
    hold on;
    plot(x, fgauss(x,mu,sig2), 'LineWidth', 2);
    hold off;
    % grids and the other elements of a professional plot
    xlabel('Value of Y');
    ylabel('Probability Density of Y');
    grid on;
    legend('Random Variable Y', 'Theoretical Value of Y');
    title(['Probability Density of Y for N = ', num2str(Nsum(k)), ' in ', num2str(Ntrials), ' trials']);
    
    
end;
% there are length(Nsum) plots to this point.

%% Problem 2.2  Uniformly distributed discrete
disp('---------------');
disp('Section 2.2');
Nsum = [2,20,40]; % checked with S22 assignment
Nsides = 8; %   checked with S22 assignment
for k = 1:length(Nsum)
    xd = randi(Nsides,Nsum(k),Ntrials); % generates Ntrials values between 1 and Nsides
    xs = sum(xd); % sum of xd values
    % 4.5 is the E[8 sided die], then times amount of dice (Nsum(k))
    mu = Nsum(k) * 4.5; % Mean of uniform distribution
    % E[X^2] = sum( (1..8)^2 / 8 ) = 25.5
    % sample var = 25.5 - 4.5^2 = 5.25, E[X^2] - (E[X])^2 = 5.25
    sig2 = Nsum(k) * 5.25; % Variance of uniform distribution
    m = mean(xs);
    S = var(xs);
    disp(['For ',int2str(Ntrials),' independent trials of the sum of ',int2str(Nsum(k)),' iid rv from U(0,1)']);
    disp(['   the theoretical mean is ',num2str(mu),' and the sample mean is ',num2str(m)]);
    disp(['   the theoretical variance is ',num2str(sig2),' and the sample variance is ',num2str(S)]);

    x = [-0.5:Nsum(k) * Nsides + 0.5];
    % new figure made
    figure();
    % scaled histogram made
    scaledHist = histogram(xs, 'BinEdges', x, 'Normalization', 'pdf');
    
    % unpack values from histogram
    [vals, numBins, binCenters] = unpackHistogram(scaledHist);
    
    % And plot the Gaussian pdf on top of the histogram with labels and
    hold on;
    plot(x, fgauss(x,mu,sig2), 'LineWidth', 2);
    hold off;
    % grids and the other elements of a professional plot
    xlabel('Value of Y');
    ylabel('Probability Density of Y');
    grid on;
    legend('Random Variable Y', 'Theoretical Value of Y');
    title(['Probability Density of Y for N = ', num2str(Nsum(k)), ' in ', num2str(Ntrials), ' trials']);

end

% Do the experiment again for a large number of trials (Ntrials) and the
% specified number of terms in the sum (Nsum)
    
%Problem 2.3  Exponentially distributed 
disp('---------------');
disp('Section 2.3');
Nsum = [5,50,150]; % checked with S22 assignment
lambda=0.5;
for k = 1:length(Nsum)
    xd = randx(Nsum(k), Ntrials, lambda);
    xs = sum(xd);
    mu = Nsum(k) / lambda;
    sig2 = Nsum(k) / (lambda^2);
    m = mean(xs);
    S = var(xs);
    disp(['For ',int2str(Ntrials),' independent trials of the sum of ',int2str(Nsum(k)),' iid rv from U(0,1)']);
    disp(['   the theoretical mean is ',num2str(mu),' and the sample mean is ',num2str(m)]);
    disp(['   the theoretical variance is ',num2str(sig2),' and the sample variance is ',num2str(S)]);
    if (k == 1)
        x = [-0.5:1:Nsum(k) * 8];
    else
        x = [-0.5:1:Nsum(k) * 4];
    end
    figure();

    scaledHist = histogram(xs, 'BinEdges', x, 'Normalization', 'pdf');
    % unpack values from histogram
    [vals, numBins, binCenters] = unpackHistogram(scaledHist);

    hold on;
    plot(x, fgauss(x,mu,sig2), 'LineWidth', 2);
    hold off;
    % grids and the other elements of a professional plot
    xlabel('Value of Y');
    ylabel('Probability Density of Y');
    grid on;
    legend('Random Variable Y', 'Theoretical Value of Y');
    title(['Probability Density of Y for N = ', num2str(Nsum(k)), ' in ', num2str(Ntrials), ' trials']);
end
%And again, with samples drawn from randx provided with Project 1.


%Problem 2.4  Sum of iid Bernoulli trials
disp('---------------');
disp('Section 2.4');
Nsum = [4,8,40];  % Checked with Project 2  40 is, in this case, a big number
p = 0.5;
for k = 1:length(Nsum)
    xd = rand(Nsum(k), Ntrials)<= p; % generates numbers between 0 and 1
    xs = sum(xd);
    mu = Nsum(k) * p; % value of p, probability of 1 appearing
    sig2 = Nsum(k) * p * (1 - p); % variance of bernoulli
    m = mean(xs);
    S = var(xs);
    disp(['For ',int2str(Ntrials),' independent trials of the sum of ',int2str(Nsum(k)),' iid rv from U(0,1)']);
    disp(['   the theoretical mean is ',num2str(mu),' and the sample mean is ',num2str(m)]);
    disp(['   the theoretical variance is ',num2str(sig2),' and the sample variance is ',num2str(S)]);

    figure();
    subplot(2,1,1);
    % subplot 1
    x = [-0.5:Nsum(k) + 0.5];
    scaledHist = histogram(xs, 'BinEdges', x, 'Normalization' , 'probability');
    hold on;
    % the 0:Nsum(k) is equivalent to x (used as k for nchoosek) without having a negative value
    nchoosek = factorial(Nsum(k)) ./ (factorial(Nsum(k) - (0:Nsum(k))) .* factorial(0:Nsum(k)));
    pmfBern =  nchoosek .* p .^ (0:Nsum(k)) .* (1-p) .^ (Nsum(k) - (0:Nsum(k)));
    stem((0:Nsum(k)), pmfBern, 'LineWidth', 2);
    hold off;
    
    % grids and the other elements of a professional plot
    xlabel('Value of Y');
    ylabel('Probability of Y');
    grid on;
    legend('Random Variable Y', 'Theoretical Value of Y');
    title(['Probability of Y for N = ', num2str(Nsum(k)), ' in ', num2str(Ntrials), ' trials']);

    % subplot 2
    subplot(2,1,2);
    dx = 1;
    x = [-0.5:dx:Nsum(k) + 0.5];
    scaledHist = histogram(xs, 'BinEdges', x, 'Normalization', 'pdf');
    % unpack values from histogram
    [vals, numBins, binCenters] = unpackHistogram(scaledHist);

    hold on;
    plot(x, fgauss(x,mu,sig2), 'LineWidth', 2);
    hold off;
    % grids and the other elements of a professional plot
    xlabel('Value of Y');
    ylabel('Probability Density of Y');
    grid on;
    legend('Random Variable Y', 'Theoretical Value of Y');
    title(['Probability Density of Y for N = ', num2str(Nsum(k)), ' in ', num2str(Ntrials), ' trials']);
end
%Don't forget to plot both the exact (PMF) answer and the approximate
%Central Limit Theorem answer.  Use good practice (i.e., stem(x,y) for the
%PMF, but use plot(x,y) for the pdf.



    
    
    
   