close all;
clear;
% Section 3.3 exponentially distributed
disp(' ');
disp('Section 3.3 exponentially distributed');

% Set up new plots as necessary.  Remember, you need ALL of the plots


Ntrials = [10 1000 100000]; %set according to the assignment
lambda = 0.5;  % as given in the assignment

theory_m =  1/lambda;   % analytical or population mean == mu
theory_v = 1/(lambda^2); % analytical or population variance == sigma^2
for ktrials = 1:length(Ntrials)   % repeat for each number of trials
    
      figure();
      % Set the bin edges
      bins_edges = [0.5: 20.5];
      
      % Generate the appropriate number of independent random trials
      data = randx(1,Ntrials(ktrials),lambda); % randx is given 
      
      % Compute (and plot) the raw histogram
      hist_raw = histogram(data, 'BinEdges', bins_edges);
      
      % Compute (and plot) the normalized histogram
      %    Use 'Normalization','pdf' because an exponential RV is a
      %    continuous random variable and has a pdf, not a pmf.
      hist_norm = histogram(data, 'BinEdges', bins_edges, 'Normalization', 'pdf');
      
      % Compute and plot the true pdf on the same axes as the normalized
      % histogram
      k = (0:20);
      theoretical_pdf = lambda * exp(-lambda * k);
      hold on;
      plot(k, theoretical_pdf);
      hold off;

      % Compute the sample mean and variance and display them
      sample_mean = mean(data);
      sample_var = var(data);

      % labels and titles
      xlabel('Input');
      ylabel('Probability');
      title(['Section 3.3: N = ', int2str(Ntrials(ktrials)), ' Probability Density Function']);
      legend('Scaled Histogram', 'Theoretical Prob Density Fnc (PDF)');

      % Either compare the sample mean and variance with the population
      % mean and variance, or save the results to output later
      disp(['For N = ', int2str(Ntrials(ktrials)), ', sample and population mean are ', num2str(sample_mean), ...
          ' & ', num2str(theory_m), ...
          ' with sample and population variance ', num2str(sample_var), ...
          ' & ', num2str(theory_v)]);
end;

%Compare the sample mean and variance to the population mean and variance
%as requested

disp('-----------');
disp(' ');