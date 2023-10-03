function [Values,NumBins,binCenters]=unpackHistogram(hgram);
%function [Values,NumBins,binCenters]=unpackHistogram(hgram)
%  This function unpacks the data structure provided by the MATLAB function
%  histogram into simple arrays that can be used for further processing.
%  By using this function, the user does not need to understand the
%  underlying structure of the histogram object, and, in fact, doesn't need
%  to know anything about structures at all.
%  
%  Calling Parameters
%      hgram:  a histogram structure created by calling the MATLAB function
%      HISTOGRAM with any of its scaling options.  
%
%  Returned parameters
%      Values:  a one dimensional array containing the value in each bin of
%      the histogram.  If the histogram is unscaled, then this is the
%      number of counts that occurred in each bin.
%
%      NumBins:  number of bins in the histogram, as determined by the
%      MATLAB function HISTOGRAM
%      
%      binCenters:  center values of the bins used by HISTOGRAM
%
%  Usage:  An unscaled histogram consisting of the count of the number of
%  data samples in each bin is created by 
%         
%        hgram = histogram(data)
%
%  The data may then be unpacked using this function
%  
%        [Values, NumBins, binCenters] = unpackHistogram(hgram);
%
%  This function DOES NOT WORK with the "old" MATLAB hist() function.
%
%  EFCL 2/27/2021
%

Values = hgram.Values;
NumBins = hgram.NumBins;
BinEdges = hgram.BinEdges;
BinWidth = hgram.BinWidth;
binCenters = BinEdges(1:NumBins)+BinWidth/2;
