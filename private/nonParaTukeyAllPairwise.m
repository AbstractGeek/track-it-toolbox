function [stats] = nonParaTukeyAllPairwise(rankedData,tiedNums,alpha)
% [stats] = nonparaTukeyAllPairwise(rankedData,tiedNums,alpha)
%
%
%
load 'Nonparametric Tukey Test Critical Values.mat';
numOfInputs = length(rankedData);
N = NaN(numOfInputs,1);
R = NaN(numOfInputs,1);
tempsum = 0;
for i=1:numOfInputs
    N(i) = length(rankedData{i,1});
    R(i) = sum(rankedData{i,1});
end
% Total N
sumN = sum(N);
t = sum(tiedNums.^3 - tiedNums);
meanRanks = R./N;
meanRanksOrdered = meanRanks;
% Arrange mean ranks and preserve group number
stats.meanRanks = NaN(length(meanRanks),2);
[stats.meanRanks(:,1),stats.meanRanks(:,2)] = sort(meanRanks,'ascend');
meanRanks = stats.meanRanks;    
% Get Q value and get the number of comparisions
Qexp=nonParaTukeyCriticalValues(numOfInputs,nonParaTukeyCriticalValues(1,:)==alpha);          %#ok<NODEF>
numOfComp = (numOfInputs*(numOfInputs-1)/2);         % Number of comparisions

% Define stats - output variable, and make it a proper table
stats.table = cell(numOfComp+1,7);
stats.table{1,1} = 'Comparision (B vs. A)';
stats.table{1,2} = 'Difference (mean(RA) - mean(RB))';
stats.table{1,3} = 'SE';
stats.table{1,4} = 'Q';
stats.table{1,5} = 'Q (obtained from the table))';
stats.table{1,6} = 'Null Hypothesis - H0';
stats.table{1,7} = 'Conclusion';

% Loop variables
currentB = numOfInputs;     % Ranks
currentA = 1;

for i=1:numOfComp
    % Get necessary parameters for B
    RB = meanRanks(currentB,1);
    Bindex = meanRanks(currentB,2);
    nB = N(Bindex);
    % Get necessary parameters for A
    RA = meanRanks(currentA,1);
    Aindex = meanRanks(currentA,2);
    nA = N(Aindex);
    
    % Get SE,Q and h
    SE = sqrt(((sumN*(sumN+1)/12) - (t/(12*(sumN-1))))*(1/nA+1/nB));
    Q = (RB - RA)/SE;
    h = (Q>Qexp);
    
    if h
        conclusion = 'Reject H0';
    else
        conclusion = 'Accept H0';
    end
    
    stats.table{i+1,1} = strcat(int2str(Bindex),'vs. ',int2str(Aindex));
    stats.table{i+1,2} = RB-RA;
    stats.table{i+1,3} = SE;
    stats.table{i+1,4} = Q;
    stats.table{i+1,5} = Qexp;
    stats.table{i+1,6} = h;
    stats.table{i+1,7} = conclusion;
    
    % Propogate the loop
    currentA = currentA+1;
    if (currentA == currentB && currentB==1)
        % loop done
        disp('Something wrong, your loop is not ending properly! Check!!');
        break;
    elseif (currentA == currentB)
        currentB = currentB-1;
        currentA = 1;
    end    
    
end

end