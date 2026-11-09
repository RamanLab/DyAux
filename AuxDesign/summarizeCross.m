function T = summarizeCross(crossWT, crossMut, names, minActive)
% crossWT, crossMut : metabolite x time matrices (CommunityResults*.CrossMets)
%                     + = species 1 -> 2, - = species 2 -> 1, 0 = no cross-feeding
% names             : cell array of exchange names (CommExchangeMets)
% minActive         : min fraction of time steps a metabolite must be cross-fed
%                     to count as sustained (default 0.5)
if nargin < 4, minActive = 0.1; end

% per-metabolite summary of one community
amt   = @(C) sum(abs(C), 2);   % total amount exchanged over the run
act   = @(C) mean(C ~= 0, 2);  % fraction of time steps with cross-feeding
n12   = @(C) sum(C > 0, 2);    % steps where species 1 feeds species 2
n21   = @(C) sum(C < 0, 2);    % steps where species 2 feeds species 1

amtWT  = amt(crossWT);         amtMut  = amt(crossMut);
actMut = act(crossMut);
s12    = n12(crossMut);        s21     = n21(crossMut);

% dominant direction in the mutant
dir = repmat({'none'}, numel(names), 1);
dir(s12 > s21) = {'1->2'};
dir(s21 > s12) = {'2->1'};
dir(s12 > 0 & s12 == s21) = {'mixed'};

gain   = amtMut - amtWT;                  % extra exchange in mutant vs WT
T = table(names(:), amtWT, amtMut, gain, actMut, dir, ...
    'VariableNames', {'met','amtWT','amtMut','gain','fracActive', ...
                      'direction'});
T = sortrows(T, 'gain', 'descend');       % biggest mutant-specific gains first
T = T(T.amtMut > 0 | T.amtWT > 0, :);     % drop metabolites never cross-fed
end
