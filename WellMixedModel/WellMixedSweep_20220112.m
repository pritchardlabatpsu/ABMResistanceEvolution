% Sweep through recruitment rates for well-mixed model

clear
close all

adj0_vec = logspace(-2,-0.5,31);

npats = 50;

for i = 1:length(adj0_vec)
    adj_frac = adj0_vec(i);

    for pat = 1:npats
        output = WellMixedModel_20220112(adj_frac,pat);
        disp([adj_frac pat])
        writematrix(output,strcat('WellMixed_adjfrac',num2str(adj_frac),'_pat',num2str(pat),'_011222.csv'))
    end

end

