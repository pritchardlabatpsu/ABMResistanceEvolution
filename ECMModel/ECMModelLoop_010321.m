%%% Loop through CAF Spatial Model

% updates for 5/11: look at degree of genetic resistance?
% updates for 5/18: added de nov mutations; allow tumor to grow and be
% treated

clear
close all

caf_vec = [0 0.01 0.05 0.1 0.25];
bdg_vec = [1 2 5 20];
mgrt_vec = [0 0.5 1];
size_vec = [1e4];
nsim = 48;

rng(100)

for k = 1:length(size_vec)
    size0 = size_vec(k);

    for i = 1:length(caf_vec)
        caf_frac = caf_vec(i);

        for j = 1:length(bdg_vec)
            bdgmax = bdg_vec(j);
            
            for k = 1:length(mgrt_vec)
                mgrt = mgrt_vec(k);
        
                parfor iter = 1:nsim
                    out = ECMModel_010321(caf_frac,bdgmax,size0,iter);
                    csvwrite(strcat('CAFModel_caffrac',num2str(caf_frac),'_bdgmax',num2str(bdgmax),'_mgrt',num2str(mgrt),'_cafeffrad',num2str(m),'_logsize',num2str(log10(size0)),'_iter',num2str(iter,'%03.f'),'_122721.csv'),out);
                    disp(iter);
                end
                
            end
            
        end
        
    end
end

