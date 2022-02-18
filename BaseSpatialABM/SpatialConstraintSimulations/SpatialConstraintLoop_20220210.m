
%%% Loop through spatial constraint simulations

clear
close all
rng(100)

size_vec = logspace(5,8,7);
bdg_vec = [1 2 5 20];

nsims = 48;

for size = size_vec
    for bdg = bdg_vec
        
        parfor n = 1:nsims
            
            out = SpatialConstraint_20220210(round(size,0),bdg);
            disp([size bdg n])
            csvwrite(strcat('SpatialConstraint_logsize',num2str(log10(size)),'_budgmax',num2str(bdg),'_iter',num2str(n)),out);
            
        end
        
    end
end