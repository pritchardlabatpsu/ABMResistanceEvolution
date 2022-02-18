%%% Spatial Agent-Based Model of Spatial Resistance in Spheroids
%%% Scott Leighow - 10/21/19

% clear
% close all
% rng(100)

% caf_frac = 0.05;
% size0 = 1e4;
% bdgmax = 1;
% m = 3; % radius of CAF effect
% iter = 1;

function output = ParacrineModel_122721(caf_frac,bdgmax,m,size0,tmax,iter)

%% Parameters

% Order of parameters: [sen res adj]

% Rates - using BT20 numbers at 1 uM letrozole; see DeathRates_041320.xlsx
% units in [/hr]
ngr = [0.0238 0.0238 0.0238];   % net growth rate
d = [0.001 0.001 0.001];        % tonic death rate
b = ngr+d;                      % division rate
a = [0.0410 0 0.0165];          % drug kill rate
                                % assume absolute genetic resistance
                                
% Population parameters
nsen0 = size0; % number of sensitive cells
nres0 = 5;  % number of resistant cells
ncaf0 = round(size0*caf_frac);  % number of CAFs (variable relative to total cells)

pop0 = [nsen0 nres0 0];
n0 = sum(pop0);

% Define neighborhood matrix
nghbrs1 = combvec(-1:1,-1:1,-1:1)';
% [~,row0s] = ismember([0 0 0],nghbrs1,'rows');
% nghbrs1(row0s,:) = [];

% Define matrix of positions within m units of a central CAF
% Note that m need not be a whole number; e.g. if m=1.5 then identify all
% cells whose centroids are <=1.5 units from CAF by Euclidean distance
m_ceil = ceil(m);
nghbrsm = combvec(-m_ceil:m_ceil,-m_ceil:m_ceil,-m_ceil:m_ceil)';
rm_rows = vecnorm(nghbrsm')'>m;
nghbrsm(rm_rows,:) = [];

% Simulation parameters
% nsim = 30;   % number of simulations
nmax = 4*n0; % maximum number of cells

% bdgmax = 2; % max budging distance during cell division

% output = NaN(nsim,4);

% Pretty colors
clr_sen = [22 193 180]/255;
clr_res = [253 13 57]/255;
clr_adj = [97 85 191]/255;

%% Simulation

% for i = 1:nsim
    
%%% Initialize state

% Randomly seed cells
len = ceil(n0^(1/3)); % determine length of initial cube
len = len + ceil(log10(size0)) + 1; % Add additional space based on size0 to make geometrically realistic
if rem(len,2)==0       % if length is even, increase by one
    len = len+1;
end
lenvec = -(len-1)/2:(len-1)/2;
cube0 = combvec(lenvec,lenvec,lenvec)';    
unocc_idx = (1:size(cube0,1))';

euc0 = sqrt(sum(cube0.^2,2)); % calculate Euclidean distance of potential points
pos0 = [0 0 0];
while size(pos0,1) < n0 % seed cells (preferentially close to center) until all cells are assigned position
    unocc_idx = (1:size(cube0,1))';
    unocc_idx(ismember(cube0,pos0,'rows')) = [];
    cells_idxi = randsample(unocc_idx,n0-size(pos0,1),true,exp(-5*euc0(unocc_idx)));
    pos_i = cube0(cells_idxi,:);
    pos0 = [pos0; pos_i];
    pos0 = unique(pos0,'rows');
end

pos0_caf = pos0(randsample(1:n0,ncaf0,false),:); % randomly choose points of spatial resistance

% Randomly assign cell phenotype
phn0 = [ones(pop0(1),1); repmat(2,pop0(2),1); repmat(3,pop0(3),1)];
phn0 = phn0(randperm(length(phn0)));

% Assign sen cells neighboring to caf cells as adj
n_adj = size(nghbrsm,1); % number of cells receiving benefit per CAF
pos0_adj = NaN(ncaf0*n_adj,3);
for k = 1:size(pos0_caf,1)
    pos0_cafk = pos0_caf(k,:);
    pos0_adjk = pos0_cafk+nghbrsm;
    pos0_adj((n_adj*(k-1)+1):(n_adj*k),:)=pos0_adjk;
end

pos0_mskd = pos0;
pos0_mskd(phn0==2,:) = NaN; % mask res cells

idx_adj = ismember(pos0_mskd,pos0_adj,'rows');
phn0(idx_adj) = 3; % assign sen cells next to caf cells as adj

state0 = [pos0 phn0];

%     % Visualize initial spheroid
%     clrs = NaN(n0,3);
%     clrs(phn0==1,:) = repmat(clr_sen,sum(phn0==1),1);
%     clrs(phn0==2,:) = repmat(clr_res,sum(phn0==2),1);
%     clrs(phn0==3,:) = repmat(clr_adj,sum(phn0==3),1);
%     scatter3(pos0(:,1),pos0(:,2),pos0(:,3),5e2,clrs,'.')
%     hold on
%     scatter3(pos0_caf(:,1),pos0_caf(:,2),pos0_caf(:,3),'k*')
%     hold off

%%% Initialize
state_curr = state0;
state_next = state_curr;
phn_curr = phn0;

state = NaN(nmax,4);
state(1:size(state_curr,1),:) = state_curr;
s = 1;

j = 1;
t = 0;
t_curr = sum(t);
% tmax = 500;

output = [j t_curr pop0];
if size0 >= 1e3
    str_intrvl = round(size0/10);   % store every x iterations in output matrix
else
    str_intrvl = 2;
end

%%% Simulate events
%     while size(state_curr,1) > 1 && size(state_curr,1) < nmax && sum(t) < 30
% while sum(phn_curr==2)>0 && sum(phn_curr==2)<size0 && size(state_curr,1) < nmax
while sum(phn_curr==2)>0 && size(state_curr,1)<nmax && t_curr<tmax

    state_curr = state_next;
    cellrows = ~isnan(state_curr(:,1));
    n_curr = sum(cellrows);
    state_curr = state_curr(cellrows,:);
    pos_curr = state_curr(cellrows,1:3);
    phn_curr = state_curr(cellrows,4);
    t_curr = sum(t);

    nS = sum(phn_curr==1);
    nR = sum(phn_curr==2);
    nA = sum(phn_curr==3);
    pop_curr = [nS nR nA];

    evts = [(b.*pop_curr)';    % cell division event
            (d.*pop_curr)';    % natural death event
            (a.*pop_curr)'];   % drug kill event

    theta = sum(evts);
    t(j) = -log(rand)/theta;

    idx_evt = randsample(length(evts),1,true,evts);

    if idx_evt<=3       % cell division

        % Draw cell
        phn_i = idx_evt; % identify phenotype of cell
        cell_i = randsample(n_curr,1,true,phn_curr==phn_i);
        pos_i = pos_curr(cell_i,:);

        % Define hollow 'cubic sphere' with radius k around cell i and scan for
        % free spaces (pos_f)
        pos_f = NaN;
        k = 1;
        prevscan = [0 0 0];

        while isnan(sum(pos_f)) && k<=bdgmax

            % Begin with cube of length 2k+1
            nghbrd_k = combvec(-k:k,-k:k,-k:k)';

            % Remove previously scanned positions and those more than k
            % units from position i
            nghbrd_k(ismember(nghbrd_k,prevscan,'rows'),:) = [];

            euc_k = sqrt(sum(nghbrd_k.^2,2));
            nghbrd_k(euc_k>k,:) = [];

            % Scan for unoccupied positions k units from position i
            nghbrs_ik = pos_i+nghbrd_k;
            free_k = nghbrs_ik(~ismember(nghbrs_ik,pos_curr,'rows'),:);
            free_ik = free_k - pos_i;

            if isempty(free_ik)
                prevscan = [prevscan; nghbrd_k];
                k = k+1;                    
            else
                % Choose free position f, favoring those nearest center
                % of spheroid
%                     pos_cntr = mean(pos_curr);
%                     euc_cntrk = sqrt(sum((free_k-pos_cntr).^2,2));
%                     pos_f = free_k(randsample(size(free_k,1),1,true,exp(-5*euc_cntrk)),:);

                % Choose free position f closest to dividing cell i
                norms_ik = sqrt(sum(free_ik.^2,2));
                min_free_k = find(norms_ik==min(norms_ik));
                if length(min_free_k)==1
                    pos_f = free_k(min_free_k,:);
                else % Given the option, choose position f to be candidate space closest to center
                    pos_cntr = mean(pos_curr);
                    nrst_k = free_k(min_free_k,:);
                    dst_cntr_nrstk = sqrt(sum((nrst_k-pos_cntr).^2,2));
                    [~,min_dcn_idx] = min(dst_cntr_nrstk);
                    pos_f = nrst_k(min_dcn_idx,:);
                end                    
            end

        end

        % If unable to find free position within bdgmax distance,
        % continue to next iteration
        if isnan(pos_f)
            continue
        end

        % Push cells between i and f until free space adjacent to i
        % Definitions:
        % position i = site of dividing cell
        % position f = site of initially free space
        % position g = site of cell that moves to position f
        % position j = site of new cell (adjacent to position i)

        nghbrsi1 = pos_i + nghbrs1;
        nghbrsi1(ismember(pos_i,nghbrsi1,'rows'),:) = [];

        while ~ismember(pos_f,nghbrsi1,'rows')

            % Find position g adjacent to f closest to line through i
            % and f

            vec_if = pos_f - pos_i;
            dstr = NaN(26,1);
            for r = 1:26
                % Candidate for migrant cell g
                pos_gr = pos_f+nghbrs1(r,:);
                vec_ig = pos_gr-pos_i;

                % Consider only positions g that are closer to i than f
                % is to i
                if norm(vec_ig) < norm(vec_if)
                    % Use parallelogram definition of cross product to
                    % find distance between candidate position g and
                    % line passing through i and f
                    dstr(r) = norm(cross(vec_ig,vec_if))/norm(vec_if);
                end
            end

            % Choose position g
            idx_mindst = find(dstr==min(dstr));
            if length(idx_mindst)>1
                pos_g = pos_f+nghbrs1(randsample(idx_mindst,1),:);
            else
                pos_g = pos_f+nghbrs1(idx_mindst,:);
            end

            % Move cell in position g to position f
            pos_curr(ismember(pos_curr,pos_g,'rows'),:) = pos_f;

            % Position g is now open to receive next pushed cell
            pos_f = pos_g;

        end

        % Once a free space is adjacent to cell i, it divides ==> cell j
        pos_j = pos_f;
        pos_curr = [pos_curr; pos_j];

        % Identify phenotype of cells after budging
        if phn_i==2 % if cell_i is res, so is cell_j
            phn_curr = [phn_curr; 2];
        else
            phn_curr = [phn_curr; 0];
        end

        % Identify phenotype of remaining cells after budging
        pos_mskd = pos_curr;
        pos_mskd(phn_curr==2,:) = NaN; % mask res cells
        phn_curr(phn_curr~=2) = 1; % initially assign all nonres cells as sen
        idx_adj = ismember(pos_mskd,pos0_adj,'rows');
        phn_curr(idx_adj) = 3; % assign sen cells next to caf cells as adj

        state_curr = [pos_curr phn_curr];

        j = j+1;
        state_next = NaN(nmax,4);
        state_next(1:length(phn_curr),:) = state_curr;
%             disp(length(phn_curr));

    elseif idx_evt>3    % cell death

        % Draw cell
        if idx_evt<=6
            phn_i = idx_evt-3;
        else
            phn_i = idx_evt-6;
        end
        cell_i = randsample(n_curr,1,true,phn_curr==phn_i);

        % Update state
        pos_curr(cell_i,:) = [];
        phn_curr(cell_i) = [];

        state_curr = [pos_curr phn_curr];

        j = j+1;
        state_next = NaN(nmax,4);
        state_next(1:length(phn_curr),:) = state_curr;
%             disp(length(phn_curr));

    end
    
    if rem(j,str_intrvl)==1
        store_idx = floor(j/str_intrvl)+1;
        output(store_idx,:) = [j t_curr pop_curr];
        
        % Update complete state: comment for efficiency; uncomment if plotting
%         state(1:size(state_curr,1),:,s) = state_curr;
%         s = s+1;
        
        disp([iter t_curr])
    end

end

output = [output;
          j t_curr pop_curr];
      
% save('ExampleSimulation_size1e4_bdgmax30_091220.mat')
% 
% if sum(phn_curr==2)==0
%     t_100 = Inf;
% elseif size(state_curr,1) >= nmax
%     t_100 = NaN;
% else
%     t_100 = sum(t);
% end
% 
% %     disp(i)
% % output(i,:) = [bdgmax ncaf0 i t_100];
% 
% time = cumsum([0 t]);
% 
% phn = state(:,4,:);
% St = reshape(sum(phn==1),[j 1]);
% Rt = reshape(sum(phn==2),[j 1]);
% At = reshape(sum(phn==3),[j 1]);

% output = [time' St Rt At];
    
% end


%% Plotting
% 
% minx = min(min(state(:,1,:)));
% maxx = max(max(state(:,1,:)));
% miny = min(min(state(:,2,:)));
% maxy = max(max(state(:,2,:)));
% minz = min(min(state(:,3,:)));
% maxz = max(max(state(:,3,:)));
% maxdim = max(abs([minx maxx miny maxy minz maxz]))+1;
% 
% clrs = NaN(nmax,3);
% 
% figure('Position',[10 10 600 600])
% for i = 1:size(state,3)
%     
%     Srows = state(:,4,i)==1;
%     Rrows = state(:,4,i)==2;
%     Arows = state(:,4,i)==3;
%     
%     clrs(Srows,:) = repmat(clr_sen,sum(Srows),1);
%     clrs(Rrows,:) = repmat(clr_res,sum(Rrows),1);
%     clrs(Arows,:) = repmat(clr_adj,sum(Arows),1);
% 
%     scatter3(state(:,1,i),state(:,2,i),state(:,3,i),5e2,clrs,'.');
%     hold on
%     scatter3(pos0_caf(:,1),pos0_caf(:,2),pos0_caf(:,3),'k*')    % plot CAFs
%     scatter3(state(:,1,i),state(:,2,i),repmat(-maxdim,size(state,1),1),5e2,[.9 .9 .9],'.'); % add shadow
%     xlim([-maxdim maxdim]);
%     ylim([-maxdim maxdim]);
%     zlim([-maxdim maxdim]);
%     view(-75+i/5,20)
%     title(strcat('Time: ',num2str(round(output(i,2),0)),' Days'))
%     grid off
%     set(gca,'TickLength',[0 0],'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);%,'visible','off')
%     hold off
%     pause(1e-2)
%     
%     F(i) = getframe(gcf);
%     drawnow    
% end
% 
% % create the video writer with 1 fps
% writerObj = VideoWriter('ExampleSim_091120b_bdgmax30.mp4','MPEG-4');
% writerObj.FrameRate = 20;
% 
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);

% time = cumsum([0 t]);
% 
% phn = state(:,4,:);
% St = reshape(sum(phn==1),[j 1]);
% Rt = reshape(sum(phn==2),[j 1]);
% At = reshape(sum(phn==3),[j 1]);
% 
% figure
% plot(time,[St Rt At])
% xlabel('Days')
% ylabel('Cell Number')
% legend({'sen' 'res' 'adj'},'Location','northeastoutside')
