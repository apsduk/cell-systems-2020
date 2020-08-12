%% %%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('all'); close('all');

% model id
modelid{ 1} = 'uber_or';
modelid{ 2} = 'frag_or';
modelid{ 3} = 'sqtr_or';
modelid{ 4} = 'op_or';
modelid{ 5} = 'var_frag_or';
modelid{ 6} = 'ortho_op_or';

% number of best controllers
N0 = 200;

% determine output threshold
outputthreshold = 1;

% max/min threshold
percdeltathreshold = 1e-3;
absdeltathreshold = 1e3;

% make figure
fsort = figure('Units', 'normalized', 'Position', [0 0 1 1]);
fplot = figure('Units', 'normalized', 'Position', [0 0 1 1]);
cmap = lines(length(modelid));

%% %%%%% ITERATE OVER MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:length(modelid)
    
    %%%%%%%% LOAD DATA AND OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % load coupling data form sims
    load(['model_',modelid{m},'_coupling.mat']);
        
    % load simualtion settings
    load(['model_',modelid{m},'_simulation_settings.mat'],'X0','X00');
        
    % load max and min
    load(['model_',modelid{m},'_max_min.mat']);
    
    % get the index list
    originalidx = 1:length(p1_end);
    
    %%%%%%%% REMOVE MODEL 1 CONTORLLERS WHERE OMEGA_g_Q AND F ARE ZERO %%%%
    if strcmp(modelid{m}, 'uber_or')
        
        % get designs
        omegagF = X00(:,3); omegagQ = X00(:,4);
        
        % where are these parameters not zero
        idx_F = find(omegagF ~= 0); idx_Q = find(omegagQ ~= 0);
        
        % indicies where omegagQ and F are not zero
        keepidx = intersect(idx_F, idx_Q);
        
    else
        
        keepidx = 1:length(p1_end)';
        
    end
        
    %%%%%%%% DETERMINE WHICH DESIGNS HAVE MAX/MIN DEFINE AS OCILLATORS %%%%
    
    % calculate difference between max and min for first induction
    delta_1 = p1maxval_1 - p1minval_1;
    
    % calculate difference between max and min for second induction
    delta_2 = p1maxval_2 - p1minval_2;
    
    % calculate the numerical threshold
    calcdeltathreshold = percdeltathreshold*max([p1maxval_1, p1maxval_2], [], 2);
    
    % if threshold is zero
    calcdeltathreshold(calcdeltathreshold == 0) = 1e-3;
    
    % if threshold > abs threshold then set to absthreshold
    calcdeltathreshold(calcdeltathreshold > absdeltathreshold) = absdeltathreshold;
    
    % find where diff between max and min is less than threshold
    maxminidx = find( (delta_1 < calcdeltathreshold) + (delta_2 < calcdeltathreshold) == 2);
        
    %%%%%%%% REMOVE RESULTS WHERE OUTPUT IS EFFECTIVELY ZERO %%%%%%%%%%%%%%
    
    % check if the expression is greater than effectively zero
    outputidx = find(p1_end > outputthreshold);
    
    %%%%%%%% CALCAULTE ANALYSIS STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % mRNA coupling
    m_delta = (m1_end - m1_ind)./m1_ind;
    m_output = m2_end;
    
    % protein coupling
    p_delta = (p1_end - p1_ind)./p1_ind;
    p_output = p2_end;
    p_inter_output = p1_ind;
    p_final_output = p2_end;
    
    % scale the coupling scores
    m_scaled = m_delta; % /max(abs(m_delta));
    p_scaled = p_delta; % /max(abs(p_delta));
    
    % determine distance from perfection (0, 0) as Euclidean distance
    score2d = sqrt(m_scaled.^2 + p_scaled.^2);
    
    %%%%%%%% SELECT THE BEST PERFORMERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % The best controllers are the top N controllers which are
    %   stable i.e. not in unstableidx
    %   have non-negligable expression i.e. are in outputidx
    %   are close to (0, 0) as determined by ridx
    %   and for model 1 are not in the discount list
    
    %%%%%%%% SELECT THE VALID CONTROLLERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get the keep controllers
    originalidx = originalidx(ismember(originalidx, keepidx));
    
    % remove controllers whose deriv is greater than odethreshold
    % originalidx = originalidx(ismember(originalidx, steadystateidx));
    
    % get the stable controllers
    % originalidx = originalidx(ismember(originalidx, stableidx));
    
    % get controllers with small max/min differences
    originalidx = originalidx(ismember(originalidx, maxminidx));
    
    % get the controllers which good expression
    originalidx = originalidx(ismember(originalidx, outputidx));
    
    % plot results
    figure(fplot.Number); box('on'); hold('on'); grid('on');
    plot3(m_delta(originalidx), p_delta(originalidx), p_output(originalidx), '.', 'Color', cmap(m,:));
    plot3([0, 0], [0, 0], [0, 4e5], 'k', 'LineWidth', 2)
    xlabel('mRNA delta')
    ylabel('Protein delta');
    zlabel('Protein output');
    xlim([-1, 1]);
    ylim([-1, 1]);
    zlim([0, 4e5]);
    view([-50 40]);
    savefig('design_performance.fig');
    
    %%%%%%%% SELECT THE BEST CONTROLLERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get Euclidean distances from (0, 0)
    topscore = score2d(originalidx);
    
    % get the indicies sorted by the controllers by distance from (0, 0)
    [~, ridx] = sort(topscore, 'ascend');
    
    if length(ridx) > N0
        N = N0;
    else
        N = length(ridx);
    end
    
    % now sort the contorllers by ridx
    originalidx = originalidx(ridx(1:N));
    
    % select best results
    bestidx = originalidx;
    bestX0 = X0(originalidx,:);
    bestX00 = X00(originalidx,:);
    m_delta = m_delta(originalidx);
    m_output = m_output(originalidx);
    p_delta = p_delta(originalidx);
    p_output = p_output(originalidx);
    topscore = score2d(originalidx);
    
    %%%%%%%% SAVE RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % save contorller designs
    save(['model_',modelid{m},'_best_controllers.mat'],'bestidx','bestX0','bestX00');
    
    % sort by output
    [~, outputsort] = sort(p_output, 'descend');
    
    % sort by coupling
    [~, couplingsort] = sort(topscore, 'ascend');
    
    % couple score
    coupling2dscore = topscore;
    
    % output score
    mrnadelta = m_delta;
    mrnaoutput = m_output;
    proteindelta = p_delta;
    proteinintermediateoutput = p_inter_output(originalidx);
    proteinfinaloutput = p_final_output(originalidx);
    
    % save all metrics
    save(['model_',modelid{m},'_different_metrics.mat'],'bestidx','bestX0','bestX00','outputsort','couplingsort', ...
        'coupling2dscore','mrnadelta','mrnaoutput','proteindelta','proteinintermediateoutput','proteinfinaloutput');
    
    %%%%%%%% PLOT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % plot orders
    figure(fsort.Number); subplot(1, length(modelid), m);
    plot(outputsort, couplingsort, '.'); xlabel('Output sort idx'); ylabel('Coupling sort idx');
    title(['modelid = ',modelid{m}]); pbaspect([1, 1, 1]);
    savefig('design_model_score_sortint');
    
end

close('all');