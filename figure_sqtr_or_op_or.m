%% %%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear('all'); close('all'); addpath('../lib','../libmodels');

hname = pwd;

%% %%%%% SET UP MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
runintmax = 120;
tind = 6;

odesolver = @ode23s;
odeoptions = odeset();

%% %%%%% SET UP MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model{1}.modelid = 'uber_or';
model{1}.Y0 = zeros(25,1);
model{1}.f = @(T, Y, X, Z) uber_or(T, Y, [X(2:9), X(2:9)], X(11:18), X(20:25), X(26:31), [X(10), X(19), Z(1), Z(2), X(44), X(32)], X(33:43), X(45:55));
model{1}.m1idx =  9;
model{1}.p1idx = 11;

model{2}.modelid = 'frag_or';
model{2}.Y0 = zeros(25,1);
model{2}.f = @(T, Y, X, Z) frag_or(T, Y, [X(2:9), X(2:9)], X(11:18), X(20:25), X(26:31), [ Z(1), Z(2), X(10), X(19), X(42)], X(32:41), X(43:53));
model{2}.m1idx = 14;
model{2}.p1idx = 16;

model{3}.modelid = 'sqtr_or';
model{3}.Y0 = zeros(23,1);
model{3}.f = @(T, Y, X, Z) sqtr_or(T, Y, [X(2:9), X(2:9)], X(11:18), X(20:25), X(26:31), [X(10), X(19), Z( 1), Z( 2), X(40), X(32)], X(33:39), X(41:51));
model{3}.m1idx =  9;
model{3}.p1idx = 11;

model{4}.modelid = 'op_or';
model{4}.Y0 = zeros(20,1);
model{4}.f = @(T, Y, X, Z) op_or(T, Y, [X(2:9), X(2:9)], X(11:18), X(20:25), X(26:31), [X(10), X(19), Z(1), Z(2), X(32)], [], X(33:43));
model{4}.m1idx =  9;
model{4}.p1idx = 11;

model{5}.modelid = 'var_frag_or';
model{5}.Y0 = zeros(25,1);
model{5}.f = @(T, Y, X, Z) var_frag_or(T, Y, [X(2:9), X(2:9)], X(11:18), X(20:25), X(26:31), [ Z(1), Z(2), X(10), X(19), X(42)], X(32:41), X(43:53));
model{5}.m1idx = 14;
model{5}.p1idx = 16;

model{6}.modelid = 'ortho_op_or';
model{6}.Y0 = zeros(20,1);
model{6}.f = @(T, Y, X, Z) ortho_op_or(T, Y, [X(2:9), X(2:9)], X(11:18), X(20:25), X(26:31), [X(10), X(19), Z(1), Z(2), X(32)], [], X(33:43));
model{6}.m1idx =  9;
model{6}.p1idx = 11;

%% %%%%% CHOOSE MODELS TO SIMULATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
topo{1}(1) = 1; mlim(1,:) = [0.995, 1.005]; plim(1,:) = [0.9, 1.1]; controllername{1}{1} = 'UBER-OR';
topo{1}(2) = 3; mlim(2,:) = [0.995, 1.005]; plim(2,:) = [0.9, 1.1]; controllername{1}{2} = 'SQTR-OR';
topo{1}(3) = 4; mlim(3,:) = [0.995, 1.005]; plim(3,:) = [0.9, 1.1]; controllername{1}{3} = 'OP-OR';

%% %%%%% ITERATE TOPOLOGIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:length(topo)
    
    topo{t}
    
    fplot = figure('Units', 'centimeters', 'Position', [5 5 30 30]); legendnames = [];
    
    sp{1} = subplot(8, 2, [1, 3]); box('on'); hold('on'); grid('on'); % pbaspect([1 1 1]);
    sp{2} = subplot(8, 2, [2, 4]); box('on'); hold('on'); grid('on'); % pbaspect([1 1 1]);
    sp{3} = subplot(8, 2, 5); hold('on'); box('on');
    sp{4} = subplot(8, 2, 6); hold('on'); box('on');
    sp{5} = subplot(8, 2, 7); hold('on'); box('on');
    sp{6} = subplot(8, 2, 8); hold('on'); box('on');
    sp{7} = subplot(8, 2, [ 9, 11]); box('on'); hold('on');
    sp{8} = subplot(8, 2, [10, 12]); box('on'); hold('on');
    
    cmap = lines(7);
    cmap = [0.5 0.5 0.5; cmap([3:4],:)];
    
    % add legend
    for m = 1:length(topo{t})
        scatter(sp{1}, -1e3, -1e3, 'o', 'MarkerFaceColor', cmap(m,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1);
    end
    
    for m = 1:length(topo{t})
                
        % make submodel
        submodel = model{topo{t}(m)}
        
        %% %%%%% PLOT ALL OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % load performance metrics
        load(['model_',submodel.modelid,'_different_metrics.mat'], 'couplingsort', 'coupling2dscore', 'proteinfinaloutput');
        couplingsort = couplingsort(1:200);
        coupling2dscore = coupling2dscore(1:200);
        proteinfinaloutput = proteinfinaloutput(1:200);
        
        scatter(sp{1}, couplingsort, coupling2dscore, 'o', 'MarkerFaceColor', cmap(m,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); % 0.33);
        set(sp{1}, 'YScale', 'log', 'XLim', [1, length(couplingsort)], 'XTick', [1, length(couplingsort)], 'YLim', [1e-3, 1e0]);
        xlabel(sp{1}, 'Design index (sorted)');
        ylabel(sp{1}, '2D score');
        scatter(sp{2}, proteinfinaloutput, coupling2dscore, 'o', 'MarkerFaceColor', cmap(m,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 1); % 0.33);
        set(sp{2}, 'XScale', 'log', 'YScale', 'log', 'XLim', [1e3, 1e6], 'YLim', [1e-3, 1e0]);
        xlabel(sp{2}, 'Final protein output (nM)');
        ylabel(sp{2}, '2D score');
        
        %% %%%%% SIMULATE PERFORMANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % load parameters
        load(['model_',submodel.modelid,'_best_controllers.mat'], 'bestX0');
        
        % select best parameters
        X = bestX0(1,:);
        
        % make ode functions
        f_0 = @(T, Y) submodel.f(T, Y, X, [   0,    0]);
        f_1 = @(T, Y) submodel.f(T, Y, X, [X(1),    0]);
        f_2 = @(T, Y) submodel.f(T, Y, X, [X(1), X(1)]);
        
        % simulate odes
        [T_0, Y_0] = odesolver(f_0, [- runintmax, 0], submodel.Y0, odeoptions);
        [T_0, Y_0] = odesolver(f_1, [- runintmax, 0], Y_0(end,:)', odeoptions);
        [T_1, Y_1] = odesolver(f_1, [0, tind], Y_0(end,:)', odeoptions);
        [T_2, Y_2] = odesolver(f_2, [max(T_1) + (1e-6), max(T_1) + 3*tind], Y_1(end,:)', odeoptions);
        T = [T_0; T_1; T_2];
        Y = [Y_0; Y_1; Y_2];
        
        hold(sp{2*m - 1}, 'on');
        if m == 1
            plot(sp{3}, T, Y(:,submodel.m1idx)./Y_1(end,submodel.m1idx), ':', 'Color', cmap(m,:), 'LineWidth', 2);
            plot(sp{5}, T, Y(:,submodel.m1idx)./Y_1(end,submodel.m1idx), ':', 'Color', cmap(m,:), 'LineWidth', 2);
        else
            plot(sp{2*m - 1}, T, Y(:,submodel.m1idx)./Y_1(end,submodel.m1idx), '-', 'Color', cmap(m,:), 'LineWidth', 2);
            set(sp{2*m - 1}, 'YLim', mlim(m,:), 'XLim', [0 max(T)]);
            ylabel(sp{2*m - 1}, '\Delta TX');
            set(sp{2*m - 1}, 'XTick', []);
        end
        if 2*m - 1 == 5
            set(sp{2*m - 1}, 'XTick', [0:12:max(T)]);
            xlabel(sp{2*m - 1}, 'Time (h)');
        end
        
        hold(sp{2*m}, 'on');
        if m == 1
            plot(sp{4}, T, Y(:,submodel.p1idx)./Y_1(end,submodel.p1idx), ':', 'Color', cmap(m,:), 'LineWidth', 2);
            plot(sp{6}, T, Y(:,submodel.p1idx)./Y_1(end,submodel.p1idx), ':', 'Color', cmap(m,:), 'LineWidth', 2);
        else
            plot(sp{2*m}, T, Y(:,submodel.p1idx)./Y_1(end,submodel.p1idx), '-', 'Color', cmap(m,:), 'LineWidth', 2);
            set(sp{2*m}, 'YLim', plim(m,:), 'XLim', [0, max(T)]);
            ylabel(sp{2*m}, '\Delta TL');
            set(sp{2*m}, 'XTick', []);
        end
        if 2*m == 6
            set(sp{2*m}, 'XTick', [0:12:max(T)]);
            xlabel(sp{2*m}, 'Time (h)');
        end
        
        bar(sp{7}, m, Y(end, submodel.m1idx), 'FaceColor', cmap(m,:));
        set(sp{7}, 'YScale', 'log', 'YLim', [1e3 10^6.5], 'XTick', 1:1:length(controllername{t}), 'XTickLabel', string(1:4)); % string(controllername{t}));
        ylabel(sp{7}, 'mRNA output (nM)');
        bar(sp{8}, m, Y(end, submodel.p1idx), 'FaceColor', cmap(m,:));
        set(sp{8}, 'YScale', 'log', 'YLim', [1e3 1e6], 'XTick', 1:1:length(controllername{t}), 'XTickLabel', string(1:4)); % string(controllername{t}));
        ylabel(sp{8}, 'Protein output (nM)');
        
    end
    
    for s = 1:length(sp)
        sp{s}.LineWidth = 1;
        sp{s}.FontSize = 12;
        sp{s}.FontWeight = 'bold';
    end
    
    legend(sp{1}, controllername{1}, 'NumColumns', 3, 'FontSize', 12);
    
end