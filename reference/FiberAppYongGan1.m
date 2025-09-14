function FiberAppYongGan
% FiberAppYongGan - 完整 GUI 程序（分三段粘贴，按顺序）
% 1) 放到 FiberAppYongGan.m
% 2) 在 MATLAB 命令行运行 FiberAppYongGan

%% ========== 用户参数（可修改） ==========
sr = 30;                 % 采样率 (Hz) — 若你的数据不同请修改
alignWin = [-2 6];       % 对齐窗口（秒）
padValue = NaN;          % 对齐窗口不足时填充值

%% ========== GUI 创建 ==========
fig = figure('Name','Fiber 越挫越勇 实验 (GUI)','Position',[100 100 1200 720],...
    'MenuBar','none','ToolBar','none');

% 右侧快速显示区
axQuick = axes('Parent',fig,'Position',[0.55 0.55 0.4 0.38]);
xlabel(axQuick,'Time (s)'); ylabel(axQuick,'dF/F (%)'); title(axQuick,'Quick plot');

% 左侧两个 table（Tube / Water Tank）
col_tube = {'TrialID','管径','障碍物','StartDrill','Escape','PushObstacle','ExitTube','Outcome'};
tblTube = uitable('Parent',fig,'Units','normalized','Position',[0.03 0.45 0.5 0.5],...
    'Data',cell(0,numel(col_tube)),'ColumnName',col_tube,'Tag','TubeTable');

col_tank = {'TrialID','水位','JumpOn','Escape','EnterWater','JumpOut','Outcome'};
tblTank = uitable('Parent',fig,'Units','normalized','Position',[0.03 0.02 0.5 0.38],...
    'Data',cell(0,numel(col_tank)),'ColumnName',col_tank,'Tag','TankTable');

% 控件按钮
uicontrol('Parent',fig,'Style','pushbutton','String','加载光纤数据',...
    'Units','normalized','Position',[0.55 0.92 0.12 0.05],'Callback',@onLoadFile);

uicontrol('Parent',fig,'Style','pushbutton','String','对齐并绘图',...
    'Units','normalized','Position',[0.69 0.92 0.12 0.05],'Callback',@onAlignAndPlot);

uicontrol('Parent',fig,'Style','pushbutton','String','导出 (当前文件)',...
    'Units','normalized','Position',[0.83 0.92 0.12 0.05],'Callback',@onExportCurrent);

uicontrol('Parent',fig,'Style','pushbutton','String','批量处理光纤数据',...
    'Units','normalized','Position',[0.55 0.86 0.4 0.05],'Callback',@onBatchProcess);

% 简易事件标注按钮（示例，写入表格）
uicontrol('Parent',fig,'Style','pushbutton','String','开始钻管','Units','normalized','Position',[0.03 0.9 0.1 0.04],...
    'Callback',@(s,e)markEvent('Tube','StartDrill'));
uicontrol('Parent',fig,'Style','pushbutton','String','逃避','Units','normalized','Position',[0.14 0.9 0.1 0.04],...
    'Callback',@(s,e)markEvent('Tube','Escape'));
uicontrol('Parent',fig,'Style','pushbutton','String','推障碍物','Units','normalized','Position',[0.25 0.9 0.1 0.04],...
    'Callback',@(s,e)markEvent('Tube','PushObstacle'));
uicontrol('Parent',fig,'Style','pushbutton','String','出管/跳出','Units','normalized','Position',[0.36 0.9 0.1 0.04],...
    'Callback',@(s,e)markEvent('Tube','ExitTube'));

%% ========== 状态变量 ==========
state = struct();
state.sr = sr;
state.alignWin = alignWin;
state.winSamples = round(alignWin * sr);
state.padValue = padValue;
state.loadedTable = [];    % 最近加载的表格 (readtable)
state.signals = struct();  % 存放计算后的 dF/F： fieldnames 为 原始列名的有效 MATLAB 名称
state.lastLoadedFile = '';

%% ========== 主回调 （在 Part 2/3 中实现子函数） ==========
    function onLoadFile(~,~)
        % 加载单个文件，自动识别光纤列，计算 dF/F 并绘图（Quick plot）
        [file,path] = uigetfile({'*.xlsx;*.csv','Data Files'},'选择单个数据文件');
        if isequal(file,0), return; end
        fullf = fullfile(path,file);
        try
            T = readtable(fullf,'PreserveVariableNames',true);
        catch ME
            errordlg(['读取文件失败: ' ME.message]);
            return;
        end
        state.loadedTable = T;
        state.lastLoadedFile = fullf;
        % 子函数 detectChannels / compute_dff 在 Part 2 定义
        chans = detectChannels(T);
        if isempty(chans)
            warndlg('未检测到可识别光纤列（如：470/410/CH1/CH2/Signal1/Signal2）');
            return;
        end
        % 计算并展示
        cla(axQuick); hold(axQuick,'on');
        state.signals = struct();
        for i = 1:numel(chans)
            col = chans{i};
            sigRaw = safeGetVar(T,col);
            dff = compute_dff(sigRaw, T, col);
            fld = matlab.lang.makeValidName(col);
            state.signals.(fld).dff = dff;
            % quick plot: CH1 红, CH2 蓝
            if i==1, plot(axQuick,(0:numel(dff)-1)/state.sr,dff,'r','LineWidth',1.2);
            else plot(axQuick,(0:numel(dff)-1)/state.sr,dff,'b','LineWidth',1.2); end
        end
        hold(axQuick,'off');
        legend(axQuick,fieldnames(state.signals),'Interpreter','none','Location','best');
        title(axQuick, sprintf('Loaded: %s', file),'Interpreter','none');
    end

    function onAlignAndPlot(~,~)
        % 对当前加载的数据以 center 对齐并在新图显示（用于快速检查）
        sFields = fieldnames(state.signals);
        if isempty(sFields)
            errordlg('请先加载光纤数据文件');
            return;
        end
        figure('Name','Align & Quick Plot','NumberTitle','off','Color','w');
        hold on;
        colors = {'r','b','m','k'};
        win = state.winSamples;
        for i = 1:numel(sFields)
            dff = state.signals.(sFields{i}).dff;
            idxCenter = round(numel(dff)/2);
            peri = getPeri(dff, idxCenter, win);
            tvec = (0:numel(peri)-1)/state.sr + state.alignWin(1);
            plot(tvec, peri, 'Color', colors{min(i,numel(colors))}, 'LineWidth', 1.5);
        end
        xlabel('Time (s)'); ylabel('dF/F (%)'); title('Aligned traces (center)'); legend(sFields,'Interpreter','none');
        hold off;
    end

    function onExportCurrent(~,~)
        % 导出当前加载文件的 Peak/AUC（以 center 对齐）
        sFields = fieldnames(state.signals);
        if isempty(sFields)
            errordlg('请先加载光纤数据');
            return;
        end
        [file,path] = uiputfile('current_results.xlsx','保存当前文件结果');
        if isequal(file,0), return; end
        fullOut = fullfile(path,file);
        rows = {};
        for i = 1:numel(sFields)
            dff = state.signals.(sFields{i}).dff;
            peri = getPeri(dff, round(numel(dff)/2), state.winSamples);
            tvec = (0:numel(peri)-1)/state.sr + state.alignWin(1);
            rows(end+1,:) = {sFields{i}, nanmax(peri), trapz(tvec, peri)}; %#ok<AGROW>
        end
        Tsave = cell2table(rows, 'VariableNames', {'Channel','Peak','AUC'});
        try
            writetable(Tsave, fullOut);
            msgbox(['已保存: ' fullOut]);
        catch
            errordlg('保存失败，请检查路径权限');
        end
    end

    function markEvent(tabName,eventName)
        % 把当前时间戳写入对应表格（简单手动标注）
        if strcmpi(tabName,'Tube')
            tbl = findobj(fig,'Tag','TubeTable');
        else
            tbl = findobj(fig,'Tag','TankTable');
        end
        colNames = tbl.ColumnName;
        newRow = cell(1,numel(colNames));
        if isempty(tbl.Data), trialID = 1; else trialID = size(tbl.Data,1)+1; end
        newRow{strcmp(colNames,'TrialID')} = trialID;
        newRow{strcmp(colNames,eventName)} = datestr(now,'yyyy-mm-dd HH:MM:SS.FFF');
        if any(strcmp(eventName,{'ExitTube','JumpOut'})), newRow{strcmp(colNames,'Outcome')}='pass';
        else newRow{strcmp(colNames,'Outcome')}='ongoing'; end
        tbl.Data = [tbl.Data; newRow];
    end

% Part 1 ends here. Continue with Part 2 and Part 3 in order.
%% =================== Part 2: 单文件处理 helper 函数 ===================
% 这些嵌套函数依赖于外面定义的 state, sr, alignWin 等

    function val = safeGetVar(T, name)
        % 从 table 安全读取列（保留原始列名）
        if ismember(name, T.Properties.VariableNames)
            val = T.(name);
        else
            % 退回到 MATLAB 合法名查找
            safe = matlab.lang.makeValidName(name);
            if ismember(safe, T.Properties.VariableNames)
                val = T.(safe);
            else
                error('列 %s 在表中不存在', name);
            end
        end
    end

    function cols = detectChannels(T)
        % 自动检测表 T 中可能代表光纤的列（返回原始列名）
        % 优先：包含 '470' 或 'ch1' 等，再考虑 '410'/'ch2'
        varnames = T.Properties.VariableNames;
        lc = lower(varnames);
        cols = {};
        % 首选 470
        idx470 = find(contains(lc,'470'),1);
        if ~isempty(idx470), cols{end+1} = varnames{idx470}; end
        % then ch1 / signal1
        idxch1 = find(contains(lc,{'ch1','signal1','channel1'}),1);
        if ~isempty(idxch1) && ~ismember(varnames{idxch1},cols), cols{end+1} = varnames{idxch1}; end
        % then 410 or ch2
        idx410 = find(contains(lc,'410'),1);
        if ~isempty(idx410) && ~ismember(varnames{idx410},cols), cols{end+1} = varnames{idx410}; end
        idxch2 = find(contains(lc,{'ch2','signal2','channel2'}),1);
        if ~isempty(idxch2) && ~ismember(varnames{idxch2},cols), cols{end+1} = varnames{idxch2}; end
        % fallback: choose first two numeric columns (exclude obvious time columns)
        if numel(cols) < 1
            for k = 1:numel(varnames)
                try
                    v = T.(varnames{k});
                catch
                    v = [];
                end
                if isnumeric(v) && numel(v) > 10 && ~contains(lower(varnames{k}), 'time')
                    cols{end+1} = varnames{k}; %#ok<AGROW>
                end
                if numel(cols) >= 2, break; end
            end
        end
        % 保证不多于两个
        if numel(cols) > 2, cols = cols(1:2); end
    end

    function dff = compute_dff(sig, T, colName)
        % 计算 dF/F：
        % - 若表中同时存在 470 (signal) 和 410 (isosbestic) 或 CH1/CH2，把 410 回归后扣除
        % - 否则使用 baseline (10th percentile) 作为基线
        sig = double(sig(:));
        varnames = T.Properties.VariableNames;
        lc = lower(varnames);
        % 尝试定位配对的 isosbestic (410/CH2等)
        paired = [];
        if contains(lower(colName),'470') || contains(lower(colName),'ch1') || contains(lower(colName),'signal1')
            % 寻找 410 / ch2
            idx410 = find(contains(lc,'410'),1);
            if isempty(idx410), idx410 = find(contains(lc,{'ch2','signal2','channel2'}),1); end
            if ~isempty(idx410), paired = varnames{idx410}; end
        end
        if ~isempty(paired)
            try
                ref = double(T.(paired)(:));
                p = polyfit(ref, sig, 1);
                fitted = polyval(p, ref);
                dff = (sig - fitted) ./ fitted * 100;
            catch
                base = prctile(sig,10);
                dff = (sig - base) ./ base * 100;
            end
        else
            base = prctile(sig,10); % 更稳健的 baseline（10%分位）
            if base == 0, base = mean(sig); end
            dff = (sig - base) ./ base * 100;
        end
    end

    function peri = getPeri(sig, centerIdx, winSamples)
        % 返回中心 index 为 centerIdx 的 [winSamples(1) .. winSamples(2)] 样本范围信号（可能短）
        s = winSamples(1); e = winSamples(2);
        n = numel(sig);
        i1 = max(1, centerIdx + s);
        i2 = min(n, centerIdx + e);
        peri = sig(i1:i2);
    end

    function vec = padToLength(vec, L, padVal)
        % 将行向量 pad 到长度 L
        vec = vec(:)';
        if numel(vec) < L
            vec = [vec repmat(padVal,1,L-numel(vec))];
        else
            vec = vec(1:L);
        end
    end

    function out = sanitizeFilename(s)
        % 把字符串变成文件名友好形式
        out = regexprep(s,'[^\w\-\.]','_');
    end
%% =================== Part 3: 批量处理与汇总保存 ===================

    function onBatchProcess(~,~)
        root = uigetdir(pwd,'选择批量处理根目录（含每只小鼠子文件夹）');
        if isequal(root,0), return; end
        resultsRoot = fullfile(root,'results');
        if ~exist(resultsRoot,'dir'), mkdir(resultsRoot); end

        mice = dir(root);
        mice = mice([mice.isdir] & ~ismember({mice.name},{'.','..','results'}));
        if isempty(mice), errordlg('根目录内未发现小鼠文件夹'); return; end

        AllTrials = {}; % rows: Mouse, Paradigm, Date, File, Channel, EventName, Peak, AUC, TrialExcelPath
        % iterate mice
        waitFig = waitbar(0,'批量处理开始...');
        totalSteps = numel(mice);
        step = 0;
        for mi = 1:numel(mice)
            mouseID = mice(mi).name;
            safeMouseID = matlab.lang.makeValidName(mouseID);
            mouseFolder = fullfile(root, mouseID);
            % 找到范式文件夹（case-insensitive匹配 'tube' 和 'water tank'）
            subdirs = dir(mouseFolder);
            subdirs = subdirs([subdirs.isdir] & ~ismember({subdirs.name},{'.','..'}));
            % attempt find tube and water tank (case-insensitive)
            paradigmNames = {};
            for sd = 1:numel(subdirs)
                n = lower(subdirs(sd).name);
                if contains(n,'tube') || contains(n,'tank') || contains(n,'water')
                    paradigmNames{end+1} = subdirs(sd).name; %#ok<AGROW>
                end
            end
            if isempty(paradigmNames), continue; end

            for pn = 1:numel(paradigmNames)
                paradigm = paradigmNames{pn};
                parFolder = fullfile(mouseFolder, paradigm);
                dateDirs = dir(parFolder);
                dateDirs = dateDirs([dateDirs.isdir] & ~ismember({dateDirs.name},{'.','..'}));
                if isempty(dateDirs)
                    % sometimes files are directly in parFolder without date subfolders
                    dateDirs = struct('name', {'.'}, 'folder', {parFolder}); % treat parFolder as single date '.'
                end

                % prepare per-mouse+paradigm summary containers
                summaryCH = struct('CH1',{{}}, 'CH2',{{}});

                for di = 1:numel(dateDirs)
                    dateName = dateDirs(di).name;
                    if strcmp(dateName,'.')
                        dataFolder = parFolder;
                    else
                        dataFolder = fullfile(parFolder, dateName);
                    end

                    files = [dir(fullfile(dataFolder,'*.xlsx')); dir(fullfile(dataFolder,'*.csv'))];
                    if isempty(files), continue; end

                    % create out folder
                    outFolder = fullfile(resultsRoot, mouseID, paradigm, dateName);
                    if ~exist(outFolder,'dir'), mkdir(outFolder); end

                    % for each file
                    for fi = 1:numel(files)
                        fname = fullfile(files(fi).folder, files(fi).name);
                        % robust read
                        try
                            T = readtable(fname,'PreserveVariableNames',true);
                        catch
                            warning('无法读取 %s，跳过', fname);
                            continue;
                        end
                        % detect channels
                        chans = detectChannels(T);
                        if isempty(chans)
                            warning('文件 %s 未找到信号列（470/410/CH1/CH2），跳过', files(fi).name);
                            continue;
                        end
                        % detect events (可能为多列)
                        evCols = detectEventColumns(T);

                        % for each channel compute dff and extract peri windows
                        for ci = 1:numel(chans)
                            chName = chans{ci};
                            sigRaw = safeGetVar(T, chName);
                            dff = compute_dff(sigRaw, T, chName);

                            if ~isempty(fieldnames(evCols))
                                evNames = fieldnames(evCols);
                                for ei = 1:numel(evNames)
                                    evName = evNames{ei};
                                    idxs = evCols.(evName);
                                    for jj = 1:numel(idxs)
                                        centerIdx = idxs(jj);
                                        peri = getPeri(dff, centerIdx, state.winSamples);
                                        L = round((state.alignWin(2)-state.alignWin(1))*state.sr) + 1;
                                        peri = padToLength(peri, L, state.padValue);
                                        tvec = linspace(state.alignWin(1), state.alignWin(2), numel(peri));
                                        pk = nanmax(peri);
                                        auc = trapz(tvec, peri);
                                        % save trial excel
                                        base = sprintf('%s_%s_%s_%s_trial%02d_%s', mouseID, sanitizeFilename(paradigm), sanitizeFilename(dateName), sanitizeFilename(files(fi).name), jj, sanitizeFilename(chName));
                                        trialT = table({files(fi).name},{evName},pk,auc,'VariableNames',{'File','Event','Peak','AUC'});
                                        dffCols = array2table(peri, 'VariableNames', strcat('dFF_', string(1:numel(peri))));
                                        trialFull = [trialT dffCols];
                                        trialFile = fullfile(outFolder, [base '.xlsx']);
                                        try
                                            writetable(trialFull, trialFile);
                                        catch
                                            warning('写入 %s 失败', trialFile);
                                        end
                                        % single-trial heatmap
                                        hfig = figure('Visible','off','Position',[100 100 600 140]);
                                        imagesc(tvec, 1, peri);
                                        colormap(redblue(256)); caxis(getSymmetricCaxis(peri)); colorbar;
                                        xlabel('Time (s)'); ylabel('Trial'); title(base,'Interpreter','none');
                                        try
                                            saveas(hfig, fullfile(outFolder, [base '_heatmap.png']));
                                        catch
                                            warning('保存 %s 失败', fullfile(outFolder, [base '_heatmap.png']));
                                        end
                                        close(hfig);
                                        % append to AllTrials summary
                                        AllTrials(end+1,:) = {mouseID, paradigm, dateName, files(fi).name, chName, evName, pk, auc, trialFile}; %#ok<AGROW>
                                        % add to per-mouse summary container
                                        if ci==1
                                            summaryCH.CH1{end+1} = peri; %#ok<AGROW>
                                        else
                                            summaryCH.CH2{end+1} = peri; %#ok<AGROW>
                                        end
                                    end
                                end
                            else
                                % center align
                                centerIdx = round(numel(dff)/2);
                                peri = getPeri(dff, centerIdx, state.winSamples);
                                L = round((state.alignWin(2)-state.alignWin(1))*state.sr) + 1;
                                peri = padToLength(peri, L, state.padValue);
                                tvec = linspace(state.alignWin(1), state.alignWin(2), numel(peri));
                                pk = nanmax(peri); auc = trapz(tvec, peri);
                                base = sprintf('%s_%s_%s_%s_center_%s', mouseID, sanitizeFilename(paradigm), sanitizeFilename(dateName), sanitizeFilename(files(fi).name), sanitizeFilename(chName));
                                trialT = table({files(fi).name},{'center'},pk,auc,'VariableNames',{'File','Event','Peak','AUC'});
                                dffCols = array2table(peri, 'VariableNames', strcat('dFF_', string(1:numel(peri))));
                                trialFull = [trialT dffCols];
                                trialFile = fullfile(outFolder, [base '.xlsx']);
                                try
                                    writetable(trialFull, trialFile);
                                catch
                                    warning('写入 %s 失败', trialFile);
                                end
                                hfig = figure('Visible','off','Position',[100 100 600 140]);
                                imagesc(tvec, 1, peri);
                                colormap(redblue(256)); caxis(getSymmetricCaxis(peri)); colorbar;
                                xlabel('Time (s)'); ylabel('Trial'); title(base,'Interpreter','none');
                                try
                                    saveas(hfig, fullfile(outFolder, [base '_heatmap.png']));
                                catch
                                    warning('保存 %s 失败', fullfile(outFolder, [base '_heatmap.png']));
                                end
                                close(hfig);
                                AllTrials(end+1,:) = {mouseID, paradigm, dateName, files(fi).name, chName, 'center', pk, auc, trialFile}; %#ok<AGROW>
                                if ci==1
                                    summaryCH.CH1{end+1} = peri; %#ok<AGROW>
                                else
                                    summaryCH.CH2{end+1} = peri; %#ok<AGROW>
                                end
                            end
                        end % channel loop
                    end % files loop
                end % date loop

                % save per-mouse/per-paradigm summary (average + heatmap + summary Excel)
                if ~isempty(summaryCH.CH1)
                    mat1 = cell2mat(summaryCH.CH1');
                    tvec = linspace(state.alignWin(1), state.alignWin(2), size(mat1,2));
                    saveSummary(mat1, tvec, resultsRoot, mouseID, paradigm, 'CH1');
                end
                if ~isempty(summaryCH.CH2)
                    mat2 = cell2mat(summaryCH.CH2');
                    tvec = linspace(state.alignWin(1), state.alignWin(2), size(mat2,2));
                    saveSummary(mat2, tvec, resultsRoot, mouseID, paradigm, 'CH2');
                end
            end % paradigms
            step = step + 1;
            waitbar(step/totalSteps, waitFig, sprintf('Processing %s (%d/%d)', mouseID, step, totalSteps));
        end % mice
        close(waitFig);

        % Save AllTrials summary table
        if ~isempty(AllTrials)
            try
                T_all = cell2table(AllTrials, 'VariableNames', {'Mouse','Paradigm','Date','File','Channel','Event','Peak','AUC','TrialFile'});
                writetable(T_all, fullfile(resultsRoot,'AllTrialsSummary.xlsx'));
            catch
                warning('保存 AllTrialsSummary.xlsx 失败');
            end
        end

        msgbox('批量处理完成，结果保存在 results 文件夹');
    end

%% =================== Part 3 Helpers ===================

    function evCols = detectEventColumns(T)
        % 自动识别事件列：包含关键词或 numeric 且稀疏非零
        vn = T.Properties.VariableNames;
        evCols = struct();
        keywords = {'event','ttl','trigger','input','start','stop','poke','lick'};
        for k = 1:numel(vn)
            name = vn{k};
            lname = lower(name);
            try
                col = T.(name);
            catch
                continue;
            end
            if any(contains(lname, keywords))
                if isnumeric(col)
                    idx = find(col > 0);
                else
                    idx = find(~cellfun(@(x) isempty(x) || (isstring(x) && strlength(x)==0), cellstr(col)));
                end
                if ~isempty(idx), evCols.(name) = unique(idx); end
                continue;
            end
            if isnumeric(col)
                nz = find(col ~= 0 & ~isnan(col));
                if ~isempty(nz) && numel(nz) < numel(col)*0.5
                    evCols.(name) = unique(nz);
                end
            end
        end
    end

    function saveSummary(mat, tvec, resultsRoot, mouseID, paradigm, chName)
        % 保存平均曲线图、总体热图和 summary excel (Peak/AUC per trial)
        outDir = fullfile(resultsRoot, mouseID);
        if ~exist(outDir,'dir'), mkdir(outDir); end
        % mean ± sem plot
        meanTrace = nanmean(mat,1);
        semTrace = nanstd(mat,0,1)./sqrt(sum(~isnan(mat),1));
        fig = figure('Visible','off','Position',[200 200 900 400]); hold on;
        fill([tvec fliplr(tvec)], [meanTrace+semTrace fliplr(meanTrace-semTrace)], [0.8 0.8 0.8], 'EdgeColor','none');
        if strcmpi(chName,'CH1'), plot(tvec,meanTrace,'r','LineWidth',2); else plot(tvec,meanTrace,'b','LineWidth',2); end
        xlabel('Time (s)'); ylabel('dF/F (%)'); title([mouseID ' - ' paradigm ' - ' chName], 'Interpreter','none');
        saveas(fig, fullfile(resultsRoot, sprintf('%s_%s_%s_meanTrace.png', mouseID, sanitizeFilename(paradigm), chName)));
        close(fig);
        % heatmap
        fig2 = figure('Visible','off','Position',[200 200 900 400]);
        imagesc(tvec, 1:size(mat,1), mat);
        colormap(redblue(256)); caxis(getSymmetricCaxis(mat)); colorbar;
        xlabel('Time (s)'); ylabel('Trial'); title([mouseID ' - ' paradigm ' - ' chName ' heatmap'], 'Interpreter','none');
        saveas(fig2, fullfile(resultsRoot, sprintf('%s_%s_%s_heatmap.png', mouseID, sanitizeFilename(paradigm), chName)));
        close(fig2);
        % summary excel
        peaks = max(mat,[],2);
        aucs = trapz(repmat(tvec,size(mat,1),1), mat, 2);
        sumT = table((1:size(mat,1))', peaks, aucs, 'VariableNames', {'Trial','Peak','AUC'});
        try
            writetable(sumT, fullfile(resultsRoot, sprintf('%s_%s_%s_summary.xlsx', mouseID, sanitizeFilename(paradigm), chName)));
        catch
            warning('写入 summary excel 失败');
        end
    end

    function cmap = redblue(n)
        if nargin<1, n=256; end
        half = floor(n/2);
        b2w = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
        w2r = [ones(n-half,1), linspace(1,0,n-half)', linspace(1,0,n-half)'];
        cmap = [b2w; w2r];
    end

    function ca = getSymmetricCaxis(mat)
        mn = nanmin(mat(:)); mx = nanmax(mat(:));
        mm = max(abs([mn mx]));
        if mm==0, mm = 1; end
        ca = [-mm mm];
    end

end % end of main function FiberAppYongGan
                           