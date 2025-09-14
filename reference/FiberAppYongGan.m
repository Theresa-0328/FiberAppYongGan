function FiberApp勇敢实验
    % Fiber photometry GUI for "越挫越勇" paradigm
    % 保存为 FiberApp勇敢实验.m，运行即可

    % 主界面
    fig = uifigure('Name','Fiber 越挫越勇实验','Position',[100 100 1000 600]);

    % 数据存储
    app.FiberData = []; % 光纤数据
    app.dFF = [];
    app.SamplingRate = 20; % 默认采样率 Hz（可改）

    % 标签页组
    tg = uitabgroup(fig,'Position',[20 20 960 560]);

    %% Tube 范式
    tab1 = uitab(tg,'Title','Tube 范式');
    colnames1 = {'TrialID','管径','障碍物','StartDrill','Escape','PushObstacle','ExitTube','Outcome'};
    data1 = cell(0,numel(colnames1));
    app.tbl1 = uitable(tab1,'Data',data1,...
        'ColumnName',colnames1,'Position',[20 150 900 350],...
        'ColumnEditable',true(1,numel(colnames1)));

    % 按钮
    uibutton(tab1,'Text','开始钻管','Position',[20 50 100 30],...
        'ButtonPushedFcn',@(src,event)markEvent(app.tbl1,'StartDrill'));
    uibutton(tab1,'Text','逃避行为','Position',[140 50 100 30],...
        'ButtonPushedFcn',@(src,event)markEvent(app.tbl1,'Escape'));
    uibutton(tab1,'Text','推障碍物','Position',[260 50 100 30],...
        'ButtonPushedFcn',@(src,event)markEvent(app.tbl1,'PushObstacle'));
    uibutton(tab1,'Text','出管','Position',[380 50 100 30],...
        'ButtonPushedFcn',@(src,event)markEvent(app.tbl1,'ExitTube'));

    %% Water tank 范式
    tab2 = uitab(tg,'Title','Water tank 范式');
    colnames2 = {'TrialID','水位','JumpOn','Escape','EnterWater','JumpOut','Outcome'};
    data2 = cell(0,numel(colnames2));
    app.tbl2 = uitable(tab2,'Data',data2,...
        'ColumnName',colnames2,'Position',[20 150 900 350],...
        'ColumnEditable',true(1,numel(colnames2)));

    uibutton(tab2,'Text','跳上水槽','Position',[20 50 100 30],...
        'ButtonPushedFcn',@(src,event)markEvent(app.tbl2,'JumpOn'));
    uibutton(tab2,'Text','逃避行为','Position',[140 50 100 30],...
        'ButtonPushedFcn',@(src,event)markEvent(app.tbl2,'Escape'));
    uibutton(tab2,'Text','进入水槽','Position',[260 50 100 30],...
        'ButtonPushedFcn',@(src,event)markEvent(app.tbl2,'EnterWater'));
    uibutton(tab2,'Text','跳出水槽','Position',[380 50 100 30],...
        'ButtonPushedFcn',@(src,event)markEvent(app.tbl2,'JumpOut'));

    %% 公共功能区
    % 加载数据按钮
    uibutton(fig,'Text','加载光纤数据','Position',[750 520 100 30],...
        'ButtonPushedFcn',@(src,event)loadFiberData);

    % 对齐绘图按钮
    uibutton(fig,'Text','对齐并绘图','Position',[870 520 100 30],...
        'ButtonPushedFcn',@(src,event)alignAndPlot);

    % 绘图区域
    app.ax = uiaxes(fig,'Position',[500 50 450 450],'Title','dF/F 信号');
    
    %% 嵌套函数：加载数据
    function loadFiberData
        [file,path] = uigetfile({'*.xlsx;*.csv','Data Files'},'选择光纤数据');
        if isequal(file,0), return; end
        T = readtable(fullfile(path,file));
        % 假设列名：F470, F415
        if all(ismember({'F470','F415'}, T.Properties.VariableNames))
            F470 = T.F470;
            F415 = T.F415;
            % 简单线性回归
            p = polyfit(F415,F470,1);
            Ffit = polyval(p,F415);
            dFF = (F470 - Ffit)./Ffit * 100;
            app.FiberData = T;
            app.dFF = dFF;
            plot(app.ax, dFF); title(app.ax,'全程 dF/F'); xlabel(app.ax,'Time'); ylabel(app.ax,'dF/F (%)');
        else
            uialert(fig,'数据必须包含 F470 和 F415 列','错误');
        end
    end

    %% 嵌套函数：标注事件
    function markEvent(tbl,eventName)
        t = now; % 当前时间
        tStr = datestr(t,'HH:MM:SS.FFF'); % 转换成字符串

        if isempty(tbl.Data)
            trialID = 1;
        else
            trialID = size(tbl.Data,1) + 1;
        end

        colNames = tbl.ColumnName;
        row = cell(1,numel(colNames));
        row{strcmp(colNames,'TrialID')} = trialID;

        if ismember(eventName,colNames)
            row{strcmp(colNames,eventName)} = tStr;
        end

        if strcmp(eventName,'ExitTube') || strcmp(eventName,'JumpOut')
            row{strcmp(colNames,'Outcome')} = 'pass';
        else
            row{strcmp(colNames,'Outcome')} = 'ongoing';
        end

        tbl.Data = [tbl.Data; row];
    end

    %% 嵌套函数：对齐绘图
    function alignAndPlot
        if isempty(app.dFF)
            uialert(fig,'请先加载光纤数据','提示');
            return;
        end

        % 简单示例：以 Tube 表格的 StartDrill 时间为事件点
        if isempty(app.tbl1.Data)
            uialert(fig,'请先在 Tube 范式里标注至少一个 StartDrill','提示');
            return;
        end

        % 随机假设时间戳对应采样点（演示用，需改成真实视频→信号时间同步）
        n = numel(app.dFF);
        idx = round(n/2); % 取信号中间点当事件
        win = -app.SamplingRate*5 : app.SamplingRate*10; % -5 到 +10 秒
        idxs = idx + win;
        idxs(idxs<1 | idxs>n) = [];
        peri = app.dFF(idxs);

        plot(app.ax,win(1:numel(peri))/app.SamplingRate, peri,'LineWidth',2);
        xlabel(app.ax,'Time (s)');
        ylabel(app.ax,'dF/F (%)');
        title(app.ax,'事件对齐 dF/F 示例');
    end
end
