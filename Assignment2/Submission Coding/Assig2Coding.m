%% Setup
clc;clear;

% Set the working directory to the place where the current file is saved
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

dt = readtable("Assig2.csv", 'VariableNamingRule','preserve');

%% (a)(i) OlS Estimator
% for regression Corruption = β1 + β2ln(gdp) + β3Socialnetshare + ε.

y= dt.corruption;
gdp = dt.lngdp;
socialnetshare = dt.socialnetshare;

X = [gdp socialnetshare];

% linear regression using fitlm
ols = fitlm(X, y);

dt2 = dt(~any(ismissing(dt),2),:);

%% (a)(i) Report Output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Please use Matlab Add-on             %
%                makeLatexTable                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

latexCode = convert_LMM2latex(...
    ols,... % array of glme results
    {},... % name of the models
    {'n','adjR2','BIC'},... % information required
    3,... % round digits
    'test_lm',... % table label
    'vertical'); % orientation (vertical or rotated)

% In case you want to replace the name of certain variables:
latexCode=replaceWords(latexCode,...
    {'Insert Title Here.','var1','var2','var3'},... 
    {'Question 4(a), Linear Regression, OLS Estimates',...
    'error','trial','difficulty'}); 

% Save the table
fid=fopen('test_lm.txt','w'); 
fprintf(fid,'%s\n',latexCode{:}); fclose(fid);

%% (a)(ii) SST, SSR and SSE

% can calculate manually and compare with fitlm's result:
% extract result to compare:
sse_fitlm = ols.SSE;
sst_fitlm = ols.SST;
ssr_fitlm = ols.SSR;
beta = ols.Coefficients.Estimate;
%
% create regressors' matrix
X = [ones(47,1) gdp socialnetshare];
% eliminate NaNs, cleaned data to be named `newdat` 
newdat = [y X];
newdat = newdat(~any(ismissing(newdat),2),:);
newdat = array2table(newdat);
newdat = renamevars(newdat,...
    ["newdat1","newdat2","newdat3","newdat4"], ...
                 ["corruption","x1","lngdp",'socialnetshare']);
%
% new regressos; matrix
Xnew =[ones(35,1) newdat.lngdp newdat.socialnetshare];
% by hand OLS estimates
b = Xnew\newdat.corruption;
    % ---> recall b = beta, 
    % so fitlm result is the same with manual computing
%
% Compute residual: y-X*\hat{\beta}
ehat = y-(X*b);
%
% eliminate omitted entries
ehat=ehat(~any(ismissing(ehat),2),:);

% SSE, Sum of squares error: \epsilon^{\prime}*\epsilon
SSE = ehat'*ehat; 
    % ---> recall SSR = ssr_fitlm = 154.9223, checked

% SST, Sum of total squares: 
SST = (newdat.corruption-mean(newdat.corruption))'*...
    (newdat.corruption-mean(newdat.corruption));
    % ---> recall SST = ssr_fitlm = 615.1223, checked

% SSR, Sum of squares explained(Regression):
yfitted = X*b;
yfitted = yfitted(~any(ismissing(yfitted),2),:);
SSR = (yfitted-mean(newdat.corruption))'*...
    (yfitted-mean(newdat.corruption));
    % ---> recall SSR = ssr_fitlm = 460.2000, checked
% Notice that also, SST = SSE + SSR

%% (a)(ii) Check R^2_{c_1} = R^2_{c_2}

%
% Recall R^2_{c_1} = 1-SSE/SST
Rsqc1 = 1- SSE/SST;
Rsqc2 = SSR/SST;
    % ---> recall R^2_{c_1} = R^2_{c_2}=0.7481, checked
%

%% Question 4(b) Calcualte the measure of leverage, h_i
%
% Initialise the leverage storage
h = zeros(height(Xnew),1);
% below the for loop computes for each obs i, the distance h_i
% and store it in h
for i = 1:height(Xnew)
    h(i) = Xnew(i,:)*(Xnew'*Xnew)*Xnew(i,:)';
end

% largest leverage hunting
[maxValue, maxIndex] = max(h);
highestleverage = dt2{maxIndex,1};

fprintf("The country that has the highest leverage among " + ...
    "observations is %s\n", highestleverage{1});
%

%% Question 4(c) (i) Partition [const lngdp]
%
Xparti = [ones(35,1) newdat.lngdp];
%

%% Question 4(c) (ii) Fwl Regression
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Theorem: FWL                              %
%    X = [const ln(gdp) Socialnetshare]        %
%    X_1 = [const ln(gdp)]                     %
%    X_2 = [Socialnetshare]                    %
%    y = Corruption                            %
%    (R1)     y = Xβ+ε,                        %
%    (R2)     M_1y = M_1βX_2 + ε,              %
%    M_1 = I_n-X(X'X)^{-1}X'                   %
%                                              %
%    Goal: 1. M_1y and M_1X_2                  %
%    2. Regress M_1y on M_1X_2                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Goal 1: M_1y (residual of corruption)
y1 = newdat.corruption;
b1 = Xparti\y1;
Corruption_R = y1-Xparti*b1;

% Goal 2: M_1X_2 (residual of socialnetshare)
y2 = newdat.socialnetshare;
b2 = Xparti\y2;
Socialnetshare_R = y2-Xparti*b2;

% Compute the coefficient
b3 = Socialnetshare_R\Corruption_R;
% Get the fitted values
fwlfitted = Socialnetshare_R*b3;

% Report

if b(3) == b3
    fprintf('The third element of b is equal to b3.\n');
else
    fprintf('The third element of b is not equal to b3.\n');
end
% here, if directly compare the result is not the same
b_string = num2str(b(3), '%.4f');
b3_string = num2str(b3, '%.4f');

if strcmp(b_string, b3_string)
    fprintf('The first four digits of b(3) and b3 are equal.\n');
else
    fprintf('The first four digits of b(3) and b3 are not equal.\n');
end


%% Replicate Plots: Setup

% Before anything, set the graph aesthetics
PS = PLOT_STANDARDS();

% Data points for plotting
dt2 = dt(~any(ismissing(dt),2),:);

% Get Country labels
obs = cellstr(dt2.obs);
short_obs = cell(length(obs),1);
for i = 1:length(obs)
    short_obs{i} = obs{i}(1:3);
end

%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Please use Matlab Add-on             %
%                Labelpoints                     %
%                   and                          %
%             Professional Plots                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf;
figure(1);
fig1_comps.fig = gca;
% turn the grid on but leave only horizontal grid:
grid on;
fig1_comps.fig.XGrid = 'off';
hold on;

fig1_comps.p1 = scatter(Socialnetshare_R, Corruption_R, ...
    'MarkerEdgeColor',PS.Red3,...
    'MarkerFaceColor',PS.Red4);
fig1_comps.p2 = plot(Socialnetshare_R, fwlfitted);

% here mark the points using the index generated in Setup
labelpoints(Socialnetshare_R, Corruption_R, short_obs);

% set aesthetics
set(fig1_comps.p2, 'Color', PS.Red5, ...
    'LineWidth', 3)
yticks(-10:5:5);
xlabel('Share of Social Network Users', ...
    'FontSize',16, 'FontName','Optima');
ylabel('Corruption', ...
    'FontSize',16, 'FontName','Optima');
title('Social Media Penetration and Corruption', ...
    'FontName', 'Optima', 'FontSize',18);
legend('Country', '$y_R = X_R \hat{\beta}$','Location','northeast',...
    'Interpreter', 'latex', ...
    'Fontsize', 12)
text(-3900,-9.5,...
    'coef.=-0.00154626, (robust) standard error = -0.00036924, t = -4.19', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize',12,'FontName','Optima');

hold off; 











