%% Setup
clc;clear;

% Set the working directory to the place where the current file is saved
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));

dt = readtable("Assig1.csv", 'VariableNamingRule','preserve');

%% Question 1 a) estimate of average weekly pay with educ=12

% +++ Recover variable `wage` from `lwklywge`: +++
%
% Recover by taking exponetial of `lwklywge`:
wage = exp(dt.lwklywge);
% Make a table of the recovered wage data:
wage1 = array2table(wage);
% Glue to a new dataset:
dt1 = horzcat(dt, wage1);

% +++ Conditional mean E(wage | educ = 12): +++
%
% Create a Sub-sample with educ == 12:
sample_12 = dt1(dt1.educ == 12,:);
% E(wage | educ = 12):
cmean_wage = mean(sample_12.wage);

%% Question 1 b) plot mean of `lwklywge` conditional on each educ group

% +++ Conditional mean: +++
%
% Get unique educ values
educ_groups = unique(dt1.educ);
% Initialize vector to store means
cmeans = zeros(length(educ_groups),1);
% Compute mean for each educ value
for i = 1:length(educ_groups)
    cmeans(i) = mean(dt1.lwklywge(dt1.educ == educ_groups(i)));
end

% +++ Replicate line chart in Fig 3.1.1: +++
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          get the add on called          %
%           `professional plot`           %
%                 first                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Before anything, set the graph aesthetics
PS = PLOT_STANDARDS();
% Plot the means as a line chart with markers
figure(1);
fig1_comps.fig = gcf;
grid off;
hold on;

% Get the conditional mean for educ == 4, 8, 12, 16
cmeans4 = cmeans(5,:);
cmeans8 = cmeans(9,:);
cmeans12 = cmeans(13,:);
cmeans16 = cmeans(17,:);

% Drawing~
% Main graph:
fig1_comps.p0 = plot(educ_groups, cmeans);
set(fig1_comps.p0, 'Color', PS.Blue5, ...
    'LineWidth', 4, 'Marker','o');
% cmeans for educ == 4, 8, 12, 16:
fig1_comps.p1 = yline(cmeans4, 'Color',PS.Grey3,...
    'LineStyle','--', LineWidth=1.5);
fig1_comps.p2 = yline(cmeans8, 'Color',PS.Red3,...
    'LineStyle','--', LineWidth=1.5);
fig1_comps.p3 = yline(cmeans12, 'Color',PS.Orange3,...
    'LineStyle','--', LineWidth=1.5);
fig1_comps.p4 = yline(cmeans16, 'Color',PS.Green3,...
    'LineStyle','--', LineWidth=1.5);

xlabel('Years of Completed Education', ...
    'FontSize',20, 'FontName','Palatino');
ylabel('Log of Weekly earnings', ...
    'FontSize',20, 'FontName','Palatino');
legend('Conditional means', ...
    '$E(\ln wage \vert educ = 4)$',...
    '$E(\ln wage \vert educ = 8)$',...
    '$E(\ln wage \vert educ = 12)$',...
    '$E(\ln wage \vert educ = 16)$',...
    'Location','best',...
    'Interpreter', 'latex', ...
    'Fontsize', 16)

hold off;

%% Question 1 d) plot mean&med of `lwklywge` conditional on each educ group

% +++ Compute the medians +++:
%
% Initialize vector to store medians
cmeds = zeros(length(educ_groups),1);
% Compute mean for each educ value
for i = 1:length(educ_groups)
    cmeds(i) = median(dt1.lwklywge(dt1.educ == educ_groups(i)));
end
%
% +++ Plot conditional means and medians +++:
%
% Drawing~
figure(2);
fig1_comps.fig = gcf;
grid on;
hold on;

fig2_comps.p0 = plot(educ_groups, cmeans);
fig2_comps.p1 = plot(educ_groups, cmeds);

set(fig2_comps.p0, 'Color', PS.Blue5, ...
    'LineWidth', 4, 'Marker','o');
set(fig2_comps.p1, 'Color', PS.Red5, ...
    'LineWidth', 4, 'Marker','Diamond');

xlabel('Years of Completed Education', ...
    'FontSize',20, 'FontName','Palatino');
ylabel('Log of Weekly earnings', ...
    'FontSize',20, 'FontName','Palatino');
legend('Conditional means', ...
    'Conditional medians',...
    'Location','best',...
    'Interpreter', 'latex', ...
    'Fontsize', 16)

hold off;
%

%% Additional Checkings for Question 1 d)

% +++ Check distribution of `lwklywge`: +++
% 
% Scott's rule for determining optimal bandwidth
bw = 3.5*std(dt1.lwklywge)*length(dt1.lwklywge)^(-1/3);
% Drawing~
figure(3);
fig1_comps.fig = gcf;
grid off;
hold on;
histogram(dt1.lwklywge, 'Normalization', 'pdf', ...
    'BinWidth', bw, 'FaceColor',PS.Grey1);

[density,xi] = ksdensity(dt1.lwklywge);

% Plot density function
plot(xi,density,'LineWidth',3.5,'Color',PS.Red5);

xlabel('Log of Weekly earnings', 'FontSize',22, 'FontName','Palatino');
ylabel('Probability Density', 'FontSize',22, 'FontName','Palatino');
legend('Distribution', 'Smoothed density',...
    'FontSize',18, 'FontName','Palatino','Location', 'Best');
hold off;
%

%% Question 1 e) polynomial fitting a line between Conditional Mean of `lwklywge` and `educ`


% Fit the polynomial of degree one
p = polyfit(educ_groups, cmeans, 1);

% Get the coefficients of the polynomial
coefs = p;

% Generate the fitted values for the independent variable
y_fit = polyval(p,educ_groups);

% Plot the data and the fitted line
figure(4);
fig1_comps.fig = gcf;
grid on;
hold on;

fig4_comps.p1 = plot(educ_groups, y_fit);
fig4_comps.p2 = plot(educ_groups, cmeans);

set(fig4_comps.p1, 'Color', PS.Blue5, ...
    'LineWidth', 3);
set(fig4_comps.p2, 'Color', PS.Red5, ...
    'LineWidth', 3, 'Marker','x');

xlabel('Years of Education', 'FontSize',20, 'FontName','Palatino'); 
ylabel('Log of Weekly earnings', 'FontSize',20, 'FontName','Palatino');
legend('Polynomial Fit of Degree one', 'Conditional Means',...
    'FontSize',18, 'FontName','Palatino','Location', 'Best');

hold off;
%%
