load('hipp_data.mat');
%% plot trajectory and spikes

figure;
plot(xN,yN);
hold on;
% get positions at spikes
[x_1, y_1] = selectTime(xN,yN, spikes);
[x_2, y_2] = selectTime(xN,yN, spikes2);
scatter(x_1,y_1,'red')
hold on;
scatter(x_2,y_2,'k')

fprintf("The spikes for both seem to be concentrated more towards the")
fprintf("negative spatial region (ex. values close to x = -.2 y=-.5 to -1)")

%% 2
% run the script
hipp_glm

sig = sum(stats.p<.05);
aic_lin = dev + 2*2;
% all are significant
fprintf("All are significant")
fprintf("This generally captures the spatial firing properties, but does not")
fprintf("capture the degree of shift befween low firing areas and high firing areas")
fprintf("and the low area around the high")

% check the fit of the graph and write it down 

%%
% quadratic = xN,yN,xN^2,yN^2, xNyN
[b,dev,stats] = glmfit([xN yN xN.*yN xN.^2 yN.^2],spikes,'poisson');

% put all of the below into a function
aic_quad = dev + 2*5;

% comput AIC
fprintf("The AIC value for the linear model is %1.8f\n", aic_lin)
fprintf("The AIC value for the quadratic model is %1.8f", aic_quad)

%% plot
% doesn't make sense to write a new function to handle both graphs 
figure;
[x_new,y_new]=meshgrid(-1:.1:1);

% compute lambda for each point on this grid using the GLM model
lambda = exp( b(1) + b(2)*x_new + b(3)*y_new + b(4)*x_new*y_new + b(5)*x_new.^2 + b(6)*y_new.^2);
lambda(find(x_new.^2+y_new.^2>1))=nan;

%plot lambda as a function of position over this grid
h_mesh = mesh(x_new,y_new,lambda,'AlphaData',0);
hold on;
plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));

fprintf("The quadratic model describes the data better")
%% Spikes Hist 

model = [xN yN xN.*yN xN.^2 yN.^2];
[b, dev, stats, aics, new_model] = addHistory(model, spikes_hist, spikes);
figure;
plot(aics)
xlabel("Value of p")
ylabel("AIC")


%% add other history

[new_b, new_dev, new_stats, new_aics, add_model] = addHistory(new_model, spikes2_hist, spikes);

num_sig = length(new_stats.p(new_stats.p(26:46) <= .05));

% set significance as .05
fprintf("There are %2f parameters that are significant", num_sig)

% the model with all spike_hist and one spike2_hist is most parsimonious 
fprintf("The model will all the spikes_hist and one value for spikes2_hist")
fprintf("is the most parsimonious because the AIC does not drop significantly")
fprintf("when adding the other spikes2_hist but the number of parameters increases")

%%


m = [spikes_hist(:,1:20) spikes2_hist(:,1:20)];


[final_b final_dev final_stats] = glmfit(m, spikes, 'poisson');


sig = sum(final_stats.p(21:41) <= .05);

fprintf("There are %2f network interaction parameters that are significant", sig)

fprintf("Theinteraction terms may be significant in this case because previously")
fprintf("the spatial location may have been triggering activity in the other neuron")
fprintf("which may have lead to the response in the modelled neuron.\n")
fprintf("As a result, the newtwork interactions would not be significant unless")
fprintf("the location was not included")











