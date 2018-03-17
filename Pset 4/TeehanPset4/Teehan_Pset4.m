
load('Kinematics.mat')

%% 1

% get a big matrix for PCA
big = [];
for i=1:length(Kinematics.Trials)
    big = vertcat(big, Kinematics.Trials{i});
end

%% PCA

[coeff,score,latent,tsquared,explained,mu] = pca(big);

%% 2

% normalize
latent = latent / sum(latent);

figure; hold on
plot(latent)
title("Normalized Eigenvalues")

fprintf("The first PC explains %2.4f percent of the variance\n", explained(1))
fprintf("The first two PCs explain %2.4f percent of the variance\n", explained(1) + explained(2)) 

total_explained = 0;
num_90 = 0;
for i=1:length(latent)
    total_explained = total_explained + explained(i);
    if total_explained >= 90
        num_90 = i;
        break
    end
end


fprintf("We need %1f PCs to explain at least 90 percent of the data\n", num_90)


%% 3
% I couldn't get it to publish with the PCplot, so uncomment this to use
% the PCplot for #3
%[~] = PCplot(coeff, mean(big), 1, [0 90], [min(min(score)) max(max(score))]);
%[~] = PCplot(coeff, mean(big), 2, [0 90], [min(min(score)) max(max(score))]);


%% 4
% PC * score for 4, not other way around
% first PC reconstruction 

s = size(Kinematics.Trials{1});
% take only the 3rd row for wr_flexion_l
trial_1 = score(1:s(1),:);

% reconstruct, PCA * trial transposed
reconstructed_1 = trial_1 * (coeff(:,1) * coeff(:,1)') + mu;
reconstructed_2 = trial_1 * (coeff(:,1:2) * coeff(:,1:2)') + mu;
reconstructed_3 = trial_1 * (coeff(:,1:3) * coeff(:,1:3)') + mu;
full_recon = trial_1 * coeff';

figure; hold on
subplot(10,1,[1 2])
plot(Kinematics.Trials{1}(3,:))
title("Measured")

subplot(10,1,[3 4])
plot(reconstructed_1(:,3))
title("1 PC")

subplot(10,1,[5 6])
plot(reconstructed_2(:,3))
title("2 PCs")

subplot(10,1,[7 8])
plot(reconstructed_3(:,3))
title("3 PCs")

subplot(10,1,[9 10])
plot(full_recon(:,3))
title("Full Reconstruction")

%% 5 

% real data * pc 
% each row in real data is a time point 

% demean the data 
demeaned_trial_1 = big(1:s(1),:) - mu;

% project onto first pc
proj_1 = demeaned_trial_1 * coeff(:,1);

% project onto second pc 
proj_2 = demeaned_trial_1 * coeff(:,2);

% project onto 3rd pc
proj_3 = demeaned_trial_1 * coeff(:,3);

figure; hold on
subplot(1,2,1)
plot(proj_1, 'r')
hold on
plot(proj_2, 'b')
hold on
plot(proj_3, 'g')
title("Projections Onto Each of the First 3 PCs")
sz = 10;
% this next part was not completely clear, but I assumed it meant plot the
% trajectory of those projections within the 3d space 
subplot(1,2,2)
l = size(proj_1, 1);
cdata = [1:l NaN];
p = patch([proj_1' NaN],[proj_2' NaN], [proj_3' NaN], 0);
set(p,'cdata', cdata, 'edgecolor','interp','facecolor','none')
view(3)
title('Wrist Trajectory PC #1, #2, and #3', 'fontsize', sz);
xlabel('PC 1', 'fontsize', sz);
ylabel('PC 2', 'fontsize', sz);
zlabel('PC 3', 'fontsize', sz);


%% analysis 

figure; hold on
for i=15:20
    hold on
    plot(demeaned_trial_1 * coeff(:,i));
end

fprintf("The amplitudes and frequencies tended to increase as the eigenvalues")
fprintf("decrease.")

