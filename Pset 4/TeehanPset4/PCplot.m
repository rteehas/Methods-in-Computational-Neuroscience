function [H,PCkin,nb,nj] = PCplot(V,mu,n,va,bounds)
%% FUNCTION DESCRIPTION
% [H,PCkin,nb,nj] = PCplot(V,mu,n,va,bounds)
% PCplot takes the output of principal components analysis applied to
% kinematic data and outputs a plot displaying the range of motion
% described by those principal components, using a simple skeletal model.
%
% INPUTS:
% V
% V is a k x q matrix of principal components, where k is the
% dimensionality of the data, and q is the maximum number of principal 
% components you wish to consider. Principal components are assumed to be
% organized such that each column corresponds to a separate PC, whereas
% each row corresponds to a different weight within that PC.
%
% Be sure that you sort your principal components matrix by their
% respective eigenvalues before using it as input for this function!
%
% mu
% mu is the 1 x k (or k x 1) vector depicting the mean kinematic position,
% which is used to center the PC projections
%
% n
% n is a m x 1 (or 1 x m) vector where each element denotes a principal
% component for which you want a kinematic plot.
%
% va
% va is the [azimuth elevation] viewing angle you want to see your movie in
%
% bounds
% bounds is the 2-element vector [min max] that specifies the range spanned
% by the PC projection. Normally, this should be the minimum and maximum
% scores of your PC, but you can instead set this to mu + 2sd or some other
% parametric value, or instead constrain your joint angles arbitrarily
% (e.g., a maximum of 90 ensures a no more than 90-degree deviation from
% mu, as PC vectors are of unit length)
% 
% OUTPUTS:
% H
% H is the set of figure handles used to construct the movie.
%
% PCkin
% PCkin gives a p x k x m matrix of sampled joint angles used to construct
% the kinematic plots. p is the number of samples (taken in this case to be
% 60), k is the dimensionality of the PCs, and m is the mth principal
% component as defined by the input vector n (i.e., if n = [1 3 2], if m =
% 2, it refers to the sampling along the third {n(2) = 3} principal
% component)
%
% nb, nj
% nb and nj are cell arrays wherein each cell contains the bodies (nb) or
% joints (nj) struct for the frame of the video corresponding to the
% similarly-indexed row of PCkin.

%% WARNINGS:
% The axes of rotation along the thumb are slightly off in this version of
% PCplot. The thumb is assumed to operate on an entirely orthogonal plane,
% instead of a more veridically-angled plane, to the rest of the digits. As
% a result, the movements of the thumb could look kinda dumb.
%
% Wrist movements are also assumed to occur along cardinal axes when in
% reality these axes are slightly tilty due to the nuanced geometry of the
% articulations between carpals, radius, and ulna. As a result, wrist
% movements might also be slightly off.
%
% Distal and intermediate phalanges are fused into a single distal segment,
% as in the implementation of this model, distal phalanges are held fixed &
% axial with the intermediate phalanges
%
% Radius and ulna are fused into a single segment. As such, it may be
% difficult to visualize pronation and supination (although the triangle
% representing the carpals should give a good idea of this axis).
%
% Note that the kinematics are unconstrained in this display function, as
% opposed to OpenSim. As such, if the ranges of your PC projections are
% really high, you might see wacky, impossible movements. Furthermore, you
% may see bodies passing through one another.
%
% If you're using this for kinematics homework, the output of this
% script should be reasonably veridical to write a 1-2 sentence qualitative
% blurb about the nature of the movements.

%% Step 0: make things clean (debug only)
% clc, close all

%% Step 1: ensure the data are in the appropriate format

pc_covmat = V'*V;
criterion = trace(pc_covmat) + sum(sum(triu(pc_covmat,1))) + sum(sum(tril(pc_covmat,-1)));

if size(V,1) ~= 22
    error('ERROR: there are 22 joints (degrees of freedom) in the kinematic dataset. Your input does not have 22 dimensions!')
elseif size(V,1) < size(V,2)
    error('ERROR: number of PCs cannot exceed dimensionality of data')
elseif abs(size(V,2) - criterion) > 1e-10
    error('ERROR: input must be a set of orthogonal principal components')
elseif min(size(n)) ~= 1
    error('ERROR: input n must be a 1-dimensional vector (or a single value)')
elseif max(n) > size(V,2)
    error('ERROR: values of n cannot exceed number of PCs')
elseif ~all(round(n) == n)
    error('ERROR: elements of n must be integer values')
elseif numel(va) ~= 2
    error('ERROR: viewing angle must be a 2-element vector [azimuth elevation]')
elseif numel(mu) ~= 22
    error('ERROR: mu must be a vector of the mean position among the 22 degrees of freedom in the model')
elseif numel(bounds) ~= 2
    error('ERROR: bounds must be a 2-element vector')
elseif bounds(1) >= bounds(2)
    error('ERROR: the first element of bounds must be less than the second element')
end


%% Step 2: for each PC desired, take 100 samples centered about the mean (base position)

% substep: define dimensions
k = size(V,1);
q = size(V,2);
m = numel(n);
p = 60;

% define kinematic bounds
min_angle = bounds(1);
max_angle = bounds(2);

PCkin = zeros(p,k,m);
for ii = 1:m % for each PC defined by n
    for jj = 1:k % for each weight in each PC
        
        % define the range of motion specified by the PC weight
        temp_min = mu(jj) + V(jj,n(m))*min_angle;
        temp_max = mu(jj) + V(jj,n(m))*max_angle;
        
        % generate the sampling vector for the range of motion
        PCkin(:,jj,ii) = linspace(temp_min,temp_max,p)';
    end
end



%% Step 3: define parameters for each body
% length:           arbitrary units
% pos:              relative position of proximal end, in base position
% each body is assumed to be a stick. This doesn't really work for the
% wrist, but we'll roll with it.
% ^  y
% |      ^ z
% |     /
% |    /
% |   /
% |  /
% | /
% ---------------------->  x   (sagittal view)
%
% Note also that the middle phalangeal and distal phalangeal segments are
% lumped under "distal", as in our kinematics, the distal phalangeal
% segment is permanently constrained to 0 and therefore always points the
% same direction as the middle phalangeal segment.
%
% these are also kinematics from the macaque, so if the proportions of the
% digits look a little off, that's why
%
% RIGHT HAND RULE assuming +x = distal, +y = superior, +z = right


bodies.humerus.length    = 8;
bodies.humerus.pos       = [0 0 0];

bodies.rad_ulna.length   = 6;
bodies.rad_ulna.pos      = bodies.humerus.pos + [0 -bodies.humerus.length 0];

bodies.wrist.length      = 1;
bodies.wrist.pos         = bodies.rad_ulna.pos + [-bodies.rad_ulna.length 0.4 0];

% make this a triangle to represent the plane, rather than suggesting that
% it's just a single linear segment
bodies.MC1.length        = 1;
bodies.MC1.pos           = bodies.wrist.pos + [-bodies.wrist.length 0.4 0];

bodies.MC2.length        = 2;
bodies.MC2.pos           = bodies.wrist.pos + [-bodies.wrist.length 0.2 0];

bodies.MC3.length        = 2;
bodies.MC3.pos           = bodies.wrist.pos + [-bodies.wrist.length 0 0];

bodies.MC4.length        = 2;
bodies.MC4.pos           = bodies.wrist.pos + [-bodies.wrist.length -0.2 0];

bodies.MC5.length        = 1.8;
bodies.MC5.pos           = bodies.wrist.pos + [-bodies.wrist.length -0.4 0];

bodies.D1_prox.length    = 0.5;
bodies.D1_prox.pos       = bodies.MC1.pos + bodies.MC1.length*[-cosd(30) 0 -sind(30)];

bodies.D2_prox.length    = 1;
bodies.D2_prox.pos       = bodies.MC2.pos + bodies.MC2.length*[-cosd(5) sind(5) 0];

bodies.D3_prox.length    = 1.2;
bodies.D3_prox.pos       = bodies.MC3.pos + bodies.MC3.length*[-1 0 0];

bodies.D4_prox.length    = 1.1;
bodies.D4_prox.pos       = bodies.MC4.pos + bodies.MC4.length*[-cosd(-5) sind(-5) 0];

bodies.D5_prox.length    = 0.8;
bodies.D5_prox.pos       = bodies.MC5.pos + bodies.MC5.length*[-cosd(-10) sind(-10) 0];

bodies.D1_dist.length    = 0.4;
bodies.D1_dist.pos       = bodies.D1_prox.pos + bodies.D1_dist.length*[-cosd(30) 0 -sind(30)];

bodies.D2_dist.length    = 0.5;
bodies.D2_dist.pos       = bodies.D2_prox.pos + bodies.D2_prox.length*[-cosd(5) sind(5) 0];

bodies.D3_dist.length    = 1.1;
bodies.D3_dist.pos       = bodies.D3_prox.pos + bodies.D3_prox.length*[-1 0 0];

bodies.D4_dist.length    = 1.1;
bodies.D4_dist.pos       = bodies.D4_prox.pos + bodies.D4_prox.length*[-cosd(-5) sind(-5) 0];

bodies.D5_dist.length    = 1;
bodies.D5_dist.pos       = bodies.D5_prox.pos + bodies.D5_prox.length*[-cosd(-10) sind(-10) 0];

% these are endpoints for the distal segments; they have length 0 for this
% reason, as they are not true segments
bodies.D1_end.length    = 0;
bodies.D1_end.pos       = bodies.D1_dist.pos + bodies.D1_dist.length*[-cosd(30) 0 -sind(30)];

bodies.D2_end.length    = 0;
bodies.D2_end.pos       = bodies.D2_dist.pos + bodies.D2_dist.length*[-cosd(5) sind(5) 0];

bodies.D3_end.length    = 0;
bodies.D3_end.pos       = bodies.D3_dist.pos + bodies.D3_dist.length*[-1 0 0];

bodies.D4_end.length    = 0;
bodies.D4_end.pos       = bodies.D4_dist.pos + bodies.D4_dist.length*[-cosd(-5) sind(-5) 0];

bodies.D5_end.length    = 0;
bodies.D5_end.pos       = bodies.D5_dist.pos + bodies.D5_dist.length*[-cosd(-10) sind(-10) 0];
%% Step 4: define your joints
% define locations & axes of rotation & all distal & coincident bodies & joints
% (all will be based on the "base" state of the arm model)

% wrist supination/pronation
joints.wr_sup_pro_l.pos   = bodies.wrist.pos;
joints.wr_sup_pro_l.u     = -[1 0 0];
joints.wr_sup_pro_l.bods  = {'wrist','MC1','MC2','MC3','MC4','MC5','D1_prox',...
    'D2_prox','D3_prox','D4_prox','D5_prox','D1_dist','D2_dist','D3_dist',...
    'D4_dist','D5_dist','D1_end','D2_end','D3_end','D4_end','D5_end'};
joints.wr_sup_pro_l.joints = {'wr_sup_pro_l','wr_adduction_l','wr_flexion_l','CMC1_flex_l',...
    'CMC1_abd_l','CMC4_ext','CMC5_ext','CMC5_abd','mcp1_flexion_l',...
    'mcp1_abduction_l','mcp2_flexion_l','mcp2_abduction_l',...
    'pm2_flexion_l','mcp3_flexion_l','mcp3_abduction_l','pm3_flexion_l',...
    'mcp4_flexion_l','mcp4_adduction_l','pm4_flexion_l','mcp5_flexion_l',...
    'mcp5_adduction_l','pm5_flexion_l'};

% wrist add/abduction
joints.wr_adduction_l.pos   = bodies.wrist.pos;
joints.wr_adduction_l.u     = [0 0 1];
joints.wr_adduction_l.bods  = {'wrist','MC1','MC2','MC3','MC4','MC5','D1_prox',...
    'D2_prox','D3_prox','D4_prox','D5_prox','D1_dist','D2_dist','D3_dist',...
    'D4_dist','D5_dist','D1_end','D2_end','D3_end','D4_end','D5_end'};
joints.wr_adduction_l.joints = {'wr_sup_pro_l','wr_adduction_l','wr_flexion_l','CMC1_flex_l',...
    'CMC1_abd_l','CMC4_ext','CMC5_ext','CMC5_abd','mcp1_flexion_l',...
    'mcp1_abduction_l','mcp2_flexion_l','mcp2_abduction_l',...
    'pm2_flexion_l','mcp3_flexion_l','mcp3_abduction_l','pm3_flexion_l',...
    'mcp4_flexion_l','mcp4_adduction_l','pm4_flexion_l','mcp5_flexion_l',...
    'mcp5_adduction_l','pm5_flexion_l'};


% wrist flex/extension
joints.wr_flexion_l.pos   = bodies.wrist.pos;
joints.wr_flexion_l.u     = [0 -1 0];
joints.wr_flexion_l.bods  = {'wrist','MC1','MC2','MC3','MC4','MC5','D1_prox',...
    'D2_prox','D3_prox','D4_prox','D5_prox','D1_dist','D2_dist','D3_dist',...
    'D4_dist','D5_dist','D1_end','D2_end','D3_end','D4_end','D5_end'};
joints.wr_flexion_l.joints = {'wr_sup_pro_l','wr_adduction_l','wr_flexion_l','CMC1_flex_l',...
    'CMC1_abd_l','CMC4_ext','CMC5_ext','CMC5_abd','mcp1_flexion_l',...
    'mcp1_abduction_l','mcp2_flexion_l','mcp2_abduction_l',...
    'pm2_flexion_l','mcp3_flexion_l','mcp3_abduction_l','pm3_flexion_l',...
    'mcp4_flexion_l','mcp4_adduction_l','pm4_flexion_l','mcp5_flexion_l',...
    'mcp5_adduction_l','pm5_flexion_l'};


% CMCs
% CMC1
joints.CMC1_flex_l.pos   = bodies.MC1.pos;
joints.CMC1_flex_l.u     = [-sind(30) 0 cosd(30)];
joints.CMC1_flex_l.bods  = {'MC1','D1_prox','D1_dist','D1_end'};
joints.CMC1_flex_l.joints = {'CMC1_flex_l','CMC1_abd_l',...
    'mcp1_flexion_l','mcp1_abduction_l'};

joints.CMC1_abd_l.pos   = bodies.MC1.pos;
joints.CMC1_abd_l.u     = [0 -1 0]; %-[1/3 3 1]/norm([3 1/3 1]);
joints.CMC1_abd_l.bods  = {'MC1','D1_prox','D1_dist','D1_end'};
joints.CMC1_abd_l.joints = {'CMC1_flex_l','CMC1_abd_l',...
    'mcp1_flexion_l','mcp1_abduction_l'};

% CMC4
joints.CMC4_ext.pos   = bodies.MC4.pos;
joints.CMC4_ext.u     = [sind(-5) cosd(-5) 0];
joints.CMC4_ext.bods  = {'MC4','D4_prox','D4_dist','D4_end'};
joints.CMC4_ext.joints = {'CMC4_ext','mcp4_flexion_l',...
    'mcp4_adduction_l','pm4_flexion_l'};

% CMC5
joints.CMC5_ext.pos   = bodies.MC5.pos;
joints.CMC5_ext.u     = [sind(-10) cosd(-10) 0];
joints.CMC5_ext.bods  = {'MC5','D5_prox','D5_dist','D5_end'};
joints.CMC5_ext.joints = {'CMC5_ext','CMC5_abd','mcp5_flexion_l',...
    'mcp5_adduction_l','pm5_flexion_l'};

joints.CMC5_abd.pos   = bodies.MC5.pos;
joints.CMC5_abd.u     = -[0 0 1];
joints.CMC5_abd.bods  = {'MC5','D5_prox','D5_dist','D5_end'};
joints.CMC5_abd.joints = {'CMC5_ext','CMC5_abd','mcp5_flexion_l',...
    'mcp5_adduction_l','pm5_flexion_l'};


% MCPs
% D1
joints.mcp1_flexion_l.pos   = bodies.D1_prox.pos;
joints.mcp1_flexion_l.u     = [-sind(30) 0 cosd(30)];
joints.mcp1_flexion_l.bods  = {'D1_prox','D1_dist','D1_end'};
joints.mcp1_flexion_l.joints = {'mcp1_flexion_l','mcp1_abduction_l'};

joints.mcp1_abduction_l.pos   = bodies.D1_prox.pos;
joints.mcp1_abduction_l.u     = [0 1 0];
joints.mcp1_abduction_l.bods  = {'D1_prox','D1_dist','D1_end'};
joints.mcp1_abduction_l.joints = {'mcp1_flexion_l','mcp1_abduction_l'};

% D2
joints.mcp2_flexion_l.pos   = bodies.D2_prox.pos;
joints.mcp2_flexion_l.u     = [sind(5) -cosd(5) 0];
joints.mcp2_flexion_l.bods  = {'D2_prox','D2_dist','D2_end'};
joints.mcp2_flexion_l.joints = {'mcp2_flexion_l','mcp2_abduction_l','pm2_flexion_l'};

joints.mcp2_abduction_l.pos   = bodies.D2_prox.pos;
joints.mcp2_abduction_l.u     = -[0 0 1];
joints.mcp2_abduction_l.bods  = {'D2_prox','D2_dist','D2_end'};
joints.mcp2_abduction_l.joints = {'mcp2_flexion_l','mcp2_abduction_l','pm2_flexion_l'};

% D3
joints.mcp3_flexion_l.pos   = bodies.D3_prox.pos;
joints.mcp3_flexion_l.u     = [0 -1 0];
joints.mcp3_flexion_l.bods  = {'D3_prox','D3_dist','D3_end'};
joints.mcp3_flexion_l.joints = {'mcp3_flexion_l','mcp3_abduction_l','pm3_flexion_l'};

joints.mcp3_abduction_l.pos   = bodies.D3_prox.pos;
joints.mcp3_abduction_l.u     = -[0 0 1];
joints.mcp3_abduction_l.bods  = {'D3_prox','D3_dist','D3_end'};
joints.mcp3_abduction_l.joints = {'mcp3_flexion_l','mcp3_abduction_l','pm3_flexion_l'};

% D4
joints.mcp4_flexion_l.pos   = bodies.D4_prox.pos;
joints.mcp4_flexion_l.u     = [sind(-5) -cosd(-5) 0];
joints.mcp4_flexion_l.bods  = {'D4_prox','D4_dist','D4_end'};
joints.mcp4_flexion_l.joints = {'mcp4_flexion_l','mcp4_adduction_l','pm4_flexion_l'};

joints.mcp4_adduction_l.pos   = bodies.D4_prox.pos;
joints.mcp4_adduction_l.u     = -[0 0 1];
joints.mcp4_adduction_l.bods  = {'D4_prox','D4_dist','D4_end'};
joints.mcp4_adduction_l.joints = {'mcp4_flexion_l','mcp4_adduction_l','pm4_flexion_l'};

% D5
joints.mcp5_flexion_l.pos   = bodies.D5_prox.pos;
joints.mcp5_flexion_l.u     = [sind(-10) -cosd(-10) 0];
joints.mcp5_flexion_l.bods  = {'D5_prox','D5_dist','D5_end'};
joints.mcp5_flexion_l.joints = {'mcp5_flexion_l','mcp5_adduction_l','pm5_flexion_l'};

joints.mcp5_adduction_l.pos   = bodies.D5_prox.pos;
joints.mcp5_adduction_l.u     = -[0 0 1];
joints.mcp5_adduction_l.bods  = {'D5_prox','D5_dist','D5_end'};
joints.mcp5_adduction_l.joints = {'mcp5_flexion_l','mcp5_adduction_l','pm5_flexion_l'};


% PIPs
% D2
joints.pm2_flexion_l.pos   = bodies.D2_dist.pos;
joints.pm2_flexion_l.u     = [sind(5) -cosd(5) 0];
joints.pm2_flexion_l.bods  = {'D2_dist','D2_end'};
joints.pm2_flexion_l.joints = {'pm2_flexion_l'};

% D3
joints.pm3_flexion_l.pos   = bodies.D3_dist.pos;
joints.pm3_flexion_l.u     = [0 -1 0];
joints.pm3_flexion_l.bods  = {'D3_dist','D3_end'};
joints.pm3_flexion_l.joints = {'pm3_flexion_l'};

% D4
joints.pm4_flexion_l.pos   = bodies.D4_dist.pos;
joints.pm4_flexion_l.u     = [sind(-5) -cosd(-5) 0];
joints.pm4_flexion_l.bods  = {'D4_dist','D4_end'};
joints.pm4_flexion_l.joints = {'pm4_flexion_l'};

% D5
joints.pm5_flexion_l.pos   = bodies.D5_dist.pos;
joints.pm5_flexion_l.u     = [sind(-10) -cosd(-10) 0];
joints.pm5_flexion_l.bods  = {'D5_dist','D5_end'};
joints.pm5_flexion_l.joints = {'pm5_flexion_l'};


%% DEBUGGING 
% function left in to help students manipulate the image to get a good
% viewing angle for their plot

[newbods ~] = jointrot(bodies,joints,mu);
plotarm(newbods,[0 90]); title('DEBUG FIGURE: base position (\mu)')

%% Step 5: start making plots!
% debug: just plot the first one for now
figure('paperposition',[0 0 4 3])
nb = cell(p,1);
nj = cell(p,1);
for ii = 1:p
    for jj = 1:size(PCkin,3)
        subplot(1,size(PCkin,3),jj)
        [nb{ii} nj{ii}] = jointrot(bodies,joints,PCkin(ii,:,jj));
        plotarm(nb{ii},va);
    end
    G(ii) = getframe(gcf);
    
    cla
end

%% Step 6: poop out a movie
close(gcf)
H = frame2im(G(1));
H = zeros([size(H),numel(G)-1]);


for ii = 2:numel(G)
    H(:,:,:,ii-1) = frame2im(G(ii));
end    

implay(H,30);

%% define functions


% small, low-level crud
    function [X Y Z] = XYZsep(XYZ)
        % take xyz coordinates and split them into separate arrays
        X = XYZ(1);
        Y = XYZ(2);
        Z = XYZ(3);
        % weird order is to accomodate MATLAB's plot3 convention
    end

    function h = plotseg(startpt,endpt)
        % plot a line based on the beginning and ending coordinates
        [Xp Yp Zp] = XYZsep(startpt);
        [Xd Yd Zd] = XYZsep(endpt);
        plot3([Xp Xd],[Yp Yd],[Zp Zd],'kx','markersize',12)
        hold all
        h = line([Xp Xd],[Yp Yd],[Zp Zd]);
    end


% arm plotting function
    function h = plotarm(bodies,va)
%         digit_colors = hsv(6);
        digit_colors = [1   0   0;
                        1   0   1;
                        1   0   1;
                        1   0   1;
                        1   0   1;
                        0   0   1;];
        % humerus
        humerus = plotseg(bodies.humerus.pos,bodies.rad_ulna.pos);
        set(humerus,'color','k')
        set(humerus,'linewidth',5)
        hold all
        
        % rad_ulna
        rad_ulna = plotseg(bodies.rad_ulna.pos,bodies.wrist.pos);
        set(rad_ulna,'color','k')
        set(rad_ulna,'linewidth',5)
        hold all
        
        % wrist
        carp_ulnar = plotseg(bodies.wrist.pos,bodies.MC5.pos);
        set(carp_ulnar,'color',digit_colors(6,:))
        set(carp_ulnar,'linewidth',5)
        hold all
        
        carp_radial = plotseg(bodies.wrist.pos,bodies.MC1.pos);
        set(carp_radial,'color',digit_colors(6,:))
        set(carp_radial,'linewidth',5)
        hold all
        
        carp_flat = plotseg(bodies.MC5.pos,bodies.MC1.pos);
        set(carp_flat,'color',digit_colors(6,:))
        set(carp_flat,'linewidth',5)
        hold all
        
        MC1 = 0;
        MC2 = 0;
        MC3 = 0;
        MC4 = 0;
        MC5 = 0;
        
        D1_prox = 0;
        D2_prox = 0;
        D3_prox = 0;
        D4_prox = 0;
        D5_prox = 0;
        
        D1_dist = 0;
        D2_dist = 0;
        D3_dist = 0;
        D4_dist = 0;
        D5_dist = 0;
        
        % metacarpals, proximal phalanges, & distal phalanges
        for i = 5:-1:1
            eval(['MC',num2str(i),' = plotseg(bodies.MC',num2str(i),'.pos,bodies.D',num2str(i),'_prox.pos);']);
            eval(['set(MC',num2str(i),',''color'',0.6*digit_colors(i,:))']);
            eval(['set(MC',num2str(i),',''linewidth'',5)']);
            hold all

            eval(['D',num2str(i),'_prox = plotseg(bodies.D',num2str(i),'_prox.pos,bodies.D',num2str(i),'_dist.pos);']);
            eval(['set(D',num2str(i),'_prox,''color'',0.8*digit_colors(i,:))']);
            eval(['set(D',num2str(i),'_prox,''linewidth'',5)']);
            hold all

            
            eval(['D',num2str(i),'_dist = plotseg(bodies.D',num2str(i),'_dist.pos,bodies.D',num2str(i),'_end.pos);']);
            eval(['set(D',num2str(i),'_dist,''color'',digit_colors(i,:))']);
            eval(['set(D',num2str(i),'_dist,''linewidth'',5)']);
            hold all
        end
        
%         g = legend('thumb','index','middle','ring','little','carpals');
%         legend boxoff
%         gg = get(g,'children');
        
        
%         for i = 1:6
%             set(gg(3*(i-1)+1),'linestyle','none');
%             set(gg(3*(i-1)+1),'marker','none');
%             set(gg(3*(i-1)+2),'linestyle','none');
%             set(gg(3*(i-1)+2),'marker','none');
%             set(gg(3*(i-1)+3),'fontweight','demi')
%             set(gg(3*(i-1)+3),'fontsize',14)
%             set(gg(3*(i-1)+3),'color',digit_colors(6-(i-1),:));
%         end
        
        
        view(va)
        box off
        xlabel('anterior (-) / posterior (+), a.u. (X)','fontsize',14)
        ylabel('inferior (-) / superior (+), a.u. (Y)','fontsize',14)
        zlabel('right (-) / left (+), a.u. (Z)','fontsize',14)

        
        axis([-15 5 -15 5 -15 5])
        grid on
        

        % set objects to front/back based on viewing angle (VERY BASIC!)
        
        % establish your viewing angle; projections onto this vector
        % determine whether objects are in front/behind
        
        % default viewing angle (0,0) is a vector defined from below:
        u0 = [0 1 0];
        
        % azimuth rotates points about the Z axis (clockwise, hence the
        % negative sign)
        u_az = [0 0 -1];
        
        % elevation is dependent on azimuth (multiplication is not
        % commutative!), but its default is to rotate points
        % along the X axis (clockwise, hence the negative sign)
        u_el = [-1 0 0];
        
        
        % define rotation matrix from base to azimuth
        R = rotmat(u_az,va(1));
        
        % define new u0 and u_el
        u0   = u0*R';
        u_el = u_el*R';
        
        % define rotation matrix from azimuth to azimuth + elevation
        R = rotmat(u_el,va(2));
        
        % update your viewing angle
        u0 = u0*R';
        
        
        % get coordinates for mean point of each body
        names = {'humerus','rad_ulna','carp_ulnar','carp_radial','carp_flat',...
            'MC1','MC2','MC3','MC4','MC5','D1_prox','D2_prox','D3_prox',...
            'D4_prox','D5_prox','D1_dist','D2_dist','D3_dist',...
            'D4_dist','D5_dist'};
        
        mX = zeros(size(names));
        mY = zeros(size(names));
        mZ = zeros(size(names));
        
        for i = 1:numel(names)
            eval(['mX(i) = mean(get(',names{i},',''Xdata''));']);
            eval(['mY(i) = mean(get(',names{i},',''Ydata''));']);
            eval(['mZ(i) = mean(get(',names{i},',''Zdata''));']);
        end
        
        % now obtain the scalar projection onto u0 (dot products)
        mu0 = u0*[mX;mY;mZ];
        
        % sort projections in descending order; larger values are further
        % away from you than smaller values are
        [~,inds] = sort(mu0,'descend');
        
        
        % NOW we get to order the display of things:
        for i = 1:numel(names)
            eval(['uistack(',names{inds(i)},', ''top'')']);
        end
        
        % this might get ugly, especially for longer segments that might be
        % aligned non-normal w.r.t. our viewing angle. But to do that
        % well, I'd need to split all the segments into smaller segments,
        % or do something clever that eludes me at the moment.
        %
        % Just as I thought, it's kinda flickery. Oh well I tried, at least
        % it's better than before, where you'd consistently have, say, a
        % thumb behind the rest of the hand coming to the fore for whatever
        % reason
        
        h = gcf;
        
    end

% rotation matrix from axis and angle
% theta will be in degrees
    function R = rotmat(u,theta)
        R = zeros(3);
        u = u/norm(u);
        
        for i = 1:3
            R(i,i) = cosd(theta) + u(i)^2 * (1-cosd(theta));
        end
        
        R(2,1) = u(2)*u(1)*(1-cosd(theta)) + u(3)*sind(theta);
        R(3,1) = u(3)*u(1)*(1-cosd(theta)) - u(2)*sind(theta);
        R(3,2) = u(3)*u(2)*(1-cosd(theta)) + u(1)*sind(theta);
        
        R(1,2) = u(2)*u(1)*(1-cosd(theta)) - u(3)*sind(theta);
        R(1,3) = u(3)*u(1)*(1-cosd(theta)) + u(2)*sind(theta);
        R(2,3) = u(3)*u(2)*(1-cosd(theta)) - u(1)*sind(theta);
    end

% joint rotation function
% update joints & bodies structs based on a 1D vector of joint rotations

    function [bodies,joints] = jointrot(bodies,joints,jr)
        jn = {'wr_sup_pro_l','wr_adduction_l','wr_flexion_l','CMC1_flex_l',...
            'CMC1_abd_l','CMC4_ext','CMC5_ext','CMC5_abd','mcp1_flexion_l',...
            'mcp1_abduction_l','mcp2_flexion_l','mcp2_abduction_l',...
            'pm2_flexion_l','mcp3_flexion_l','mcp3_abduction_l','pm3_flexion_l',...
            'mcp4_flexion_l','mcp4_adduction_l','pm4_flexion_l','mcp5_flexion_l',...
            'mcp5_adduction_l','pm5_flexion_l'};
        
        u = 0;
        pos = 0;
        bods = 0;
        jts = 0;
        c = 0;
        newpts = 0;
        
        for i = 1:numel(jr)
            eval(['u = joints.',jn{i},'.u;']);
            eval(['pos = joints.',jn{i},'.pos;']);
            eval(['bods = joints.',jn{i},'.bods;']);
            eval(['jts = joints.',jn{i},'.joints;']);
            R = rotmat(u,jr(i));
            
            for j = 1:numel(bods)
                eval(['c = bodies.',bods{j},'.pos;']);
                eval(['bodies.',bods{j},'.pos = pos + (c-pos)*R'';']);
            end
            
            for j = 1:numel(jts)
                eval(['c = joints.',jts{j},'.pos;']);
                eval(['u = joints.',jts{j},'.u;']);
                newpts = ([c;c+u]-[pos;pos])*R';
                newpts = newpts + [pos;pos];
                newpts(2,:) = newpts(2,:) - newpts(1,:);
                eval(['joints.',jts{j},'.pos = newpts(1,:);']);
                eval(['joints.',jts{j},'.u = newpts(2,:);']);
            end 
        end
    end








end