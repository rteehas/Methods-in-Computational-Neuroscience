
% Creates a colormap from an array of directions that associates each
% direction with a color for data associated with it when plotted 

function colormap = GetColormap(directions)
% get the unique directions

uq_dir = unique(directions);


% create a colormap
% to let the keys be doubles, we need to create a nonempty map with 
% double key

colormap = containers.Map(double(5), [0 0 0]);
colormap.remove(5);
rgb_range = linspace(0,1,16);
for i=1:length(uq_dir)
    if ~colormap.isKey(uq_dir(i))
        mod_val = mod(i,3) + 1;
        base = [0 0 0];
        base(mod_val) = rgb_range(i);
        colormap(uq_dir(i)) = base;
    end
end

end