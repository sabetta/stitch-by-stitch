
fnBase = "C:\\Users\\Sam\\Documents\\GitHub\\uniaxial\\ignore\\Rib_x_pearl\\"; % Folder containing video1.mp4 through video5.mp4
% description = fileread(fnBase + "description.txt");

forces = {};
colors = {'#ac92eb', '#4fc1e8', '#a0d568', '#ffce54','#ed5564','#000000'};
for l = 1:5
    fnHead = sprintf("video%d",l);

    obj = load(sprintf("%sforce_%s_orig.mat",fnBase,fnHead));
    forces{l} = obj.forces;
    
    while true
        cla
        scatter(1:numel(forces{l}),forces{l})

        rect = drawrectangle;
        minBound = min(rect.Vertices,[],1);
        maxBound = max(rect.Vertices,[],1);

        if maxBound == minBound
            break
        end
        
        ixs = 1:numel(forces{l});
        
        killIx = (ixs < maxBound(1)) & (ixs > minBound(1)) & (forces{l} >= minBound(2)) & (forces{l} < maxBound(2));
        
        forces{l}(killIx) = NaN;
    end
    
    obj.forces = forces{l};
    save(sprintf("%sforce_%s_modify.mat",fnBase,fnHead),'-struct','obj');
end