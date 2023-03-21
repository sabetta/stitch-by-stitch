
fnBase = "C:\\Users\\Sam\\Documents\\GitHub\\uniaxial\\ignore\\Rib_y_pearl\\"; % Folder containing video1.mp4 through video5.mp4
description = fileread(fnBase + "description.txt");
doInterp = true;
pinMode = 2;
requireQuasistatic = false;


figure(1)
clf
figure(2)
clf
vertStretch = {};
horStretch = {};
traj = {};
forces = {};
colors = {'#ac92eb', '#4fc1e8', '#a0d568', '#ffce54','#ed5564','#000000'};
for l = 1:5
    fnHead = sprintf("video%d",l);
%     fnHead = "Seed_y";
%     description = "Baby " + fnHead;

    obj = load(sprintf("%sforce_%s.mat",fnBase,fnHead));
    forces{l} = obj.forces;
    
%     forces{l}(forces{l} > 30) = NaN;
%     firstPart = forces{l}(1:1200);
%     firstPart(firstPart > 7) = NaN;
%     forces{l}(1:1200) = firstPart;

    if isfield(obj,'frameReduceFactor')
        frameReduceFactor = obj.frameReduceFactor;
    else
        frameReduceFactor = 1;
    end
    obj = load(sprintf(fnBase + "measurements_video1.mat",l));
%     obj = load(sprintf(fnBase + "measurements_%s.mat",fnHead));
    clampSize = obj.clampSize;
    fabricSize = obj.fabricSize;
    tracks = readmatrix(sprintf(fnBase + "%sautotracks.csv",fnHead));

    [~,I] = sort(tracks(:,8));

    trackId = tracks(I,3);
    x = tracks(I,5);
    y = tracks(I,6);
    t = frameReduceFactor .* tracks(I,8);

    track = unique(trackId);
    
    if numel(track) ~= pinMode && pinMode ~= 4 && pinMode ~= 2
        error("Wrong pinMode!")
    end

    figure(1)
%     clf
    hold on
    trajs = nan(numel(forces{l}),2,numel(track));
    for k = 1:numel(track)
        ts = 1 + t(trackId == track(k));
        trajs(ts,1,k) = x(trackId == track(k));
        trajs(ts,2,k) = y(trackId == track(k));

        scatter(trajs(:,1,k),-trajs(:,2,k),'.','MarkerEdgeColor',colors{l})
        text(trajs(1,1,k),-trajs(1,2,k)+25,string(l))
    end
    
    if pinMode == 7
        [~,CTpoint] = min(min(trajs(:,2,:),[],1),[],3);
        [~,CBpoint] = max(max(trajs(:,2,:),[],1),[],3);

        ix = setdiff(1:7,[CTpoint,CBpoint]);
    elseif pinMode == 5
       ix = 1:5; 
    elseif pinMode == 4 || pinMode == 2
       ix = 1:numel(track);
    end
    
    if pinMode > 2
        [~,Lpoint] = min(min(trajs(:,1,ix),[],1),[],3);
        [~,Rpoint] = max(max(trajs(:,1,ix),[],1),[],3);
    end

    [~,Tpoint] = min(min(trajs(:,2,ix),[],1),[],3);
    [~,Bpoint] = max(max(trajs(:,2,ix),[],1),[],3);
    
%     [Bpoint] = setdiff(1:4,[Lpoint,Rpoint,Tpoint]);
%     Bpoint = Bpoint(1);
    if pinMode == 7
        clampStretch{l} = vecnorm(trajs(:,:,CTpoint) - trajs(:,:,CBpoint),2,2);
    end
    vertStretch{l} = vecnorm(trajs(:,:,ix(Tpoint)) - trajs(:,:,ix(Bpoint)),2,2);
    
    if pinMode > 2
        horStretch{l} = vecnorm(trajs(:,:,ix(Lpoint)) - trajs(:,:,ix(Rpoint)),2,2);
    end
%     poisson = diff(horStretch{l}) ./ diff(vertStretch{l});
%     vertStretch{l} = (vertStretch{l}) ./ vertStretch{l}(1);
%     horStretch{l} = (horStretch{l}) ./ horStretch{l}(1);
    
    if pinMode == 7
        traj{l} = trajs(:,:,[ix(Rpoint),ix(Lpoint),ix(Bpoint),ix(Tpoint),ix(setdiff(1:5,[Lpoint,Rpoint,Tpoint,Bpoint])),CBpoint,CTpoint]);
    elseif pinMode == 5
        traj{l} = trajs(:,:,[ix(Rpoint),ix(Lpoint),ix(Bpoint),ix(Tpoint),ix(setdiff(1:5,[Lpoint,Rpoint,Tpoint,Bpoint]))]);
    elseif pinMode == 4
        traj{l} = trajs(:,:,[ix(Rpoint),ix(Lpoint),ix(Bpoint),ix(Tpoint)]);
    elseif pinMode == 2
        traj{l} = trajs(:,:,[ix(Bpoint),ix(Tpoint)]);
    end
        
    figure(2)
    % clf
    hold on
    ms = 17;
    if pinMode == 7
        clamp = scatter(forces{l},clampStretch{l},ms,'s','filled','MarkerEdgeColor','k','MarkerFaceColor',colors{l});
    end
    vert = scatter(forces{l},vertStretch{l},ms,'o','filled','MarkerEdgeColor','k','MarkerFaceColor',colors{l});
    if pinMode > 2
        hor = scatter(forces{l},horStretch{l},ms,'d','filled','MarkerEdgeColor','k','MarkerFaceColor',colors{l});
    end
    %     ylabel("Relative Strain")
    ylabel("Extension [pixels]")
    xlabel("Force Gauge [Newtons]")
    
%     figure(3)
%     hold on
%     plot(forces{l}(2:end),poisson,'Color',colors{l});
%     ylabel("Differential Poisson Ratio")
%     xlabel("Stress [Newtons]")
end
figure(2)

if pinMode == 7
    legend([clamp,vert,hor],{'Clamp','Vertical','Horizontal'},'location','east')
elseif pinMode == 5 || pinMode == 4
    legend([vert,hor],{'Vertical','Horizontal'},'location','east')
elseif pinMode == 2
    legend([vert],{'Vertical'},'location','east')
end

ax = gca; ax.FontSize = 20;
xlim([0,30])

figure(3)
clf
hold on
for l = 1:numel(forces)
    scatter((1:numel(forces{l}))/30,forces{l},'.','MarkerEdgeColor',colors{l})
end
ylim([0,30])
xlabel("Time [seconds, roughly]")
ylabel("Force Gauge Reading")

% title("Sample labelled 'Sam'")

hName = fnBase + 'outputData.h5';
if isfile(hName)
    delete(hName)
end

trajInterp = traj;

for l = 1:numel(forces)
    if doInterp
        for k = 1:pinMode
            firstFrame = find(all(~isnan(traj{l}(:,:,k)),2),1);
            lastFrame = find(all(~isnan(traj{l}(:,:,k)),2),1,'last');
            xq = firstFrame:lastFrame;
            nn = find(all(~isnan(traj{l}(:,:,k)),2));
            xqAll = 1:size(traj{l},1);
            trajInterp{l}(xq,1,k) = interp1(xqAll(nn),traj{l}(nn,1,k),xq);
            trajInterp{l}(xq,2,k) = interp1(xqAll(nn),traj{l}(nn,2,k),xq);
        end
    end
    traj = trajInterp;
    fprintf("\nClean pin tracking fraction: %f\n",numel(find(all(all(~isnan(traj{l}),3),2))) / size(traj{l},1));
    fprintf("Clean force gauge reading fraction: %f\n",numel(find(~isnan(forces{l}'))) / numel(forces{l}));
    threePointDiff = abs(forces{l} - [forces{l}(2:end), NaN]) + abs(forces{l} - [NaN,forces{l}(1:end-1)]);% + abs(forces{l} - [NaN,NaN,forces{l}(1:end-2)]) + abs(forces{l} - [forces{l}(3:end),NaN,NaN]);
    fprintf("Quasistatic (3 in a row constant force) fraction: %f\n",numel(find(threePointDiff' == 0)) / numel(forces{l}));
    
    goodFrame = ~isnan(forces{l}') & all(all(~isnan(traj{l}),3),2);
    
    if requireQuasistatic
        goodFrame = goodFrame & (threePointDiff' == 0);
    end
    fullFrames = find(goodFrame);
    fprintf("First frame: %d\n",fullFrames(1));
    fprintf("Accepted frames fraction: %f\n",numel(fullFrames) / numel(forces{l}));
    
    h5create(hName,sprintf('/video%d/frames',l),[Inf],'ChunkSize',[100]);
    h5write(hName,sprintf('/video%d/frames',l),fullFrames,[1],[numel(fullFrames)]);
    
    h5create(hName,sprintf('/video%d/traj',l),[Inf 2 pinMode],'ChunkSize',[100 2 pinMode]);
    h5write(hName,sprintf('/video%d/traj',l),traj{l}(fullFrames,:,:),[1 1 1],size(traj{l}(fullFrames,:,:)));
    
    h5create(hName,sprintf('/video%d/forces',l),[Inf],'ChunkSize',[100]);    
    h5write(hName,sprintf('/video%d/forces',l),forces{l}(fullFrames),[1],[numel(fullFrames)]);
    
    figure(4)
    hold on
    if pinMode == 7
        scatter(forces{l}(fullFrames),clampStretch{l}(fullFrames),'.','MarkerEdgeColor',colors{l})
        ylabel("Clamp Extension [pixels]")
    elseif pinMode == 5 || pinMode == 4 || pinMode == 2
        scatter(forces{l}(fullFrames),vertStretch{l}(fullFrames),'.','MarkerEdgeColor',colors{l+1})
        ylabel("Vertical Pin Extension [pixels]")
    end
    xlabel("Force Gauge Reading [N]")
    
end
h5writeatt(hName,'/','description', description);
h5create(hName,'/clampSize', 1); h5write(hName,'/clampSize',clampSize);
h5create(hName,'/fabricSize', 1); h5write(hName,'/fabricSize',fabricSize);

if isfile(sprintf('experimentData/outputData_%s.h5',description))
    fprintf("Output file already exists! Skipping.\n")
    return
else
    copyfile(hName,sprintf('experimentData/outputData_%s.h5',description))
end

% save(fnBase + 'outputData','forces','vertStretch','horStretch','traj','description')
