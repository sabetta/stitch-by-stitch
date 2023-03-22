close all
fnBase = "C:\\Users\\Sam\\Documents\\GitHub\\uniaxial\\ignore\\Rib_x_acr_new\\"; % Folder containing video1.mp4 through video5.mp4
inFileExt = ".mp4";

for v = 1:1
    fnHead = "video1";
%     fnHead = "Seed_y_fixed";
    fn = sprintf("%s%s.mp4",fnBase,fnHead);
    vidRead = VideoReader(fn);

    frame = read(vidRead,Inf);
    f = figure;
    f.WindowState = 'maximized';
    imshow(frame)
    ax = gca;
    ax.YLim = [0 (size(frame,1)*3/5)];
    ax.XLim = [(.0 * size(frame,2)) (.9 * size(frame,2))];
    % Include: Fabric width (corners, at the final frame of stretching).
    % Clamp "height" (short way, top plane of thing only).
    fprintf("1. Draw line representing top of fabric\n2. Draw line across clamp\n")
    fabricLine = drawline;
    fabricPos = fabricLine.Position;
    fabricSize = vecnorm(diff(fabricPos,1));
    
    clampLine = drawline;
    clampPos = clampLine.Position;
    clampSize = vecnorm(diff(clampPos,1));
    
    close

    save(sprintf("%smeasurements_%s",fnBase,fnHead),'fabricSize','clampSize')
end