close all
% fnBase = "C:\\Users\\Sam\\Documents\\GitHub\\uniaxial\\ignore\\Rib_x_acr_new\\"; % Folder containing video1.mp4 through video5.mp4
fnBase = "C:\\Users\\Sam\\Downloads\\";
frameReduceFactor = 4;
% Import .MOV with command:
%
%   for($i=1; $i -le 5; $i++){ffmpeg -i video$i.MOV -vcodec h264 -pix_fmt yuv420p -acodec aac video$i.mp4}
%


% processFile(fnBase,"graphiteTrack",frameReduceFactor)

for v = 5:5
    fnHead = "stockxout";
    processFile(fnBase,fnHead,frameReduceFactor)
end

function processFile(fnBase,fnHead,frameReduceFactor)
    inFileExt = ".mp4";
    fn = sprintf("%s%s%s",fnBase,fnHead,inFileExt);
    vidRead = VideoReader(fn);
    vidWrite = VideoWriter(sprintf("%s%scrop.avi",fnBase,fnHead),'Grayscale AVI');

    frameOne = read(vidRead,1);
    frame = read(vidRead,Inf);
    imshow(frame/2 + frameOne/2)
    
    % Include: Fabric width (corners, at the final frame of stretching).
    % Clamp "height" (short way, top plane of thing only).
    fprintf("1. Draw box around entire area of sample\n2. Draw box around entire area of force gauge screen\n")
    croptangle = drawrectangle;
    cropPos = croptangle.Position;
    
    
    forcetangle = drawrectangle;
    forcePos = forcetangle.Position;
    
    close
    vidRead = VideoReader(fn);
    open(vidWrite)

    framesKeep = 1:1:vidRead.NumFrames;
    forces = zeros(1,numel(framesKeep));

    tick = 0;
    for k = 1:numel(framesKeep)
        if mod(tick,50) == 0
            fprintf("%d/%d\n",framesKeep(k),framesKeep(end))
        end
        tick = tick + 1;
        I = read(vidRead,framesKeep(k));
    %     imshow(I)
        forces(k) = getForce(imcrop(I,forcePos));
        if mod(k-1,frameReduceFactor) == 0
            I = imcrop(I,cropPos);

%             I = imcomplement(I(:,:,2)/2 + I(:,:,3)/2);
            I = I(:,:,1) - I(:,:,2);
%             I = imcomplement(im2gray(I));
            writeVideo(vidWrite,I);
        end
    end
    close(vidWrite)
    save(sprintf("%sforce_%s",fnBase,fnHead),'forces','frameReduceFactor')
    fprintf("Force recovery rate: %f\n",numel(find(~isnan(forces))) / numel(forces))
end

function force = getForce(I)
    manualFiltering = false;
    % I = imgaussfilt(I,10);

%     BW = edge(rgb2gray(I),'canny');

%     clf

%     Ired = I(:,:,1);
    
%     imshow(Ired)
%     drpt = drawpoint;
%     pt = drpt.Position;

%     if Ired(1457,625) < 3.8e4
%         pt = [625,1457];
%     else
%         pt = [625,1557];
%     end
%     mask = bwselect(Ired < 4e4,pt(1),pt(2));
    
%     CC = bwconncomp(Ired < 5);
    CC = bwconncomp(I(:,:,2) - I(:,:,1) > 50);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    mask = zeros(size(I(:,:,1)));
    mask(CC.PixelIdxList{idx}) = 1;
    maskOrn = regionprops(mask,'Orientation');
    rotI = imrotate(I,-maskOrn.Orientation);
    rotMask = imrotate(mask,-maskOrn.Orientation);
    
    maskBB = regionprops(rotMask,'BoundingBox');

%     ssDisp = drawrectangle;
%     ssDispI = rot90(imcrop(I,ssDisp.Position),-2);

    ssDispI = rot90(imcrop(rotI,maskBB.BoundingBox),-2);

%     ssDispBW = rgb2gray(ssDispI) < 70;
%     ssDispBW = rgb2gray(ssDispI) < 80;
    ssDispBW = imcomplement(imbinarize(ssDispI(:,:,2)));
    
    %%
%     figure(2)
%     figure(1)
%     imshow(ssDispI)
    ssDispBW = imerode(ssDispBW,strel('disk',3,8));
%     figure(1)
%     imshow(ssDispBW)
    [L,n] = bwlabel(ssDispBW);
    regProps = regionprops(L,'BoundingBox','Centroid','Area');

    imHeight = size(ssDispBW,1);
    imWidth = size(ssDispBW,2);

    relHeight = zeros(1,n);
    for k = 1:n
%         bBox = regProps(k).BoundingBox;
        cent = regProps(k).Centroid;
        relHeight(k) = cent(2)/imHeight;
        ar = regProps(k).Area;
        if ar < 3
            relHeight(k) = NaN;
        end
    end
%     [~,triROI] = max(relHeight);
%     [~,topROI] = min(relHeight);

%     images.roi.Rectangle(gca,'Position',regProps(triROI).BoundingBox);
%     images.roi.Rectangle(gca,'Position',regProps(topROI).BoundingBox);

    % segWidth = regProps(topROI).BoundingBox(3);
%     triX = regProps(triROI).Centroid(1);
%     triY = regProps(triROI).Centroid(2);
%     triWidth = regProps(triROI).BoundingBox(3);
%     triHeight = regProps(triROI).BoundingBox(4);

    slant = 0.022 * imWidth;
    
%     digitSpacing = 43/12 * triWidth;
%     digitSpacing = 47/12 * triWidth;
    
    digitSpacing = 0.175 * imWidth;
    
%     topBottomOffset = 63/16 * triHeight;
%     topBottomOffset = 69/16 * triHeight;
    topBottomOffset = 0.4830 * imHeight;
%     leftRightOffset = 22/12 * triWidth;
    leftRightOffset = 0.0896 * imWidth;

%     triToTopSegCenters = 83.1/16 * triHeight;
%     triToTopSegCenters = 92.1/16 * triHeight;
%     triToTopSegCenters = 0.6447 * imHeight;
    triToTopSegCenters = 0.6047 * imHeight;
%     segThickness = 10/23 * triHeight;
    segThickness = 0.0487 * imHeight;

    decoder = [1,1,1,1,1,1,0;
        0,1,1,0,0,0,0;
        1,1,0,1,1,0,1;
        1,1,1,1,0,0,1;
        0,1,1,0,0,1,1;
        1,0,1,1,0,1,1;
        1,0,1,1,1,1,1;
        1,1,1,0,0,0,0;
        1,1,1,1,1,1,1;
        1,1,1,1,0,1,1];


    segs = zeros(4,7);
    outDig = nan(4,1);
    for dig = 1:4
        digRef = [(dig - 2) * digitSpacing + (0.44 * imWidth),-triToTopSegCenters + (0.88 * imHeight)];
        for seg = 1:7
            segOff = [0,0];
            switch seg
                case 1
                    segOff = [slant,0];
                case 2
                    segOff = [leftRightOffset/2 + 3*slant/4,topBottomOffset/4];
                case 3
                    segOff = [leftRightOffset/2 + 1*slant/4,3*topBottomOffset/4];
                case 4
                    segOff = [0,topBottomOffset];
                case 5
                    segOff = [-leftRightOffset/2 + 1*slant/4,3*topBottomOffset/4];
                case 6
                    segOff = [-leftRightOffset/2 + 3*slant/4,topBottomOffset/4];
                case 7
                    segOff = [slant/2,topBottomOffset/2];
            end
            target = segOff + digRef;
            tolerance = vecnorm([10,10])/1.8;
%             hold on
%             scatter(target(1),target(2),"rx")
%             images.roi.Circle(gca,'Center',target,'Radius',tolerance,'Color','r');
            for l = 1:n
                actual = regProps(l).Centroid;
                if seg == 1
%                     scatter(actual(1),actual(2),"bo")
%                     images.roi.Circle(gca,'Center',actual,'Radius',tolerance,'Color','b');
                end
                if vecnorm(target - actual) < tolerance
                    segs(dig,seg) = 1;
                end
            end
        end
        for d = 1:10
            if segs(dig,:) == decoder(d,:)
                outDig(dig) = d-1;
                break
            end
        end
        if all(segs(dig,:) == 0)
            outDig(dig) = 0;
        end
        if isnan(outDig(dig)) && manualFiltering
            th = segThickness * 4;
            figure(2)
            imshow(imcrop(ssDispI,[digRef(1) - leftRightOffset/2 - th, digRef(2) - th, leftRightOffset + 2 * th, topBottomOffset + th *2]));
            digSelect = inputdlg("What is this number?");
            outDig(dig) = str2double(digSelect);
            close
        end
    end
%     outDig
%     drawnow
    outVal = outDig' * (10 .^ [1,0,-1,-2]');
    force = outVal;
    if outDig(2) == 8
%         imshow(ssDispI)
%         for dig = 1:4
%             digRef = [(dig - 2) * digitSpacing + (0.4469 * imWidth),-triToTopSegCenters + (0.88 * imHeight)];
%             for seg = 1:7
%                 segOff = [0,0];
%                 switch seg
%                     case 1
%                         segOff = [slant,0];
%                     case 2
%                         segOff = [leftRightOffset/2 + 3*slant/4,topBottomOffset/4];
%                     case 3
%                         segOff = [leftRightOffset/2 + 1*slant/4,3*topBottomOffset/4];
%                     case 4
%                         segOff = [0,topBottomOffset];
%                     case 5
%                         segOff = [-leftRightOffset/2 + 1*slant/4,3*topBottomOffset/4];
%                     case 6
%                         segOff = [-leftRightOffset/2 + 3*slant/4,topBottomOffset/4];
%                     case 7
%                         segOff = [slant/2,topBottomOffset/2];
%                 end
%                 target = segOff + digRef;
%                 hold on
%                 scatter(target(1),target(2),"rx")
%             end
%         end
%         outDig
%         ;
    end
end

