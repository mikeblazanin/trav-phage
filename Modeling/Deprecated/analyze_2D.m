% load('onedimensional3.mat');
    
video = VideoWriter('wavemotion.avi');
        video.FrameRate = 10;
        open(video);
        for t = 15:15:1620
            frame_file = ['wavemotion_',num2str(t),'.jpg'];
            if ~exist(frame_file, 'file')
                break
            end
            img = imread(frame_file);
            writeVideo(video, img);
        end
        close(video)

% figure
% hold on
% yyaxis right
% plot(LX(locs_min:2*locs_max-locs_min),f(locs_min:2*locs_max-locs_min))
% xlabel('Distance (mm)')
% ylabel('Perceived Signal')
% figure
% hold on
% plot(LX(locs_min:2*locs_max-locs_min),gradient(locs_min:2*locs_max-locs_min),'--')
% 
% gradient = [0 ((f(3:end)-f(1:end-2))/2/dx) 0];
% 
% for i = 1:length(CWBias)
%     grad_at_peak(i) = gradient(max_ind(i)+locs_min);
% end
% c = MFTmodel.chi .* grad_at_peak;
% figure
% hold on
% plot(CWBias,grad_at_peak)
% figure
% hold on
% plot(CWBias,grad_at_peak)
% plot(CWBias,MFTmodel.chi)
% figure
% hold on
% plot(CWBias,c)