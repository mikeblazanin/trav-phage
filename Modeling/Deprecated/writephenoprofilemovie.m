function writephenoprofilemovie(SimParams, Names)
% This function assembles frame images created in drawphenoprofilemovieframe.m to create a movie.

disp('Making movie ...')
if ~exist(Names.phenoProfileMovieDir, 'dir')
    mkdir(Names.phenoProfileMovieDir)
end
video = VideoWriter([Names.phenoProfileMovieDir, 'pheno_profile.avi']);
video.FrameRate = 10;
open(video);
for t = 1:SimParams.nT
    disp(['Loading movie frame ', num2str(t), ' ...'])
    frameFile = [Names.phenoProfileMovieDir, 'pheno_profile_',num2str(t*SimParams.outDt),'.jpg'];
    if ~exist(frameFile, 'file')
        disp(['Moving ends at frame ' num2str(t)])
        break
    end
    img = imread(frameFile);
    writeVideo(video, img);
end
disp('Saving movie ...')
close(video)
end