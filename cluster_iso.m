function cluster_iso(filename,connectivity)
%% This code is to separate clusters from a complicated microstructure. This code will generate N (cluster#) .mat files with one cluster each
% NU_GT project
% Author: Xiaolin Li
% Advisor: Prof. Brinson & Prof. Chen
% Date: Jun 22, 2015

%Var:
%filename: file name (contains 3D matrix, filler is labeled as 1, matrix is
%           labeled as 0)
%connectivity: The number of neighboring points to be considered as connected.
%               used in the bwlabel function. could be one of [6,18,26]
%% Code starts here
clearvars -except filename connectivity
if (connectivity ~= 6) && (connectivity ~=18) && (connectivity~=26)
    disp('Error: The connectivity you entered is invalid!');
    exit(0);
end
folder = strcat(filename(1:end-4),'_clusters');
mkdir(folder)
file = strcat(filename);
load(file);
Bimg = Bimg_coarse;
[L, num] = bwlabeln(Bimg,connectivity);
for i= 1: num
    cluster = double(L==i);
    savefile = strcat('cluster_',int2str(i));
    save(strcat('./',folder,'/',savefile),'cluster');
end
exit();
end
