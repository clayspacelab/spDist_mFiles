% spDist_loadRoot.m
%
% returns the root directory - easier for distributing code, etc

function root = spDist_loadRoot

rootdir = 'spDist';

if ismac
    % updated 9/24/2020 - for Grace
    root = sprintf('/share/data/%s/',rootdir);
else
    this_computer = char(java.net.InetAddress.getLocalHost.getHostName);
    if strcmpi(this_computer,'tcs-compute-1')
        root = '/mnt/LabShare/projects/nyu/spDist/';
    elseif strcmpi(this_computer,'tcs-precision')
        root = 'Z:/projects/nyu/spDist/';
    else % for vader.psych.nyu.edu
        root = '/deathstar/data/spDist/';
    end
end

return