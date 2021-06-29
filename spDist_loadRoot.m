% spDist_loadRoot.m
%
% returns the root directory - easier for distributing code, etc

function root = spDist_loadRoot

rootdir = 'spDist';

if ismac
    % updated 9/24/2020 - for Grace
    root = sprintf('/share/datb/%s/',rootdir); %updated 23 march 2021, after spDist was moved to datb for storage reasons
else
    this_computer = char(java.net.InetAddress.getLocalHost.getHostName);
    if strcmpi(this_computer,'tcs-compute-1')
        root = '/home/local/PSYCH-ADS/sprague/labshare/projects/nyu/spDist/';
    elseif strcmpi(this_computer,'tcs-precision')
        root = 'Z:/projects/nyu/spDist/';
    else % for vader.psych.nyu.edu
        root = '/deathstar/datb/spDist/'; %updated 05/04/21 to point to new storage location
    end
end

return