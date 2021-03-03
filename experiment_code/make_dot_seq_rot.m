% make_dot_seq_rot.m
%
% Generate sequence of dot coordinates, each {pair of rows? 3rd dim?}
% submitted as a 'frame' to DrawDots in a dot-stimulus animation
%
% Coordinates span [-1, 1], and speed is normalized to unit range (so need
% to convert to [unit]/fr. Coords can be multiplied by aperture_radius_pix
% outside of this function. Can provide a random seed, though rng should be
% initialized outside of this function when used in psychophysics, etc.
%
% (this could/should be used for timing-intensive applications, where dot
% sequences need to be generated rapidly and presented w/ high-precision
% timing. this will be better than updating dot positions on each frame,
% especially for VPixx Pro presentation)
%
% DIRECTION is wrt clockwise tangential motion, in deg. 0 is clockwise, 180
% is ccw, 90 is inward, 270 is outward. If DIRECTION is n_frames long, then
% different direction can be used to update on each frame. <FOR NOW> motion
% is wrt coords [0,0]; if different 'origin' is necessary, use
% aperture_shape below. 
%
% DIRECTION, COHERENCE, and SPEED can all be n_frames long; if so, then
% sequence is generated with dynamic properties. Otherwise, it's generated
% with constant values for each. 
%
% APERTURE_SHAPE is either 0 (circle), 1 (square), or dimensions of a
% rectangle, where 1 is treated as above (for coordinate span). Maximum of 
% APERTURE_SHAPE should be 1. Oval not quite supported yet; could try scaling 
% the coords outside of this function?
%
% COHERENCE: 1-coherence is % of dots with uniform random distribution of
% directions, updated at end of their life. For scenarios where coherence
% should be defined as % of dots moving in one of 2 directions (cw/ccw, for
% example), call this function twice w/ 2 different #'s of dots and draw
% both. 
%
% NOTE: NDOTS will appear inside the aperture; not including inner
% apertures here. So, for fixation-centered stimuli employing an aperture
% around fixation, NDOTS won't be precise. Could add this later on?
%
% NOTE: only circular aperture shape supported now.... (not easy to
% intuit how to replot dots for square/rectangular apperture; that only
% makes sense, for now, for planar motion)
%
% Overall very similar to Graziano et al, 1994 motion stimuli, including
% scaling for distance from center by 0.3*radius, and randomly choosing new
% position when going inside/outside
% http://www.jneurosci.org/content/jneuro/14/1/54.full.pdf
%
% TCS 3/14/2017 (pi)
% - made it
%
% USAGE: [dot_seq] = make_dot_seq_rot(ndots, direction, nframes [,speed] [,coherence] [,dot_life] [aperture_shape])
%
% DEFAULTS:
% - SPEED: 0.1/fr
% - COHERENCE: 1 (100%)
% - DOT_LIFE: 5 fr
% - APERTURE_SHAPE: 0 (circle)

% TODO: fix this error!!!!!
% 
% Index in position 2 exceeds array bounds.
% 
% Error in make_dot_seq_rot (line 224)
%     new_pos(:,to_reloc) = init_xy(:,1:sum(to_reloc));
% 
% Error in fspri_pilot (line 294)
%         p.ds{tt} = make_dot_seq_rot(p.ndots,p.motion_cond(p.conditions(tt,5)),p.fps*p.stim_dur,p.dot_speed,p.motion_coh,p.dot_life_frames);


function [dot_seq] = make_dot_seq_rot(ndots,direction,nframes,speed,coherence,dot_life)

%% CHECK ARGUMENTS <<<<NOTE: update to check whether the argument exists>>>
if nargin < 3
    error('dotTools:make_dot_seq_rot:insufficentArguments','Insufficient arguments defined; need to define # of dots, direction, and # frames');
    return
end

% default speed
if nargin < 4
    speed = 0.1*ones(nframes,1);
end

% default coherence (100%)
if nargin < 5
    coherence = 1*ones(nframes,1);
end

% default dot_life
if nargin < 6
    dot_life = 5;
end


%% CHECK ARGUMENT SIZE
%  if direction, speed, coherence are scalars, multiply by ones(nframes,1);

if length(direction)==1
    direction = direction*ones(nframes,1);
end

if length(speed)==1
    speed = speed*ones(nframes,1);
end

if length(coherence)==1
    coherence = coherence*ones(nframes,1);
end

%% INITIALIZE VARIABLES

% full sequence of dot positions, 2 x n_dots x n_frames
dot_seq = nan(2,ndots,nframes);

% initialize dot age; if mod(age,ndots)==0, ensure fully uniform;
% otherwise, pad w/ mod(age,ndots) rounded random values
% (NOTE: this would be different # of rand operations; should be ok...
% could also just randomize like: dot_age = ceil(rand(ndots,1)*dot_life);

dot_age = [repmat(1:dot_life,1,floor(ndots/dot_life)) ceil(rand(1,  mod(ndots,dot_life)) * dot_life)  ];

% initialize dot direction; will round coherence*ndots up
%dot_dir = [direction*ones(1,ceil(coherence*ndots)];
% FOR NOW: if both coherence & direction are constant, repmat; otherwise,
% loop over frames
%
% NOTE: directions random from 0:360, but exactly uniformly sampled; could
% add a random_dir flag that makes these actually random/continuous, but
% not guaranteed to be uniform
%
% NOTE: when coherence changes, the direction of non-coherent dots, right
% now, will instantaneously change subtly, giving them curved trajectory.
% can add a mode that checks whether there are nframes distinct values of
% coherence, in which case can do something more complicated, one value
% (easy), or fewer than ndots discrete coherences (shuffle directions each
% time coherence changes? base this on dot life?)

% dot_dir: n_frames x n_dots
if length(unique(coherence)) == 1 && length(unique(direction))==1
    
    dot_dir = [repmat(direction,1,ceil(coherence(1)*ndots)) repmat( linspace(360/floor( (1-coherence(1))*ndots ), 360, floor((1-coherence(1))*ndots)), nframes, 1) ];   %repmat([ direction * ones(1,ceil( ],nframes,1);

    
else % SLOWER, but still pretty fast
    
    % initialize
    dot_dir = nan(nframes,ndots);
    
    for ii = 1:nframes
        dot_dir(ii,:) = [ direction(ii) * ones(ceil(coherence(ii))*ndots,1) linspace(360/floor( (1-coherence(ii))*ndots ), 360, floor((1-coherence(ii))*ndots) )   ];
    end

end


% first dot positions:
% - randomly generate double the number of dots needed, that way we can
%   cull them based on aperture shape
init_xy = 2*rand(2,ndots*5)-1;  % spans [-1 1]

% cull dots based on aperture shape:
% if length(aperture_shape)>1 % RECTANGLE
%     % make sure abs(x) & abs(y) both fit within width and height
%     init_xy = abs(init_xy(1,:))<=aperture_shape(1) & abs(init_xy(2,:))<=aperture_shape(2);
    
% elseif aperture_shape == 0 % CIRCLE
    % make sure radius <= 1
    init_xy = init_xy(:,sqrt(sum(init_xy.^2,1))<=1);
    
% elseif aperture_shape == 1 % SQUARE
%     % don't need to cull
% end


if size(init_xy,2) < ndots % if we screwed up and generated too few dots...
    error(); % SHOULD NEVER GET HERE!!!!
end

dot_seq(:,:,1) = init_xy(:,1:ndots);

clear init_xy;

%% COMPUTE DOT SEQUENCE

for ff = 2:nframes
    
    % replot dead dots
    dead_dots = dot_age == dot_life; % or >????
    
    % new dot positions
    init_xy = 2*rand(2,sum(dead_dots)*5)-1;
    % cull new dots based on aperture shape (adapted from above...):
%     if length(aperture_shape)>1 % RECTANGLE
%         init_xy = abs(init_xy(1,:))<=aperture_shape(1) & abs(init_xy(2,:))<=aperture_shape(2);
%     elseif aperture_shape == 0 % CIRCLE
        init_xy = init_xy(:,sqrt(sum(init_xy.^2,1))<=1);
%         if size(init_xy,2) < sum(dead_dots) % if we screwed up and generated too few dots...
%             error(); % SHOULD NEVER GET HERE!!!!
%         end
%     elseif aperture_shape == 1 % SQUARE
%         % nothing to do here
%     end
    
    if sum(dead_dots)>0
        dot_seq(:,dead_dots,ff) = init_xy(:,1:sum(dead_dots));
    end
    dot_age(dead_dots) = 0; % or 1????
    
    clear init_xy;
    
    % update position of alive dots
    
    alive_dots = ~dead_dots; % use these indices to insert new dot pos, update life
    
    % get polar angle of alive dots (sum(alive_dots) long)
    dot_ang = atan2d(dot_seq(2,alive_dots,ff-1),dot_seq(1,alive_dots,ff-1));
    
    dot_r = sqrt(sum(dot_seq(:,alive_dots,ff-1).^2,1));
    new_r = dot_r + speed(ff).*sind(-1*dot_dir(ff,alive_dots));
    
    % BELOW: no scaling w/ r
    %dxy = speed(ff).*  [ cosd(dot_ang + dot_dir(ff,alive_dots) + 90); sind(dot_ang + dot_dir(ff,alive_dots) + 90) ];
    
    dxy = speed(ff) .*  0.3 .* dot_r .* [ cosd(dot_ang + dot_dir(ff,alive_dots) + 90); sind(dot_ang + dot_dir(ff,alive_dots) + 90) ];
    
    new_pos = dot_seq(:,alive_dots,ff-1) + dxy;
    
    % find new_r > 1, set to 0; find new_r < 0, set to 1
    to_reloc = new_r>1 | new_r<=0;
        
    init_xy = 2*rand(2,sum(to_reloc)*10)-1;
    init_xy = init_xy(:,sqrt(sum(init_xy.^2,1))<=1);
    new_pos(:,to_reloc) = init_xy(:,1:sum(to_reloc));
    %new_pos(:,new_r>1)  = 0.001 * [cosd(dot_ang(new_r>1));sind(dot_ang(new_r>1))]; % dont' redraw at 0,0 because there's undefined atan2 there
    %new_pos(:,new_r<=0) = 1 * [cosd(dot_ang(new_r<=0));sind(dot_ang(new_r<=0))];
    

    dot_seq(:,alive_dots,ff) = new_pos;
    
    dot_age = dot_age+1;
    clear dead_dots alive_dots dxy dot_ang new_pos dot_r new_r;
    
end

return