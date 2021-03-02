% spDist_scanner.m
%
% ??? 0.75 s TRs
%
% 8/29/2018 - TCS & GH, adapted from MR code
%
% Stimulus presentation script for spatial distractor experiment. 
% 2 conditions, pre-cued (magenta & cyan) - distractor, no-distractor
% Each trial begins with a pre-cue indicating whehter a distractor will be
% presented during the delay period (4.5 s). If a distractor is presented,
% participants must attend to the distractor location and make a
% challenging discrimination (variable-coherence spiraling random dots).
% Then, after a second delay period, participants report the remembered
% location with a memory-guided saccade. Distractor position is
% pseudo-randomized (one of 7 relative positions compared to target
% position, with jitter). 
%
% TIMING INFO:
%
%
%
% XDAT INFO:
%
%
%
%
%
%
% updated TCS 9/13/2017 - reacts to variable TR (for mb experiments, p.TR
% should usually be some # of TRs - here, 4x0.75=3; could also use shorter,
% but not all that necessary pending how task timing is set up)
%
% TCS 6/9/2017
% TCS 9/26/2017 - 500 Hz sample rate
%
% TODO: stop recording w/ ESC
% TODO: fixation as dot within circle


function spDist_scanner(subj,run,task_coh)

try
p.expt_name = 'spDist_scanner';

p.do_et = 1;

p.TR = 3; % 4x multiband, so measured TR is 0.75, but "TR" for stim is 3

p.subj = subj;
p.run = run;
p.task_coh = task_coh;

 p.scanner = 1;
%p.scanner = 0
if ~exist('./data','dir')
    mkdir('./data');
end

p.filename = sprintf('./data/%s_r%02.f_%s_%s.mat',p.subj,p.run,p.expt_name,datestr(now,30));
if p.do_et == 0
    p.eyedatafile = sprintf('%s_D%02.f',p.subj(1:min(length(p.subj),3)),p.run);
end

p.rng_seed = cputime*1000;
rng(p.rng_seed);


% define geometry of WM stimuli
% ------ size of relevant stim features, etc ------ %
p.wm_ecc = 12;     % deg [in behavioral room, max ecc of circle 15 deg]
p.cue_size = 0.55; % deg
p.wm_size = 0.65;  % deg, size of WM dots

% how far to be away from meridians (target pos); polar angle deg
p.merid_sep = 5; % must be this far (polar angle deg) away from vert/horiz

p.aperture_size = 15; % [or max ecc?]

%p.n_pos = 16; % and we'll use 2 offset angles

p.fix_size_in  = 0.075; % radius, deg
p.fix_size_out = 0.30; % radius, deg
p.fix_pen = 1.5;

p.fix_size_mult = 1.25; % for pre/post experiment, fix is bigger



% ------ color of relevant stim features ----------- %
% TODO: adjust??
p.bg_color  = 60;%20*[1 1 1];
p.fix_color = 90;%%75*[1 1 1];%[150 150 150];        % during trial/delay/etc

p.wm_color = p.fix_color;

p.go_color = p.fix_color;% 130*[1 1 1];%[255 255 255]; % when subj should choose, color of fix

p.dim_amt = 0.9; % multiply fix_color by this during ITI, start, end periods

% magenta:cyan RGB ratio: 2.46/1.18 (cyan must be 1.18/2.46x as intense as
% magenta)
p.cue_colors = [1 0 1;0 1 1]; % TODO: equiluminant! (magenta = no dist, cyan = dist)
p.cue_rel_lum = 0.5; % relative intensity, both colors scaled by this; magenta is this*255, cyan is this*1.18/2.46*255
p.cue_colors(1,:) = round(p.cue_colors(1,:)*p.cue_rel_lum*255);
p.cue_colors(2,:) = round(255*p.cue_rel_lum*(1.18/2.46)*p.cue_colors(2,:));


p.dist_feedback_colors = round([0 0.44*255 0;255 0 0;.305*255 .305*255 0]); % green, red, yellow for corr, incorr, miss




% define geometry, dot properties, etc of distractors
% -------- distractor properties ------------ %

p.dist_jitter = 12; % ABSOLUTE randomize distractor position by + or - this (polar ang deg)

p.dist_rad = 1;         % 1 dva radius
p.dist_dotsize = 0.075; % dva?
p.dist_density = 17;    % dots/dva^2
p.dist_speed = 3;       % not exact units, due to rotation... (can convert to deg polar ang/s)

p.dist_dotlife = .1; % sec (12 120 Hz frames)

p.dist_dir_proto = [0 180]; % CW (resp=1), CCW (resp=2) (NOTE: inconsistent w/ docs...)

p.dist_dotcolor = [2 2 2;0 0 0]*p.bg_color; % so these are at 100% contrast relative to bg


% ------ conditions ------ %

p.p_dist = 0.7; % probability that distractor occurs (used to compute # of trials)
p.n_d_offsets = 7; % # of relative distractor locations
p.dist_step = 360/p.n_d_offsets;
p.repetitions = 1; % only 1 repetition of each distractor location

d_offsets = (-1*(p.n_d_offsets-1)/2):1:((p.n_d_offsets-1)/2); % distractor offset steps (centered at 0; only overlaps w/ 0 if odd)

% generate blank p.conditions that's this big, then consider blanks
% no-distractor
p.ntrials = p.repetitions*p.n_d_offsets/p.p_dist; % # of trials w/ no distractor

% column 1: distractor/no-distractor (2/1)
% column 2: relative distractor position (-3:3) or NaN
p.conditions = [ones(p.ntrials,1) nan(p.ntrials,1)];
% fill in distractors
idx = 1;
for ii = 1:p.repetitions
    for dd = 1:length(d_offsets)
        p.conditions(idx,1) = 2;
        p.conditions(idx,2) = d_offsets(dd);
        idx=idx+1;
    end
end
clear idx;

% shuffle these
p.rnd_idx = randperm(p.ntrials);
p.conditions = p.conditions(p.rnd_idx,:);


% generate list of WM, distractor positions
qtmp = floor(4*rand(p.ntrials,1)); % first, pick a quadrant (any quadrant)
atmp = p.merid_sep + (90-2*p.merid_sep)*rand(p.ntrials,1); % randomize pos w/in quadrant
% above in polar ang

p.wm_ang = qtmp*90+atmp; clear qtmp atmp;

p.dist_ang = nan(p.ntrials,1); % set up distractor angles
p.jitter_amt = nan(p.ntrials,1); % save out jitter amount too (CARTESIAN)
dsttmp = p.conditions(:,1)==2; % to clean up indexing below

p.jitter_amt(dsttmp) = p.dist_jitter*2*rand(sum(dsttmp),1)-p.dist_jitter;

p.dist_ang(dsttmp) = p.wm_ang(dsttmp) + p.conditions(dsttmp,2) * p.dist_step + p.jitter_amt(dsttmp);

% generate distractor directions
p.dist_dir = nan(p.ntrials,1);
p.dist_dir(dsttmp) = 1+round(rand(sum(dsttmp),1)); % 1 = CCW = 180; 2 = CW = 0


clear dsttmp;



%%
%p.TR = 3;
 
% ------ timing of trial events --------- %
% ITIs should be int_num * TR + 0.9 + 0.9 so that beginning of delay is
% locked to TR
 
% THINGS THAT HAPPEN BEFORE DELAY
p.cue_dur  = 1.0; % pre-cue: whether it's a distractor trial or not
p.targ_dur = 0.5; % WM target (every trial)


% WE WANT DELAY ONSET TO LOCK TO A TR
p.delay1_dur = 4.5; % first delay: pre-distractor
p.distractor_dur = 1.0; % distractor duration (or blank-screen duration)
p.delay2_dur = 6.5; % second delay: post-distractor
 
% THINGS THAT HAPPEN AFTER DELAY (end of delay also locked to TR)
p.resp_dur = 0.8; % MGS period (time to make MGS response) LOOK UP FROM MASIH's
p.MGS_feedback_dur = 1.0; % MGS feedback stimulus presentation
p.post_feedback_dur = 0.7; % time after MGS feedback before distractor feedback
p.distractor_feedback_dur = 1.0; % duration of distractor feedback stimulus (fix point change)



% NOT included in trial time:
p.distractor_resp_window = 2.5; % can respond up to this amt after distractor onset

% try to be in range of ~ 6-12 s (avg 9)


% short/medium/long ITIs...

ITI_TRs = [1*ones(floor(p.ntrials/4),1); 2*ones(2*ceil(p.ntrials/4),1); 3*ones(floor(p.ntrials/4),1)];


p.itis = p.TR-(p.targ_dur+p.cue_dur) + ...
         p.TR*ITI_TRs  + ...
     mod(p.TR-(p.resp_dur+p.MGS_feedback_dur+p.post_feedback_dur+p.distractor_feedback_dur),p.TR);

p.itis = p.itis(randperm(length(p.itis)));     

% -------- timing of experiment events ------- %
p.start_wait = 1 * p.TR + p.TR-(p.cue_dur+p.targ_dur);%2.4 + 0.9; % after first trigger [after dummys]
p.end_wait = 1*p.TR+mod(p.TR-(p.resp_dur+p.MGS_feedback_dur+p.post_feedback_dur+p.distractor_feedback_dur+p.itis(1)),p.TR);

p.trial_dur = p.cue_dur + p.targ_dur + p.delay1_dur + p.distractor_dur + p.delay2_dur + ...
              p.resp_dur + p.MGS_feedback_dur + p.post_feedback_dur + p.distractor_feedback_dur + p.itis;

p.exp_dur = p.start_wait + p.end_wait + sum(p.trial_dur);


% ------- things to save -------- %
p.wm_coords = p.wm_ecc .* [cosd(p.wm_ang) sind(p.wm_ang)];
p.dist_coords = p.wm_ecc .* [cosd(p.dist_ang) sind(p.dist_ang)];

% being a bit redundant here...
p.trial_start = nan(p.ntrials,1);
p.targ_start  = nan(p.ntrials,1);
p.delay_start = nan(p.ntrials,1);
p.dist_start  = nan(p.ntrials,1);
p.MGS_start   = nan(p.ntrials,1);
p.MGS_feedback_start = nan(p.ntrials,1);
p.iti_start   = nan(p.ntrials,1);
p.trial_end   = nan(p.ntrials,1);

p.behind_by = nan(p.ntrials,1); % keep track of timing errors, they should be basically tiny


% ------- keyboard stuff --------------------------- %
if ismac == 1
    p.esc_key = KbName('escape'); % press this key to abort
else
    try
        p.esc_key = KbName('esc');
    catch
        p.esc_key = KbName('escape');
    end
end
p.start_key = [KbName('5%') KbName('5')];  % should be lower-case t at prisma? (or %5, which is top-row 5, or 5, which is numpad 5)
p.space = KbName('space');
p.resp_keys = [KbName('1!') KbName('2@')]; % LEFT = CCW, RIGHT = CW




% ------- Screen setup, optics --------- %

% TODO: be a bit more clever here so that we can use this for practice
% outside scanner too....
if p.scanner == 0
    Screen('Preference', 'SkipSyncTests', 1) %changed from 21 on 3/12 
    p.desired_resolution = [3000 2000]; % desired resolution, compare to the 'actual' resolution below [this is for mbp laptop]
    p.desired_refresh_rate = 120;  % change this to 120 313

else
    p.desired_resolution = [1280 1024]; % scanner
    p.desired_refresh_rate = 120; %changed from 120  on 312
end

FlushEvents('KbName')
% TODO: assert correct refresh rate, resolution; bail otherwise

tmp_res = Screen('Resolution',max(Screen('Screens')));
p.resolution = [tmp_res.width tmp_res.height];
p.refresh_rate = tmp_res.hz;
p.screen_height = 36.3; % cm

if sum(p.resolution==p.desired_resolution)~=2 || p.refresh_rate~=p.desired_refresh_rate
    sprintf('Unexpected resolution/RR: expected %i, %i, %i; found %i, %i, %i',p.desired_resolution(1),p.desired_resolution(2),...
        p.desired_refresh_rate,p.resolution(1),p.resolution(2),p.refresh_rate)
    error('spDist_scanner:displayError','Unexpected resolution/RR: expected %i, %i, %i; found %i, %i, %i',p.desired_resolution(1),p.desired_resolution(2),...
        p.desired_refresh_rate,p.resolution(1),p.resolution(2),p.refresh_rate);
end


p.screen_width = p.screen_height * p.resolution(1)/p.resolution(2); %
p.viewing_distance = 63; % cm
p.screen_height_deg = 2*atan2d(p.screen_height/2,p.viewing_distance);
p.screen_width_deg  = 2*atan2d(p.screen_width/2, p.viewing_distance);
p.ppd = p.resolution(2)/p.screen_height_deg;  % used to convert rects, positions later on

p.center = p.resolution/2;  % could do offset centers, etc?

[w, p.scr_rect] = Screen('OpenWindow',max(Screen('Screens')),[0 0 0]); HideCursor;
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference','TextRenderer',1);

p.ifi = Screen('GetFlipInterval',w);

% ------------- generate distractor dot sequences ------------------ %
% note: we're doing all unit conversions here - so above definitions are in
% sec, deg, etc

p.dist_ndots = round(p.dist_rad^2*pi*p.dist_density);

p.dist_dot_seq = cell(p.ntrials,1);
for tt = 1:p.ntrials
    if p.conditions(tt,1)==2 % if a distractor trial, make a dot sequence
        
        p.dist_dot_seq{tt} = make_dot_seq_rot(p.dist_ndots, ...
                   p.dist_dir_proto(p.dist_dir(tt)), ...
                   p.distractor_dur*p.refresh_rate, ...
                   p.dist_speed/(p.refresh_rate*p.distractor_dur),...
                   p.task_coh,...
                   round(p.dist_dotlife*p.refresh_rate));
        
    end
end

% ------------------ generate rects we use later ---------------------- %
p.aperture_rect = CenterRectOnPoint([0 0 2 2] * p.ppd * p.aperture_size,p.center(1),p.center(2));
p.fix_rect_out  = CenterRectOnPoint([0 0 2 2] * p.ppd * p.fix_size_out, p.center(1),p.center(2));
p.fix_rect_in   = CenterRectOnPoint([0 0 2 2] * p.ppd * p.fix_size_in,  p.center(1),p.center(2));

% store onset of each frame of distractor
p.dist_fr_time = nan(p.ntrials,p.distractor_dur*p.refresh_rate);

% store distractor responses
p.dist_resp = nan(p.ntrials,1);
p.dist_RT = nan(p.ntrials,1);


% create black/white colors
dist_dotcolors_idx = repmat(1:size(p.dist_dotcolor,1),1,ceil(p.dist_ndots/size(p.dist_dotcolor,1)));
dist_dotcolors_idx = dist_dotcolors_idx(1:p.dist_ndots);

dist_dotcolors = nan(3,length(dist_dotcolors_idx));
for cc = 1:size(p.dist_dotcolor,1)
    dist_dotcolors(:,dist_dotcolors_idx==cc) = repmat(p.dist_dotcolor(cc,:).',1,sum(dist_dotcolors_idx==cc));
end

% --------- eyetracking ----------- %
if p.do_et == 1
    
    if p.scanner == 1
        Eyelink('SetAddress','192.168.1.5')
    end
    
    el=EyelinkInitDefaults(w);
    
    el.backgroundcolour=p.bg_color(1);  % TODO: fix this?
    el.calibrationtargetcolour=p.fix_color(1);
    el.calibrationtargetsize=2*p.wm_size;
    el.calibrationtargetwidth=1;
    el.msgfontcolour=p.fix_color(1);
    p.foregroundcolour=p.fix_color(1);
     % 192.168.1.5
    % sca
    EyelinkUpdateDefaults(el);

    
    Eyelink('Initialize','PsychEyelinkDispatchCallback') % initialises the eyetracker
   % SCANNER: right eye!!!!!!
    Eyelink('command','calibration_type=HV13'); % updating number of callibration dots
    s=Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,AREA');% (,GAZERES,HREF,PUPIL,STATUS,INPUT');
    s=Eyelink('command', 'sample_rate = 500');
    s=Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
    s=Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
    
    
    
    % make sure that we get gaze data from the Eyelink
    
    
    
    %------ calibrate the eye tracker --------
    EyelinkDoTrackerSetup(el);
    if s~=0
        error('link_sample_data error, status: ',s)
    end
    Eyelink('openfile',p.eyedatafile);

end

draw_aperture = @() Screen('FillOval',w,p.bg_color,p.aperture_rect);

Screen('FillRect',w,[0 0 0]);
draw_aperture();

Screen('DrawDots',w,[0;0],p.fix_size_mult*p.fix_size_out*p.ppd*2+p.fix_pen,p.dim_amt*p.fix_color,p.center,2);
Screen('DrawDots',w,[0;0],p.fix_size_mult*p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2);

Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.dim_amt*p.fix_color,p.center,2);
Screen('Flip',w);


% check for esc, space.... 

resp = 0;
while resp == 0
    [resp, ~] = checkForResp(p.start_key, p.esc_key);
    if resp == -1
        Screen('CloseAll'); ShowCursor;
        if p.do_et == 1
            Eyelink('ShutDown');
        end
        return;
    end
end
clear resp;
p.expt_start = GetSecs;


% blank screen - fix back to normal size
Screen('FillRect',w,[0 0 0]);
draw_aperture();
Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.dim_amt*p.fix_color,p.center,2); 
Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.dim_amt*p.fix_color,p.center,2); 

Screen('Flip',w);

if p.do_et == 1
    Eyelink('Message','xDAT %i', 49);
    Eyelink('StartRecording'); % make 1 big edf file (save time)
end


% ------ initial wait time (check for esc) ---------

resp = 0;
while (GetSecs-p.expt_start) < p.start_wait
    [resp, ~] = checkForResp([], p.esc_key);
    if resp == -1
        sca;
        fprintf('ESC pressesd during post-trigger wait time\n');
        ShowCursor;
        if p.do_et == 1
            Eyelink('ShutDown');
        end
        return;
    end
end
clear resp;

for tt = 1:p.ntrials
    
    
    % distractor cue (XDAT 1) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % NOTE: I'm making this a "trial should have started" marker, to keep
    % timing from accumulating errors & to keep below code more readable
    trial_start = p.expt_start + p.start_wait + sum(p.trial_dur(1:(tt-1)));
    p.trial_start(tt) = GetSecs;
    p.behind_by(tt) = p.trial_start(tt)-trial_start;
    
    if p.do_et == 1

        Eyelink('Message','TarX %s', num2str(0));
        Eyelink('Message','TarY %s', num2str(0));
        
        
        Eyelink('Message','xDAT %i',1);
        
        Eyelink('command', 'record_status_message "TRIAL %d of %d"', tt, p.ntrials);
        
    end

    
    while GetSecs < trial_start + p.cue_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        
        % fixation
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2);         
        % cue (filled fix)
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.cue_colors(p.conditions(tt,1),:),p.center,2);         
        
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    
    
    % targets (XDAT 2) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    p.targ_start(tt) = GetSecs;
    
    
    if p.do_et == 1
        
        Eyelink('Message','xDAT %i',2);
        
    end
    
    
    
    while GetSecs < trial_start + p.cue_dur + p.targ_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        % target 1
        Screen('DrawDots',w,p.ppd*[1;-1].*p.wm_coords(tt,:).', p.wm_size*p.ppd, p.wm_color, p.center, 2);
        
        
        
        % fixation
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        if isnan(p.targ_start(tt))
            p.targ_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    
    % Delay 1 (XDAT 3) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    if p.do_et == 1
        Eyelink('Message','xDAT %i',3);    
    end
    
    
    while GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay1_dur
        
        % draw aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        if isnan(p.delay_start(tt)) % only do this once
            p.delay_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    
    % DISTRACTOR (or not) (XDAT 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if p.do_et == 1
        Eyelink('Message','xDAT %i',4);
    end
    
    % {TODO}
    % if p.conditions(tt,1) == 1, just hang out
    % otherwise, show our distractor stimulus & poll keyboard for response
    if p.conditions(tt,1) == 1 % if it's a no-distractor trial
       
        while GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay1_dur + p.distractor_dur 
            
            % draw aperture
            Screen('FillRect',w,[0 0 0]);
            draw_aperture();
            
            % draw fixation
            Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2);
            Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2);
            Screen('Flip',w);
            
            if isnan(p.dist_start(tt)) % only do this once
                p.dist_start(tt) = GetSecs;
            end
            
            % check for esc....
            [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
            if resp == -1
                Screen('CloseAll'); ShowCursor;
                Eyelink('StopRecording');
                Eyelink('ShutDown');
                save(p.filename,'p');
                return;
            end
            
            
        end
        
        
    elseif p.conditions(tt,1) == 2 % if there's a distractor, draw it!
        
        
        % TODO: store frame timing!
        % loop over frames
        for ff = 1:size(p.dist_dot_seq{tt},3)
            
            % draw aperture
            Screen('FillRect',w,[0 0 0]);
            draw_aperture();
            
            % draw fixation
            Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2);
            Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2);
                        
            % draw distractor
            
            % first start w/ a proto-distractor
            % Screen('DrawDots',w,p.dist_coords(tt,:).*[1 -1]*p.ppd,p.fix_size_out*p.ppd*2,[255 0 0],p.center,2);
            Screen('DrawDots',w,p.dist_dot_seq{tt}(:,:,ff).*[1;-1]*p.ppd,p.dist_dotsize*p.ppd,dist_dotcolors,p.dist_coords(tt,:).*[1 -1]*p.ppd+p.center,2);
            
            % Flip screen
            Screen('Flip',w);
            p.dist_fr_time(tt,ff) = GetSecs;
            
            % TODO: if first frame, save start_resp
            if ff == 1
                % use this to decide if delay-2 response is in time
                resp_start_time = p.dist_fr_time(tt,1);
                
                p.dist_start(tt) = p.dist_fr_time(tt,1);
                
            end
            
            % check for escape & response
            % check for esc....
            [resp] = checkForResp(p.resp_keys, p.esc_key); % TODO: maybe turn on/off gaze indicator?
            if resp == -1
                Screen('CloseAll'); ShowCursor;
                Eyelink('StopRecording');
                Eyelink('ShutDown');
                save(p.filename,'p');
                return;
            elseif resp~=0
                % handle responses (TODO)
                if isnan(p.dist_resp(tt))
                    p.dist_resp(tt) = find(p.resp_keys==resp);
                    p.dist_RT(tt) = GetSecs-resp_start_time;
                end
                
            end
            
            
            
        end
        
    end
    
    % Delay 2 (XDAT 5) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if p.do_et == 1
        Eyelink('Message','xDAT %i',5);
    end
    
    while GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay1_dur + p.distractor_dur + p.delay2_dur
        
        % draw aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2);
        Screen('Flip',w);
        
        if isnan(p.delay_start(tt)) % only do this once
            p.delay_start(tt) = GetSecs;
        end
        
        
        
        % always check for response
        [resp] = checkForResp(p.resp_keys, p.esc_key); % TODO: maybe turn on/off gaze indicator?
        
        
        % if -1, do escape stuff
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        elseif resp ~= 0 && GetSecs < (resp_start_time+p.distractor_resp_window)
            % try to record responsees
            % handle responses (TODO)
            if isnan(p.dist_resp(tt))
                p.dist_resp(tt) = find(p.resp_keys==resp);
                p.dist_RT(tt) = GetSecs-resp_start_time;
            end
            
        end

        
    end

    clear resp_start_time;
    
    % MGS response (XDAT 6) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if p.do_et == 1
        Eyelink('Message','xDAT %i',6);
    end

    
    
    
    while GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay1_dur + p.distractor_dur + p.delay2_dur + p.resp_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

        Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.go_color,p.center,2);
        %Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.bg_color,p.center,2); 

        Screen('Flip',w);
        
        if isnan(p.MGS_start(tt))
            p.MGS_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
    end
    
    % feedback (XDAT 7, tarx, tary) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if p.do_et == 1
        Eyelink('Message','TarX %s', num2str(p.wm_coords(tt,1)));
        Eyelink('Message','TarY %s', num2str(p.wm_coords(tt,2)));
        Eyelink('Message','xDAT %i',7);
        
    end
    
    while GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay1_dur + p.distractor_dur + p.delay2_dur + p.resp_dur + p.MGS_feedback_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        % target 1
        Screen('DrawDots',w,p.ppd*[1;-1].*p.wm_coords(tt,:).', p.wm_size*p.ppd, p.wm_color, p.center, 2);
        
        % fixation
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

        Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.go_color,p.center,2);
        %Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.bg_color,p.center,2); 
        
        %Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 

        Screen('Flip',w);
        
        if isnan(p.MGS_feedback_start(tt))
            p.MGS_feedback_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
 
    end
    
    % Post-feedback wait (XDAT 8) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % (wait for subj to return to fixation, then show behavioral feedback
    % on next epoch)
    
    if p.do_et == 1
        Eyelink('Message','xDAT %i',8);    
        Eyelink('Message','TarX %s', num2str(0));
        Eyelink('Message','TarY %s', num2str(0));
        
    end
    
    
    while GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay1_dur + p.distractor_dur + p.delay2_dur + p.resp_dur + p.MGS_feedback_dur + p.post_feedback_dur
        
        % draw aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        if isnan(p.delay_start(tt)) % only do this once
            p.delay_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    
    % Distractor response feedback (XDAT 9) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if p.do_et == 1
        Eyelink('Message','xDAT %i',9);    
    end
    
    
    while GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay1_dur + p.distractor_dur + p.delay2_dur + p.resp_dur + p.MGS_feedback_dur + p.post_feedback_dur + p.distractor_feedback_dur
        
        % draw aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        
        % draw response feedback (if distractor trial)
        if p.conditions(tt,1) == 2
            if isnan(p.dist_resp(tt))
                % draw yellow (feedback_colors 3)
                this_color = p.dist_feedback_colors(3,:);
            else
                if p.dist_resp(tt) == p.dist_dir(tt)
                    % draw green (1)
                    this_color = p.dist_feedback_colors(1,:);
                else
                    this_color = p.dist_feedback_colors(2,:);
                end
            end
            Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,this_color,p.center,2);
            clear this_color;
        end
        
        
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        if isnan(p.delay_start(tt)) % only do this once
            p.delay_start(tt) = GetSecs;
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    
    % ITI (XDAT 10) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if p.do_et == 1

        Eyelink('Message','xDAT %i',10);
    end

    % TODO: base this on sum of all events up til then? or do that
    % everywhere?
    while GetSecs < trial_start + p.cue_dur + p.targ_dur + p.delay1_dur + p.distractor_dur + p.delay2_dur + p.resp_dur + p.MGS_feedback_dur + p.post_feedback_dur + p.distractor_feedback_dur + p.itis(tt)
        
        % draw aperture
        Screen('FillRect',w,[0 0 0]);
        draw_aperture();
        
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.dim_amt*p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.dim_amt*p.fix_color,p.center,2); 
        
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        
        Screen('Flip',w);
        
        if isnan(p.iti_start(tt))
            p.iti_start(tt) = GetSecs;
            % save [note: in scanner, do this at beginning of ITI after first flip]
            save(p.filename,'p');
        end
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if p.do_et==1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    
    

    
end

% ------- wait for p.end_wait -------- %
end_tmp = GetSecs;

resp = 0;
while (GetSecs-end_tmp) < p.end_wait
    [resp, ~] = checkForResp([], p.esc_key);
    if resp == -1
        sca;
        fprintf('ESC pressesd during post-experiment wait time\n');
        ShowCursor;
        return;
    end
end
clear resp;
p.end_expt = GetSecs;

p.correct = p.dist_resp==p.dist_dir;
p.n_misses = sum(isnan(p.dist_resp(p.conditions(:,1)==2)));
p.acc = nanmean(p.correct(~isnan(p.dist_resp)));

save(p.filename,'p');




% END OF EXPERIMENT - TEXT
if p.do_et == 1
    Eyelink('Message','xDAT %i',50);
end

Screen('FillRect',w,[0 0 0]);
draw_aperture();
txt = sprintf('End of run %i',p.run);
txt2 = sprintf('Accuracy: %0.02f%%, %i misses',p.acc*100,p.n_misses);
DrawFormattedText(w,txt, 'center',p.center(2)-5*p.ppd,p.fix_color);
DrawFormattedText(w,txt2,'center',p.center(2)-2*p.ppd,p.fix_color);

Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.dim_amt*p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.dim_amt*p.fix_color,p.center,2); 

Screen('Flip',w);

if p.do_et == 1
    Eyelink('StopRecording');
    Eyelink('ReceiveFile',[p.eyedatafile '.edf'],[p.eyedatafile '.edf']);
    
    p.eyedatafile_renamed = [p.filename(1:(end-3)) 'edf'];
    movefile([p.eyedatafile '.edf'],p.eyedatafile_renamed);
    
    Eyelink('ShutDown');
end

resp = 0;
while resp == 0
    [resp, ~] = checkForResp(p.space, p.esc_key);
end
clear resp;



Screen('CloseAll');
ShowCursor;
catch
   Screen('CloseAll');
   ShowCursor;
end


return