function DotMotion(motion, direction, density, time, rand_percent)
% DotMotion([window, rect, fps, screenNumber, doublebuffer], motion, density, direction, time, rand_percent, draw_mask)
% dot motion using Screen('DrawDots') subfunction from PsychToolbox 3
% 3 different motions: rotational motion = 'rot'
%                      radial motion = 'rad'
%                      translational motion = 'trans'
% density: number of dots per degree^2
% direction: 'cw' or 'ccw' for rotational motion
%            'in' or 'out' for radial motion
%            'left', 'right', 'up' or 'down' for translational motion
% time: on-screen time in millisecond
% rand_percent: percentage of points going in random directions (0 to 1)
% 
% Author: Paul Fricker, CerCo-UMR5549
% Contact: paul.fricker@cnrs.fr
% Made with the PTB3 toolbox on Matlab R2017a
% Some functions used may not work in older versions


% Screen('Preference', 'SkipSyncTests', 0); 
Screen('Preference', 'VisualDebugLevel', 0);

InitializeMatlabOpenGL;
AssertOpenGL;

if nargin < 6
    rand_percent = [];
end

if isempty(rand_percent)
    rand_percent = 0;
end

try
    %% Set dot field parameters
    mon_width   = 497;      % horizontal dimension of viewable screen/monitor (cm)
    v_dist      = 180;      % viewing distance from monitor (cm)
    dot_speed   = 7;        % dot speed (deg/sec)
    max_d       = 100;      % maximum radius of annulus (degrees)
    min_d       = 0;        % minumum
    dot_w       = 0.20;     % width of dot (deg)
    life_limit  = 200;      % lifetime of each dot (ms)
    waitframes  = 1;        % Show new dot-images at each waitframes monitor refresh.

    % Bool flag for motion selection
    switch motion
        case 'rot'
            rot_flag   = 1;
            rad_flag   = 0;
            trans_flag = 0;
        case 'rad'
            rot_flag   = 0;
            rad_flag   = 1;
            trans_flag = 0;
        case 'trans'
            rot_flag   = 0;
            rad_flag   = 0;
            trans_flag = 1;
        otherwise
            rot_flag   = 0;
            rad_flag   = 0;
            trans_flag = 0;
    end

    %% Open the screen

    doublebuffer=1; % ensure a good display by computing the wanted data first and then displaying it
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    [w, rect] = Screen('OpenWindow', screenNumber, 0,[], 32, doublebuffer+1);
    
    % Enable alpha blending with proper blend-function. We need it
    % for drawing of smoothed points:
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    [center(1), center(2)] = RectCenter(rect);
    fps=Screen('FrameRate',w);          % frames per second
    ifi=Screen('GetFlipInterval', w);   % inter-frame interval (1/fps)
    if fps==0
        fps=1/ifi;
    end
    
    HideCursor;	% Hide the mouse cursor
    Priority(MaxPriority(w));

    % Do initial flip
    vbl=Screen('Flip', w);


    %% Initialize dot positions and velocities

    ifi = 1/fps;                                % inter-frame interval
    [center(1), center(2)] = RectCenter(rect);  % Center coordinates of the screen
    white = WhiteIndex(w);                      % Define color white
    nframes = (time * 1e-3) * fps + 20;         % number of frames to draw given the wanted on-screen display time +20 frames for the randomness to be properly set

    v_angle = 2 * atan((mon_width/2)/v_dist);                                           % vision angle formula
    ppd = ((rect(3)-rect(1)) / v_angle) * (pi/180);                                     % pixels per degree
    pfs = dot_speed * ppd / fps;                                                        % dot speed (pixels/frame)
    s = dot_w * ppd;                                                                    % dot size (pixels)
    r_speed_factor = 0.3 / ppd;                                                         % factor for radius depending velocities

    ndots =  floor(density * ((rect(3) - rect(1)) * (rect(4) - rect(2)) * (1/ppd)^2));  % number of dots from density(dots/degree^2)
    random_dir = zeros(ndots, 1);                                                       % vector defining if a dot will have a random direction (1) or not (0)
    rand_n = ndots * rand_percent;                                                      % number of dots having a random direction
    random_dir(randperm(numel(random_dir), floor(rand_n))) = 1;                         % update vector

    rmax = max_d * ppd;                                                                 % maximum radius of annulus (pixels from center)
    rmin = min_d * ppd;                                                                 % minimum
    r = rmax * sqrt(rand(ndots,1));                                                     % position from the center
    r(r<rmin) = rmax;                                                                   % redraw when in the dead zone
    t = 2*pi*rand(ndots,1);                                                             % theta polar coordinate
    cs = [cos(t), sin(t)];                                                              % cosinus and sinus for coordinates conversion
    xy = [r r] .* cs;                                                                   % dot positions in Cartesian coordinates (pixels from center)
    lifetimes = floor(life_limit * rand(ndots, 1));

    % motion direction (in/ccw (-) or out/cw (+)) for each dot
    if strcmp(direction ,'in') || strcmp(direction, 'ccw')
        mdir = 2 * floor(rand(ndots,1)) - 1;
    end
    if strcmp(direction ,'out') || strcmp(direction, 'cw')
        mdir = 2 * floor(rand(ndots,1)) + 1;
    end
    if strcmp(direction, 'in')
        r_speed_factor = r_speed_factor * 0.7;  % Adjustement for contraction motion 
    end

    % movement update for each dot
    if rot_flag
        dt = pfs * mdir;                                    % change in theta per frame (radians)
        dr = pfs * mdir ./ dt;      
        mdir_x = -mdir;
        mdir_y = mdir;
        dxdy = [dr dr] .* [pfs*mdir_x(1)*sin(t), pfs*mdir_y(1)*cos(t)];
    end
    if rad_flag
        dt = pfs * mdir .* 0;                               % initial null change of theta for random movement
        dr = pfs * mdir .* (r_speed_factor * r);            % change in radius per frame (pixels)
        dxdy = [dr dr] .* cs;                               % change in x and y per frame (pixels)
    end
    if trans_flag
        x = ((rect(3) - rect(1)) * ppd) * rand(ndots,1);    % change to Cartesian coordinates for movement
        y = ((rect(4) - rect(2)) * ppd) * rand(ndots,1);
        xy = [x y] - center;                                % translation of the drawing zone (origin at [0 0])
        switch direction
            case 'left'
                mdir_x = 2 * floor(rand(ndots,1)) - 1;      % left (-) or right (+)
                mdir_y = 2 * floor(rand(ndots,1));          % up (-) or down (+)
            case 'right'
                mdir_x = 2 * floor(rand(ndots,1)) + 1;
                mdir_y = 2 * floor(rand(ndots,1));
            case 'up'
                mdir_x = 2 * floor(rand(ndots,1));
                mdir_y = 2 * floor(rand(ndots,1)) - 1;
            case 'down'
                mdir_x = 2 * floor(rand(ndots,1));
                mdir_y = 2 * floor(rand(ndots,1)) + 1;
            otherwise
                mdir_x = 2 * floor(rand(ndots,1));
                mdir_y = 2 * floor(rand(ndots,1));
        end
        dxdy = [pfs * mdir_x pfs * mdir_y];
    end

    % Clamp point sizes to range supported by graphics hardware:
    [minsmooth,maxsmooth] = Screen('DrawDots', w);
    s = min(max(s, minsmooth), maxsmooth);

    %% Animation loop
    for i = 1:nframes
        if (i>1)
            lifetimes = lifetimes + (ifi * 1e3);
            expansion_point = [mx my]; % origin of field of expansion
            Screen('DrawDots', w, xymatrix, s, white, expansion_point, 1);  % draw white dots at any x and y from xymatrix
            
            if (i<20)
                Screen('FillRect', w, [0 0 0], rect);  % Black screen for the first 20 frames to let the random process to be fully active
            end
            Screen('DrawLines', w, [center(1)-.5*ppd center(1)+.5*ppd center(1) center(1); 
                                    center(2) center(2) center(2)-.5*ppd center(2)+.5*ppd],...
                                    3, [255 0 0]);  % Center cross
            Screen('DrawLines', w, [expansion_point(1)-.5*ppd expansion_point(1)+.5*ppd expansion_point(1) expansion_point(1);
                                    expansion_point(2) expansion_point(2) expansion_point(2)-.5*ppd expansion_point(2)+.5*ppd],...
                                    1, [0 255 0]);  % Fixation cross
            Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
        end
        
        [mx, my, buttons]=GetMouse(screenNumber); % monitor the position of the mouse and its button click
        if KbCheck || any(buttons) % break out of loop when any key or button is pressed
            break;
        end

        % movement update
        if rot_flag
            r = r + dr;
            t = t + dt;                         % update theta
            xy = xy + dxdy;                     % compute new positions
        end
        if rad_flag
            r = r + dr;							% update polar coordinates too
            t = t + dt;
            xy = xy + dxdy;						% move dots

        end
        if trans_flag
           xy = xy + dxdy;                      % move dots
        end
        
        % check which dots will randomly move
        r_rand = find(random_dir == 1);
        nrand = length(r_rand);
        rand_mov_x = rand(ndots, 1);
        rand_mov_y = rand(ndots, 1);

        if nrand
            random_dir(r_rand) = 0;
            if rot_flag
                mdir_x(r_rand) = 2*rand_mov_x(r_rand)-1;
                mdir_y(r_rand) = 2*rand_mov_y(r_rand)-1;
                dt(r_rand) = pfs * ((mdir_x(r_rand) + mdir_y(r_rand)) / 2);
                dr(r_rand) = pfs * ((mdir_x(r_rand) + mdir_y(r_rand)) / 2) ./ dt(r_rand);
                dxdy(r_rand,:) = [dr(r_rand) dr(r_rand)] .* [pfs*mdir_x(r_rand).*sin(t(r_rand)) pfs*mdir_y(r_rand).*cos(t(r_rand))];
            end
            if rad_flag
                mdir(r_rand) = 2*rand_mov_x(r_rand)-1;
                dt(r_rand) = pfs * (2*rand_mov_y(r_rand)-1).* (r_speed_factor*r(r_rand));
                dxdy(r_rand,:) = [dr(r_rand) dr(r_rand)] .* [cos(t(r_rand) + dt(r_rand)) sin(t(r_rand) + dt(r_rand))];
            end
            if trans_flag
               mdir_x(r_rand) = 2*rand_mov_x(r_rand)-1;
               mdir_y(r_rand) = 2*rand_mov_y(r_rand)-1;
            end
        end
        
        % check to see which dots have gone beyond the borders of the
        % drawing zone or reached their life end
        r_out = find(r > rmax | r < rmin | lifetimes > life_limit);	 % dots to reposition
        nout = length(r_out);
        
        if nout
            % choose new coordinates
            r(r_out) = rmax * sqrt(rand(nout,1));
            r(r<rmin) = rmax;
            t(r_out) = 2*pi*(rand(nout,1));
            
            % now convert the polar coordinates to Cartesian
            cs(r_out,:) = [cos(t(r_out)), sin(t(r_out))];
            xy(r_out,:) = [r(r_out) r(r_out)] .* cs(r_out,:);
            
            % reset lifetime
            lifetimes(r_out) = 0;
            
            % compute the new cartesian velocities
            if rot_flag
                dxdy(r_out,:) = [dr(r_out).* (r_speed_factor*r(r_out)) dr(r_out).* (r_speed_factor*r(r_out))] .* [(pfs/4)*mdir_x(r_out).*sin(t(r_out)) (pfs/4)*mdir_y(r_out).*cos(t(r_out))];
            end
            if rad_flag
                dr(r_out) = (pfs/4) * mdir(r_out) .* (r_speed_factor*r(r_out));
                dxdy(r_out,:) = [dr(r_out) dr(r_out)] .* [cos(t(r_out) + dt(r_out)) sin(t(r_out) + dt(r_out))];
            end
            if trans_flag
               dxdy(r_out,:) = [pfs * mdir_x(r_out) pfs * mdir_y(r_out)];
            end
        end
        xymatrix = transpose(xy);

        if (doublebuffer==1)
            vbl=Screen('Flip', w, vbl + (waitframes-0.5)*ifi, 0);
        end
    end
    Priority(0);
    ShowCursor
    sca;
catch
    Priority(0);
    ShowCursor
    sca;
end
