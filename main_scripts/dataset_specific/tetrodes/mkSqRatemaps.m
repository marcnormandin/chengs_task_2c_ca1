% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [sq] = mkSqRatemaps(ratemaps)

v = homography_solve([20 0 0 20; 0 0 30 30], [20 0 0 20; 0 0 20 20]); % transform
x_edges = 0:1:20;
y_edges = 0:1:20;
kern = ratemaps.cfg.smoothingKernel;

for a = 1:length(ratemaps.data)
    sq(a).animal = ratemaps.data(a).animal;
    sq(a).group = ratemaps.data(a).group;
    for s = 1:length(ratemaps.data(a).session)
        sq(a).session(s).name = ratemaps.data(a).session(s).name;
        for tr = 1:length(ratemaps.data(a).session(s).trial)
            sq(a).session(s).trial(tr).trialNum = ratemaps.data(a).session(s).trial(tr).trialNum;
            sq(a).session(s).trial(tr).context = ratemaps.data(a).session(s).trial(tr).context;
            sq(a).session(s).trial(tr).dig = ratemaps.data(a).session(s).trial(tr).dig;
            pos = [ratemaps.data(a).session(s).trial(tr).x; ratemaps.data(a).session(s).trial(tr).y];
            xy = homography_transform(pos, v);
            tx = xy(1,:); ty = xy(2,:); t = ratemaps.data(a).session(s).trial(tr).t_sec;
            sq(a).session(s).trial(tr).t = t;
            sq(a).session(s).trial(tr).x = tx;
            sq(a).session(s).trial(tr).y = ty;
            
            tspeed = ratemaps.data(a).session(s).trial(tr).runSpeed;
            spInd = tspeed > 2;
            
            occ_binned = histcounts2(ty(spInd), tx(spInd), y_edges, x_edges) ;
            TimeMap = occ_binned .* median(diff(t));
            smoothedOcc = conv2(TimeMap, kern, 'same');
            no_occ_idx = find(smoothedOcc < 0.05); % NaN out bins never visited
            smoothedOcc(no_occ_idx) = nan;
            sq(a).session(s).trial(tr).label = ratemaps.data(a).session(s).trial(tr).label;
            
            for c = 1:length(sq(a).session(s).trial(tr).label)
                spk_t = ratemaps.data(a).session(s).trial(tr).spk_t{c};
                spk_x = interp1(t, tx, spk_t)';
                spk_y = interp1(t, ty, spk_t)';
                spk_speed = interp1(t, tspeed, spk_t);
                spk_ind = spk_speed > 2;
                
                if sum(spk_ind) < ratemaps.cfg.minSpike
                    fprintf('um theres still slow spikes');
                    continue;
                end
                
                SpikeMap = histcounts2(spk_y(spk_ind), spk_x(spk_ind), y_edges, x_edges);
                smoothedSpk = conv2(SpikeMap, kern, 'same');
                smoothedSpk(no_occ_idx) = nan;
                switch ratemaps.cfg.smoothProtocol
                    case 'before'
                        tc = smoothedSpk./smoothedOcc;
                        tc(~isfinite(tc)) = NaN;
                        map = tc;
                        sq(a).session(s).trial(tr).maps{c} = map;
                        sq(a).session(s).trial(tr).mfr{c} = ratemaps.data(a).session(s).trial(tr).mfr{c};
                        
                    case 'after'
                        tc = SpikeMap ./ TimeMap;
                        tc(~isfinite(tc)) = 0;
                        
                        map = conv2(tc, kern, 'same');
                        map(no_occ_idx) = nan;
                        sq(a).session(s).trial(tr).maps{c} = map;
                        sq(a).session(s).trial(tr).mfr{c} = ratemaps.data(a).session(s).trial(tr).mfr{c};
                end
            end
        end
    end
end




