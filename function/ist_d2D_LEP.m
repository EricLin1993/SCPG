function X = ist_d2D_LEP(y, mask, options)

% INPUT
% y - 2D NUS measurement vectorization (time domain, dimensions M x 1)
% mask - sampling scheme (logical)


% PARAMETRES
% options.res ∈(0; ...) - residue norm (should be equal to noise norm)
% options.threshold ∈(0; 1) - relative to the highest peak in the spectrum with artifacts
% options.max_iter - limit of iterations

% OUTPUT
% X (frequency domain)

% POSSIBLE MODIFICATIONS:
% linear change of the threshold from iteration to iteration: options.change_thr = 1
%==============================================================================
% Enping Lin
% enpinglin2022@126.com
% 20220530

%%%%%%%%%%%% parameters %%%%%%%%%%%%
mask = logical(mask);


res = options.res;

if isfield(options, 'max_iter')
    max_iter = options.max_iter;
else
    max_iter = 500;
end

if isfield(options, 'threshold')
    threshold = options.threshold;
else
    threshold = 0.99;
end

%%%%%%%%%%%% modifications %%%%%%%%%%%%
if isfield(options, 'change_thr')
    change_thr = options.change_thr; % = 1
else
    change_thr = 0;
end

%%%%%%%%%%%% initialization %%%%%%%%%%%%

Temp = zeros(size(mask));
Temp(mask) = y;
sp_art = fft2(Temp)/sqrt(prod(size(Temp)));
 % spectrum with artifacts
t = threshold*max(abs(sp_art(:)));
X = zeros(size(mask));
r = y;

%%%%%%%%%%%% main loop %%%%%%%%%%%%
for i = 1:max_iter
%     disp(horzcat('iteration', ' ', int2str(i)))
    if norm(r) >= res
        i
        X = X + delta(sp_art, t);
        Temp = ifft2(X)*sqrt(prod(size(X)));
        Temp = Temp(mask);
        r = y(:) - Temp;
        fX = zeros(size(mask));
        fX(mask) = r;
        sp_art = fft2(fX)/sqrt(prod(size(fX)));
        t = threshold*max(abs(sp_art(:)));
        if change_thr == 1
            t = t*threshold*(max_iter - i)/max_iter;
        end
    else
        disp('stopping criterion has been reached')
        break
    end
end
end
