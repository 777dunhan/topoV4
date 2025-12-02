clearvars -except heatmaps

% Generate grid coordinates for a 224x224 grid
xx = zeros(224*224, 1);
yy = zeros(224*224, 1);

n = 1;
for i = 1:224
    for j = 1:224
        xx(n) = i;
        yy(n) = j;
        n = n + 1;
    end
end

% Initialize the fit structure for 256 heatmaps
pnum = 512; % Number of heatmaps
init = zeros(pnum, 1);
Fit = struct('a0', init, 'a1', init, 'b1', init, 'b2', init, 'c1', init, 'c2', init, 'k', init, 'rsquare', init);
rsquares = zeros(pnum, 1);
% Define the 2D Gaussian model
g = fittype( @(a0,a1,b1,b2,c1,c2,k, x, y) a0+a1*exp(-0.5*(((x-b1)./c1).^2-k*(x-b1).*(y-b2)+((y-b2)./c2).^2)), 'independent', {'x', 'y'}, 'dependent', 'z' );
coeffnames(g)

% Fit the Gaussian for each heatmap
for pi = 1:pnum
    tmp1 = squeeze(heatmaps(pi, :, :)); % Extract heatmap of size (224, 224)
    if min(tmp1(:)) >= 0
        tmp2 = tmp1'; % Transpose the heatmap for fitting
        tmp3 = imgaussfilt(tmp1, 1); % Smooth the heatmap
        [~, I] = max(tmp3(:)); % Find the peak response
        [II, JJ] = ind2sub([224, 224], I); % Convert index to (row, col)

        zz = tmp2(:); % Flatten the heatmap into a vector
        
        if any(isnan(zz)) || any(isinf(zz)) || std(zz) == 0
            disp(['Skipping heatmap ', num2str(pi), ' due to invalid data.']);
            continue;
        end
        
        try
            % Fit the 2D Gaussian to the heatmap
            [fitobject, gof] = fit([xx, yy], zz, g, ...
                'StartPoint', [0, 1, II, JJ, 10, 10, 0]);
    
            % Store the fitted parameters and goodness of fit
            Fit.a0(pi) = fitobject.a0;
            Fit.a1(pi) = fitobject.a1;
            Fit.b1(pi) = fitobject.b1;
            Fit.b2(pi) = fitobject.b2;
            Fit.c1(pi) = fitobject.c1;
            Fit.c2(pi) = fitobject.c2;
            Fit.k(pi) = fitobject.k;
            Fit.rsquare(pi) = gof.rsquare;
            % Display progress
            disp(['Progress: ', num2str(pi / pnum * 100, '%.2f'), '% | R^2: ', num2str(gof.rsquare)]);
            rsquares(pi) = gof.rsquare;
        catch
            Fit.a0(pi) = 0.0;
            Fit.a1(pi) = 0.0;
            Fit.b1(pi) = 0.0;
            Fit.b2(pi) = 0.0;
            Fit.c1(pi) = 0.0;
            Fit.c2(pi) = 0.0;
            Fit.k(pi) = 0.0;
            Fit.rsquare(pi) = 0.0;
            fprintf('fitting error occurred');
        end
    end
end

% Save the fit results
% save Fit Fit

% a0 = Fit.a0; % Baseline offset
% a1 = Fit.a1; % Amplitude
% b1 = Fit.b1; % X-center
% b2 = Fit.b2; % Y-center
% c1 = Fit.c1; % X-standard deviation
% c2 = Fit.c2; % Y-standard deviation
% k = Fit.k;   % Correlation parameter

% Save each parameter into its own .mat file
% save('41a0.mat', 'a0');
% save('41a1.mat', 'a1');
% save('41b1.mat', 'b1');
% save('41b2.mat', 'b2');
% save('41c1.mat', 'c1');
% save('41c2.mat', 'c2');
% save('41k.mat', 'k');

%% Calculate param

center = [480,639];
LT = [405,560];
x0 = (center(1)-LT(1))/2;
y0 = (center(2)-LT(2))/2;

init = zeros(pnum,1);
param = struct('theta',init,'r',init,'a',init,'b',init);

for pi = 1:pnum
    if Fit.rsquare(pi)>0
        b1 = Fit.b1(pi);
        b2 = Fit.b2(pi);
        c1 = Fit.c1(pi);
        c2 = Fit.c2(pi);
        k = Fit.k(pi);
        xx = b1-x0;
        yy = b2-y0;
        [theta, r] = cart2pol(xx,yy);
        A = 1/c1^2;
        B = -k;
        C = 1/c2^2;
        aa = (2/(A+C-sqrt(B^2+(C-A)^2)))^0.5;
        bb = (2/(A+C+sqrt(B^2+(C-A)^2)))^0.5;
        param.theta(pi) = theta/3.1415926*180;
        param.r(pi) = r;
        param.a(pi) = aa;
        param.b(pi) = bb;
    end
end
%%
x = param.r/15;
y = (2*param.a.*param.b*3.1415926).^0.5/15;
scatter(x,y,'Marker','.')
%%
Sigma = zeros(pnum,2);
Sigma(:,1) = param.a;
Sigma(:,2) = param.b;
save Sigma Sigma -v7.3