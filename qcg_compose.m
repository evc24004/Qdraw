L = load(fullfile('output', 'glyph_HUSKY.mat'));
qv = L.qv;
pv = L.pv;
crop = @(H, ql, pl) flipud(H(abs(pv) <= pl, abs(qv) <= ql));
norm1 = @(H) max(H/max(H(:)) - 0.05, 0) / 0.95;
up4 = @(H) interp2(H, 2);

husky = up4(norm1(crop(L.H, 3.4, 3.4)));
letters = cell(1, 5);
names = {'U', 'C', 'O', 'N', 'N2'};
for i = 1:5
    L = load(fullfile('output', sprintf('glyph_%s.mat', names{i})));
    letters{i} = up4(norm1(crop(L.H, 2.9, 3.1)));
end

[lh, lw] = size(letters{1});
[hh, hw] = size(husky);
gap = 28;
marg = 48;
rowgap = 12;
W = 5*lw + 4*gap + 2*marg;
Hgt = marg + hh + rowgap + lh + marg;
canvas = zeros(Hgt, W);

x0 = round((W - hw)/2);
canvas(marg+(1:hh), x0+(1:hw)) = husky;
y0 = marg + hh + rowgap;
for i = 1:5
    x = marg + (i-1)*(lw+gap);
    canvas(y0+(1:lh), x+(1:lw)) = max(canvas(y0+(1:lh), x+(1:lw)), letters{i});
end

navy = [0 14 47]/255;
t = (linspace(0,1,256)').^0.85;
cmap = navy + t .* ([1 1 1] - navy);
imwrite(ind2rgb(1 + round(canvas*255), cmap), ...
    fullfile('output', 'uconn_quantum.png'));

anchors = [0.0015 0.0005 0.0139; 0.0873 0.0444 0.2244; ...
           0.2588 0.0386 0.4064; 0.4166 0.0903 0.4328; ...
           0.5783 0.1481 0.4044; 0.7355 0.2154 0.3297; ...
           0.8654 0.3168 0.2262; 0.9541 0.4433 0.1201; ...
           0.9877 0.6521 0.2113; 0.9647 0.8434 0.4321; ...
           0.9884 0.9984 0.6449];
cmapInf = interp1(linspace(0,1,size(anchors,1)), anchors, linspace(0,1,256));
imwrite(ind2rgb(1 + round(canvas*255), cmapInf), ...
    fullfile('output', 'uconn_quantum_inferno.png'));
fprintf('wrote uconn_quantum.png and uconn_quantum_inferno.png (%dx%d)\n', W, Hgt);
