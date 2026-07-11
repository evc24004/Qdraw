rng(7);

N = 4;
d = 2^N;
shots = 4096;

A = diag(sqrt(1:d-1), 1);
Nop = diag(0:d-1);
Dop = @(al) expm(al*A' - conj(al)*A);
Rop = @(th) expm(1i*th*Nop);
Sop = @(z) expm((conj(z)*A*A - z*(A')*(A'))/2);
qp2al = @(q0,p0) (q0 + 1i*p0)/sqrt(2);

Ubody = Dop(qp2al(0, -0.9)) * Sop(-log(1.12));
Ustem = Dop(qp2al(0, 2.2)) * Sop(log(1.6));
Uleaf = Dop(qp2al(1.55, 2.0)) * Rop(-55*pi/180) * Sop(log(1.6));

objs = {0, Ubody, 0.30; 1, Ubody, 0.40; 0, Ustem, 0.10; 0, Uleaf, 0.17};
gates = {unitaryGate(1:N, Ubody), unitaryGate(1:N, Ubody), ...
         unitaryGate(1:N, Ustem), unitaryGate(1:N, Uleaf)};

nSet = 3^N;
nPauli = 4^N;
sig = {eye(2), [0 1;1 0], [0 -1i;1i 0], [1 0;0 -1]};

res = 128;
lim = 3.4;
qv = linspace(-lim, lim, res);
pv = linspace(-lim, lim, res);
[Q, P] = meshgrid(qv, pv);
AL = (Q(:) + 1i*P(:)) / sqrt(2);
nn = 0:d-1;
Cmat = exp(-abs(AL).^2/2) .* (AL .^ nn) ./ sqrt(factorial(nn));

H = zeros(res*res, 1);
for o = 1:size(objs, 1)
    k = objs{o,1};
    est = zeros(nPauli, 1);
    cnt = zeros(nPauli, 1);
    kb = dec2bin(k, N) - '0';
    prep = arrayfun(@(q) xGate(q), find(kb), 'UniformOutput', false);
    for setIdx = 0:nSet-1
        setting = zeros(1, N);
        v = setIdx;
        for q = 1:N, setting(q) = mod(v,3) + 1; v = floor(v/3); end
        basis = {};
        for q = 1:N
            if setting(q) == 1
                basis{end+1} = hGate(q); %#ok<*SAGROW>
            elseif setting(q) == 2
                basis{end+1} = rzGate(q, -pi/2);
                basis{end+1} = hGate(q);
            end
        end
        gl = [prep(:); gates(o); basis(:)];
        c = quantumCircuit(vertcat(gl{:}), N);
        m = randsample(simulate(c), shots);
        bits = char(m.MeasuredStates) - '0';
        w = m.Counts / shots;
        for mask = 1:2^N-1
            mb = bitget(mask, N:-1:1);
            val = ((-1).^(mod(bits*mb', 2)))' * w;
            pIdx = polyval(setting .* mb, 4) + 1;
            est(pIdx) = est(pIdx) + val;
            cnt(pIdx) = cnt(pIdx) + 1;
        end
    end
    rho = eye(d) / d;
    for pIdx = 2:nPauli
        if cnt(pIdx) == 0, continue; end
        v = pIdx - 1;
        dig = zeros(1, N);
        for q = N:-1:1, dig(q) = mod(v,4); v = floor(v/4); end
        Pm = 1;
        for q = 1:N, Pm = kron(Pm, sig{dig(q)+1}); end
        rho = rho + (est(pIdx)/cnt(pIdx)) * Pm / d;
    end
    H = H + objs{o,3} * max(real(sum((conj(Cmat) * rho) .* Cmat, 2)), 0);
    fprintf('component %d/%d done\n', o, size(objs,1));
end

H = reshape(H, res, res);
img = flipud(H / max(H(:)));
anchors = [0.0015 0.0005 0.0139; 0.0873 0.0444 0.2244; ...
           0.2588 0.0386 0.4064; 0.4166 0.0903 0.4328; ...
           0.5783 0.1481 0.4044; 0.7355 0.2154 0.3297; ...
           0.8654 0.3168 0.2262; 0.9541 0.4433 0.1201; ...
           0.9877 0.6521 0.2113; 0.9647 0.8434 0.4321; ...
           0.9884 0.9984 0.6449];
cmap = interp1(linspace(0,1,size(anchors,1)), anchors, linspace(0,1,256));
rgb = ind2rgb(1 + round(img*255), cmap);
[xg, yg] = meshgrid(linspace(1, res, 512));
up = zeros(512, 512, 3);
for ch = 1:3, up(:,:,ch) = interp2(rgb(:,:,ch), xg, yg, 'linear'); end
if ~exist('output', 'dir'), mkdir('output'); end
imwrite(up, fullfile('output', 'apple_qcg_hw.png'));
fprintf('wrote output/apple_qcg_hw.png\n');
