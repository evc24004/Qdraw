glyphs = ["U" "C" "O" "N" "N2" "HUSKY"];

N = 5;
d = 2^N;
shots = 4096;

A = diag(sqrt(1:d-1), 1);
Nop = diag(0:d-1);
Dop = @(al) expm(al*A' - conj(al)*A);
Rop = @(th) expm(1i*th*Nop);
Sop = @(z) expm((conj(z)*A*A - z*(A')*(A'))/2);
qp2al = @(q0,p0) (q0 + 1i*p0)/sqrt(2);
Ieye = eye(d);

vbar = @(q0,p0,s) Dop(qp2al(q0,p0)) * Sop(log(s));
hbar = @(q0,p0,s) Dop(qp2al(q0,p0)) * Sop(-log(s));
dbar = @(q0,p0,s,deg) Dop(qp2al(q0,p0)) * Rop(deg*pi/180) * Sop(log(s));

spec = struct();
spec.U = { ...
    0, vbar(-1.4, 0.7, 2.1), 1; ...
    0, vbar( 1.4, 0.7, 2.1), 1; ...
    0, hbar( 0, -1.9, 2.0),  1};
spec.C = { ...
    0, hbar( 0.3,  2.1, 1.8), 1; ...
    0, vbar(-1.4,  0,   2.1), 1; ...
    0, hbar( 0.3, -2.1, 1.8), 1};
spec.O = { ...
    2, Ieye, 0.45; ...
    3, Ieye, 0.55};
spec.N = { ...
    0, vbar(-1.5, 0, 2.2), 1; ...
    0, vbar( 1.5, 0, 2.2), 1; ...
    0, dbar( 0,   0, 2.9, 44), 1};
spec.N2 = spec.N;
spec.HUSKY = { ...
    0, Dop(qp2al(0, 0.3)), 0.16; ...
    1, Dop(qp2al(0, 0.3)), 0.22; ...
    0, dbar(-1.5, 2.5, 1.6,  15), 0.11; ...
    0, dbar( 1.5, 2.5, 1.6, -15), 0.11; ...
    0, vbar( 0, -1.2, 1.3), 0.16; ...
    0, hbar( 0, -2.3, 1.3), 0.08; ...
    0, Dop(qp2al(-1.7, -0.6)), 0.08; ...
    0, Dop(qp2al( 1.7, -0.6)), 0.08};

nSet = 3^N;
nPauli = 4^N;
sig = {eye(2), [0 1;1 0], [0 -1i;1i 0], [1 0;0 -1]};
res = 128;
lim = 4.5;
qv = linspace(-lim, lim, res);
pv = linspace(-lim, lim, res);
[Q, P] = meshgrid(qv, pv);
AL = (Q(:) + 1i*P(:)) / sqrt(2);
nn = 0:d-1;
Cmat = exp(-abs(AL).^2/2) .* (AL .^ nn) ./ sqrt(factorial(nn));

for g = 1:numel(glyphs)
    name = char(glyphs(g));
    objs = spec.(name);
    H = zeros(res*res, 1);
    for o = 1:size(objs, 1)
        k = objs{o,1};
        gate = unitaryGate(1:N, objs{o,2});
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
            gl = [prep(:); {gate}; basis(:)];
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
        H = H + objs{o,3} * max(real(sum((conj(Cmat)*rho).*Cmat, 2)), 0);
        fprintf('%s: component %d/%d done\n', name, o, size(objs,1));
    end
    H = reshape(H, res, res);
    if ~exist('output', 'dir'), mkdir('output'); end
    save(fullfile('output', sprintf('glyph_%s.mat', name)), 'H', 'qv', 'pv');
    imwrite(flipud(H/max(H(:))), ...
        fullfile('output', sprintf('glyph_%s_preview.png', name)));
    fprintf('saved output/glyph_%s.mat\n', name);
end
