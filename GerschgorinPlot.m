% to plot Gerschgorin graph, R must be square martix
% R = [1+1i, 2-2i; 1+2i, 3-1i];
% GerschgorinPlots(R)
function GerschgorinPlot(R)
    len = size(R, 1);
    stardand_circle = exp(1i * 2 * pi * (0: 1e-3: 1)');
    r = zeros(len, 1);
    for k = 1: 1: len - 1   % can be len\, don't be (len - 1)\
        r(k) = sum(abs(R(k, :))) - abs(R(k, k));
        cir = stardand_circle * r(k);
        Rcir = real(cir) + real(R(k, k));
        Icir = imag(cir) + imag(R(k, k));
        plot(Rcir, Icir); hold on;
        plot(real(R(k, k)), imag(R(k, k)), 'Marker', '*', 'MarkerSize', 10);
        hold on;
    end
    xlabel('Real axis'); ylabel('Imag axis');
    title('Gerschogorin circle'); grid on
end