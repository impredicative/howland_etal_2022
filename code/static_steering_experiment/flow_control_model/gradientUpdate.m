function [ x, m, v ] = gradientUpdate( optimizer, learning_rate, x, f, m, v, maxOrMin );
% Vanilla and ADAM gradient-based optimization update
beta1 = 0.9; beta2 = 0.999; eps = 10^(-8);

if strcmp(maxOrMin, 'min')==1;
    if strcmp(optimizer, 'vanilla') == true;
        % Vanilla
        x = x - learning_rate .* f;
    elseif strcmp(optimizer, 'adam') == true;
        % Adam
        m = beta1*m + (1-beta1)*f;
        v = beta2*v + (1-beta2)*(f.^2);
        x = x - learning_rate .* m(1:end) ./ (sqrt(v(1:end)) + eps);
    end
elseif strcmp(maxOrMin, 'max')==1;
    if strcmp(optimizer, 'vanilla') == true;
        % Vanilla
        x = x + learning_rate .* f;
    elseif strcmp(optimizer, 'adam') == true;
        % Adam
        m = beta1*m + (1-beta1)*f;
        v = beta2*v + (1-beta2)*(f.^2);
        x = x + learning_rate .* m(1:end) ./ (sqrt(v(1:end)) + eps);
    end
end

end

