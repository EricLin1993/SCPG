function y_thresholded = delta(y, threshold)
    s = abs(y) - threshold;
    s = (s + abs(s))/2;
    y_thresholded = sign(y).*s;
end

