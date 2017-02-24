function [estimates, model] = fitcurvedemo(xdata, ydata)
% Call fminsearch with a random starting point.

start_point = rand(1, 2b);
model = @expfun;
estimates = fminsearch(model, start_point);

% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for A * exp(-lambda * xdata) - ydata, 
% and the FittedCurve. FMINSEARCH only needs sse, but we want to 
% plot the FittedCurve at the end.

    function [sse, FittedCurve] = expfun(params)
        A = params(1);
	B = params(2);
        FittedCurve = A + 1./(xdata).^B;
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end
end