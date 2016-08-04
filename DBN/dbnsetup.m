function dbn = dbnsetup(dbn, x, opts)
n = size(x, 2);
dbn.sizes = [n, dbn.sizes];
if (strcmp(opts.visibleDist,'Gauss')==1)
    for u = 1 : numel(dbn.sizes) - 1
        if (u==1)
            dbn.rbm{u}.alpha    = opts.alpha;
            dbn.rbm{u}.momentum = opts.momentum;

            dbn.rbm{u}.W  = zeros(dbn.sizes(u + 1), dbn.sizes(u));
            dbn.rbm{u}.vW = zeros(dbn.sizes(u + 1), dbn.sizes(u));

            dbn.rbm{u}.a  = zeros(dbn.sizes(u), 1);
            dbn.rbm{u}.va = zeros(dbn.sizes(u), 1);

            dbn.rbm{u}.c  = zeros(dbn.sizes(u + 1), 1);
            dbn.rbm{u}.vc = zeros(dbn.sizes(u + 1), 1);
        else
            dbn.rbm{u}.alpha    = opts.alpha;
            dbn.rbm{u}.momentum = opts.momentum;

            dbn.rbm{u}.W  = zeros(dbn.sizes(u + 1), dbn.sizes(u));
            dbn.rbm{u}.vW = zeros(dbn.sizes(u + 1), dbn.sizes(u));

            dbn.rbm{u}.b  = zeros(dbn.sizes(u), 1);
            dbn.rbm{u}.vb = zeros(dbn.sizes(u), 1);

            dbn.rbm{u}.c  = zeros(dbn.sizes(u + 1), 1);
            dbn.rbm{u}.vc = zeros(dbn.sizes(u + 1), 1);
        end
    end
else
    for u = 1 : numel(dbn.sizes) - 1
        dbn.rbm{u}.alpha    = opts.alpha;
        dbn.rbm{u}.momentum = opts.momentum;
        
        dbn.rbm{u}.W  = zeros(dbn.sizes(u + 1), dbn.sizes(u));
        dbn.rbm{u}.vW = zeros(dbn.sizes(u + 1), dbn.sizes(u));
        
        dbn.rbm{u}.b  = zeros(dbn.sizes(u), 1);
        dbn.rbm{u}.vb = zeros(dbn.sizes(u), 1);
        
        dbn.rbm{u}.c  = zeros(dbn.sizes(u + 1), 1);
        dbn.rbm{u}.vc = zeros(dbn.sizes(u + 1), 1);
    end
end

end
