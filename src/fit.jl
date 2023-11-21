# using AffineInvariantMCMC

function find_RT(sp::SeriesPoint; nwalkers=100)
    numdims = 5
    thinning = 10
    numsamples_perwalker = 1000
    burnin = 100

    stds = exp(5 * randn(numdims))
    means = 1 + 5 * rand(numdims)
    llhood = x->begin
        retval = 0.
        for i = eachindex(x)
            retval -= .5 * ((x[i] - means[i]) / stds[i]) ^ 2
        end
        return retval
    end
    x0 = rand(numdims, numwalkers) * 10 - 5
    chain, llhoodvals = AffineInvariantMCMC.sample(llhood, numwalkers, x0, burnin, 1)
    chain, llhoodvals = AffineInvariantMCMC.sample(llhood, numwalkers, chain[:, :, end], numsamples_perwalker, thinning)
    flatchain, flatllhoodvals = AffineInvariantMCMC.flattenmcmcarray(chain, llhoodvals)
end
