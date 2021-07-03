{
    l <- with(params, season_mov(INPUT,
               rFUN = smooth_wWHIT,
               wFUN = wTSM,
               iters = 2, wmin = 0.1,
               IsOptim_lambda = FALSE,
               lambda = lambda,
               # nf = nf, frame = frame,
               maxExtendMonth = maxExtendMonth,
               # r_max = r_max, r_min = r_min,
               # r_minPeakHeight = r_minPeakHeight,
               # calendarYear = calendarYear,
               # ...,
               # IsPlot.vc = FALSE,
               # plotdat = INPUT, print = TRUE,
               # titlestr = "")
               IsPlot = TRUE,
               IsPlot.OnlyBad = FALSE,
               # minpeakdistance = nptperyear/36*2, # 20 days
               MaxPeaksPerYear = 3,
               MaxTroughsPerYear = 4))
}

## 1. good:
{
    l <- divide_seasons(dat_gpp, 365, is.plot = TRUE, lambda = 200,
                        maxExtendMonth = 12,
                        .v_curve = TRUE
    )
}

## 2. error: process_phenofit version

