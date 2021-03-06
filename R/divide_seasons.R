#' divide_seasons
#'
#' @param d data.frame of vegetation index
#' @param sp one line data.frame
#' @param sitename character
#' @param .movmean moving mean by phenofit::wSG
#'
#' @note site-year may be not continuous.
#'
#' @export
divide_seasons <- function(d, sp, nptperyear = 23,
    iters = 2,
    lambda = 100,
    wFUN = wTSM,
    r_max = 0.2, r_min = 0.0,
    r_minPeakHeight = 0.05,
    calendarYear = FALSE,
    .movmean = TRUE,
    .v_curve = FALSE,
    ...)
{
    sp = sp[1, , drop = FALSE]
    sitename = d$site[1]

    ## moving average first
    # check_season(sitename, df_raw, stations)
    # d     <- df[site == sitename & scale == "0m", .(t, date, y = EVI, w, SummaryQA)] #%T>% plotdata(365)
    # sp    <- st[site == sitename, ]

    south <- sp$lat < 0
    dnew  <- add_HeadTail(d, south = south, nptperyear)

    INPUT <- check_input(dnew$t, dnew$y, dnew$w, QC_flag = NULL, nptperyear,
        maxgap = ceiling(nptperyear/12*1.5),
        south = sp$lat < 0,
        date_start = d$t[1],
        date_end = last(d$t))

    frame = floor(INPUT$nptperyear/8) * 2 + 1

    if (.v_curve) {
        lg_lambdas <- seq(3.3, 5, 0.1) # 2000-
        r <- v_curve(INPUT, lg_lambdas, d = 2, IsPlot = FALSE)
        lambda = r$lambda
    }
    # if (.movmean) {
    #     r_wSG <- smooth_wSG(INPUT$y, INPUT$w, nptperyear = nptperyear, INPUT$ylu, iters = 2, frame = frame)
    #     INPUT$y = r_wSG$zs %>% last()
    # }
    # plot_input(INPUT)
    # browser()
    # parameters for season_mov
    threshold_max = 0.1
    nf = 1

    FUN_fit        <- "smooth_wWHIT"
    threshold_max  <- ifelse(cv_coef(d$y)[3] >= 1, 0.1, 0.2) # empirical param
    # FUN_fit <- ifelse(sp$IGBP %in% IGBP_forest, "wHANTS", "wWHIT")
    # "wBisquare"
    maxExtendMonth <- ifelse(sp$IGBP == "EBF", 2, 2)

    # wFUN <- "wBisquare", "wTSM", threshold_max = 0.1, IGBP = CSH
    # INPUT <- get_input(df, st, sitename)
    brks2  <- season_mov(INPUT,
        rFUN = get(FUN_fit),
        wFUN = wFUN,
        iters = iters, wmin = 0.1,
        IsOptim_lambda = FALSE,
        lambda = lambda, nf = 3, frame = frame,
        maxExtendMonth = 12,
        r_max = r_max, r_min = r_min,
        r_minPeakHeight = r_minPeakHeight,
        calendarYear = calendarYear,
        # ...,
        # IsPlot.vc = FALSE,
        # plotdat = INPUT, print = TRUE,
        # titlestr = "")
        IsPlot = FALSE, IsPlot.OnlyBad = FALSE,
        minpeakdistance = nptperyear/36*2, # 20 days
        MaxPeaksPerYear = 3,
        MaxTroughsPerYear = 4,
        ypeak_min = 0.08,
        ...
    )
    listk(INPUT, brks = brks2, lambda = lambda)
}
