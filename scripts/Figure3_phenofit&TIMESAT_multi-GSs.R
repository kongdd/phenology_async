# source("test/main_pkgs.R")
library(lubridate)
library(phenofit)
library(ggplot2)
library(data.table)
library(Ipaper)   # rpkgs/Ipaper
library(rTIMESAT) # eco-hydro/rTIMESAT
library(patchwork)

devtools::load_all("../phenofit.R")
devtools::load_all("../rTIMESAT.R")
devtools::load_all()

simu_VI <- function(SOS = 50, EOS = 100, rate = 0.1, mx = 0.6, mn = 0.1, year = 2010, wmin = 0.2) {
    par <- c(mn, mx, SOS, rate, EOS, rate)

    t <- seq(1, 365, 8)
    w <- rep(1, length(t))

    noise <- rnorm(n = length(t), mean = 0, sd = 0.05)
    I_noise <- noise < 0
    noise[!I_noise] <- 0
    w[I_noise] <- wmin
    y0 <- doubleLog.Beck(par, t)
    y  <- y0 + noise
    data.table(year, doy = t, t = as.Date(sprintf("%d%03d", year, t), "%Y%j"),
               y, y0, w)
}

# two growing season
{
    # d2_1 <- simu_VI(50, 120, 0.05, mn = 0.1, year = 2012)
    # d2_2 <- simu_VI(180, 250, 0.1, mn = 0.2, year = 2012)
    # d2_a = rbind(d2_1[doy < 150, ], d2_2[doy >= 150, ])
    d2_a <- simu_VI(100, 250, 0.1, year = 2010)
    d2_a[doy >= 250, y := y + 0.2]

    t = d2_a$doy
    tout = 1:366
    r <- with(d2_a, FitDL.AG2(y, doy, tout))
    d =  c(list(t = r$tout), r$zs) %>% as.data.table()

    ggplot(d2_a, aes(doy, y)) +
        geom_line(aes(y = y0), color = "black") +
        geom_line(aes(y = y), color = "green") +
        geom_line(data = d, aes(t, iter2), color = "red")
}

{
    set.seed(0)
    d1_a <- simu_VI(150, 250, 0.1, year = 2010)
    d1_b <- simu_VI(150, 250, 0.15, year = 2011)

    # two growing season
    d2_1 <- simu_VI(50, 120, 0.05, year = 2012)
    d2_2 <- simu_VI(180, 250, 0.1, year = 2012)
    d2_a = rbind(d2_1[doy < 150, ], d2_2[doy >= 150, ])

    d2_1 <- simu_VI(50, 120, 0.1, year = 2013)
    d2_2 <- simu_VI(180, 250, 0.05, year = 2013)
    d2_b = rbind(d2_1[doy < 150, ], d2_2[doy >= 150, ])

    # triple growing season
    d3_1 <- simu_VI(25, 75, 0.05, year = 2014)
    d3_2 <- simu_VI(100, 150, 0.1, year = 2014)
    d3_3 <- simu_VI(200, 300, 0.1, year = 2014)
    d3_a = rbind(d3_1[doy < 85, ],
                 d3_2[doy %in% c(85:175), ],
                 d3_3[doy > 175, ])

    d3_1 <- simu_VI(25, 75, 0.05, year = 2015)
    d3_2 <- simu_VI(100, 150, 0.1, year = 2015)
    d3_3 <- simu_VI(200, 300, 0.1, year = 2015)
    d3_b = rbind(d3_1[doy < 85, ],
                 d3_2[doy %in% c(85:175), ],
                 d3_3[doy > 175, ])

    dat = rbind(d1_a, d1_b, d2_a, d2_b, d3_a, d3_b)
    # dat$w %<>% as.factor()
    nptperyear = 46
    ggplot(dat, aes(t, y)) +
        geom_line(aes(y = y0), color = "black") +
        geom_line(aes(y = y), color = "green")
    # geom_point(aes(color = w, shape = w))
}

{
    devtools::load_all("../phenofit.R/")
    r_pheno = process_phenofit(dat$y, dat$t, dat$w,
                        nptperyear = nptperyear,
                        nextend = 0,
                        yticks = c(0, 0.3, 0.6),
                        # wFUN = wBisquare_julia,
                        wFUN = wTSM,
                        # outfile = "Figure3_phenofit_multi-GSs.pdf",
                        overwrite = T, show = TRUE)
}

{
    meths = c("wWHIT", "AG", "Zhang")
    d_obs = r_pheno$data
    d_season = r_pheno$brks$dt
    d_fit <- rbind(
        r_pheno$brks$fit %>% select(t, y, starts_with("ziter")) %>% mutate(meth = "wWHIT"),
        r_pheno$dfit %>% select(t, y, starts_with("ziter"), meth)
    ) %>% .[meth %in% meths] %>%
        dplyr::rename(z = ziter2)
    d_fit$meth %<>% factor(meths)

    mytheme = theme(
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 14, family = "Times"),
        axis.title = element_text(size = 14),
        # legend.position = "bottom",
        legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size = 12),
        legend.margin =  margin(1, 1, 1, 1),
        plot.margin = margin(0, 1, -4, 1)*2)

    p_kong = plot_phenofit2(d_obs, d_season, d_fit) + mytheme +
        labs(x = "Time", y = "EVI",
             color = "Fitting:",
             title = expression("(a)"~italic(phenofit))) +
        theme(legend.key.size = unit(0.6, 'cm'))

    # devtools::load_all("../rTIMESAT.R/")
    r = TIMESAT_process(dat, nptperyear, half_win = 4, p_trs = 0.02,
                        # methods = "SG",
                        seasonpar = 0.0,
                        cache = FALSE)
    p_TS = plot_TIMESAT(dat, r) + mytheme +
        labs(x = "Time", y = "EVI", color = "Fitting:",
             title = expression("(b) TIMESAT"))
    # write_fig(p_TS, "Figure3b_TIMESAT_multi-GSs.pdf", 10, 3)
    write_fig( p_kong / p_TS, "Figure3_phenofit&TIMESAT_multi-GSs_v2.pdf", 10, 6.5)
}
