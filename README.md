
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Ridge Regression for Paired Comparisons

This document illustrates how to use the R (R Core Team, 2024) code
accompanying [Varin and Firth (2024)](https://arxiv.org/abs/2406.09597)
to replicate the analysis of the last 28 seasons of English Premier
League. The figures in the paper differ in minor respects from those
produced in this document due to some manual editing for inclusion in
the paper.

<!-- badges: start -->
<!-- badges: end -->

Load some libraries that will be used in this document:

``` r
library("dplyr")
library("tidyr")
library("ggplot2")
library("signs")
```

Load the code that implements the methodology discussed in [Varin and
Firth (2024)](https://arxiv.org/abs/2406.09597):

``` r
source("peb.R")
```

Ridge regression for paired comparisons with the pairwise empirical
Bayes method is implemented in the`peb` function:

``` r
peb <- function(tournament, lambda = NULL, home_effect = FALSE, a = 2 * nlevels(tournament$home), eps = 1e-4, ...)
```

The main argument of `peb` is `tournament`. The code assumes that
`tournament` is a `data.frame` with at least three columns called
`home`, `away` and `outcome`. Columns `home` and `away` are factors that
identify the teams playing at home and away, respectively. Column
`outcome` is the outcome of the match coded as `1` in case of victory of
the away team, `2` in case of draw and `3` in case of victory of the
home team. Further arguments of `peb` are:

- the ridge tuning parameter `lambda`. The default value is `NULL` in
  which case $\lambda$ is estimated by the pairwise empirical Bayes
  method;
- the logical variable `home_effect` that indicates whether a home-field
  effect parameter has to be included (`TRUE`) or not (`FALSE`);
- the small-sample adjustment parameter `a`. The default value for `a`
  is `2 * levels(tournament$data)` which corresponds in adding one
  imaginary concordant pair and one imaginary discordant pair to the
  counts for each of the teams, as discussed in the paper. Estimation
  without the small-sample adjustment is obtained with `a=0`.

Please refer to [Varin and Firth
(2024)](https://arxiv.org/abs/2406.09597) for further details on the
meaning of the quantities indicated above.

Load English Premier League results for the seasons 1995-1996 to
2022-2023:

``` r
load("premier_1995-2023.RData")
```

The loaded object is a `data.frame` with the results of all 10640
matches in the 28 seasons from 1995-1996 to 2022-2023:

``` r
class(premier_data)
#> [1] "data.table" "data.frame"
```

``` r
head(premier_data)
#>         date week          home            away outcome    season
#> 1 1995-08-19    1   Southampton Nott'ham Forest       1 1995-1996
#> 2 1995-08-19    1 Newcastle Utd   Coventry City       3 1995-1996
#> 3 1995-08-19    1     Wimbledon          Bolton       3 1995-1996
#> 4 1995-08-19    1     Liverpool  Sheffield Weds       3 1995-1996
#> 5 1995-08-19    1      West Ham    Leeds United       1 1995-1996
#> 6 1995-08-19    1     Blackburn             QPR       3 1995-1996
```

The long-term proportions of home wins, draws and away wins are

``` r
options(digits = 2)
ref_probs <- prop.table(table(premier_data$outcome))
ref_probs
#> 
#>    1    2    3 
#> 0.29 0.25 0.46
```

For illustration, let us fit the paired comparison model to season
2015/2016 without home-field effect:

``` r
season_1516 <- premier_data %>%
  filter(season == "2015-2016") 
fit <- peb(season_1516, home_effect = FALSE)
```

The fitted object is a `list` with various components:

``` r
class(fit)
#> [1] "list"
```

``` r
names(fit)
#>  [1] "ties"        "home_effect" "delta"       "gamma"       "tau"        
#>  [6] "lambda"      "X"           "lik"         "score"       "mu"
```

Among the components of `fit` are the model parameters $\gamma$ and
$\delta$:

``` r
fit$gamma
#> [1] 0.37
```

``` r
fit$delta
#> [1] 0
```

The estimate of $\delta$ is set equal to zero since no home-field effect
is assumed. The tuning parameter $\lambda$ estimated with the pairwise
empirical Bayes method is

``` r
fit$lambda
#> [1] 8.1
```

The most relevant component of `fit` is the vector $\mu$ of the
estimated strength parameters, below ordered from best to worst rated
team:

``` r
sort(fit$mu, decreasing = TRUE)
#>  Leicester City         Arsenal       Tottenham  Manchester Utd Manchester City 
#>           0.562           0.368           0.365           0.256           0.252 
#>        West Ham     Southampton       Liverpool         Chelsea       Blackburn 
#>           0.209           0.191           0.168           0.011           0.000 
#>          Bolton   Coventry City    Leeds United   Middlesbrough Nott'ham Forest 
#>           0.000           0.000           0.000           0.000           0.000 
#>             QPR  Sheffield Weds       Wimbledon    Derby County        Barnsley 
#>           0.000           0.000           0.000           0.000           0.000 
#>    Charlton Ath   Bradford City    Ipswich Town          Fulham Birmingham City 
#>           0.000           0.000           0.000           0.000           0.000 
#>      Portsmouth          Wolves  Wigan Athletic         Reading   Sheffield Utd 
#>           0.000           0.000           0.000           0.000           0.000 
#>       Hull City         Burnley       Blackpool    Cardiff City        Brighton 
#>           0.000           0.000           0.000           0.000           0.000 
#>    Huddersfield       Brentford      Stoke City         Everton    Swansea City 
#>           0.000           0.000          -0.022          -0.052          -0.083 
#>       West Brom         Watford     Bournemouth  Crystal Palace      Sunderland 
#>          -0.129          -0.145          -0.193          -0.196          -0.225 
#>   Newcastle Utd    Norwich City     Aston Villa 
#>          -0.278          -0.361          -0.697
```

Now, we refit the model with the inclusion of the home-field effect:

``` r
fit2 <- peb(season_1516, home_effect = TRUE)
```

The estimates of the model parameters and the tuning parameter become

``` r
fit2$gamma
#> [1] 0.37
```

``` r
fit2$delta
#> [1] 0.14
```

``` r
fit2$lambda
#> [1] 9.1
```

The corresponding vector of estimated strength parameters is

``` r
sort(fit2$mu, decreasing = TRUE)
#>  Leicester City         Arsenal       Tottenham Manchester City  Manchester Utd 
#>          0.5491          0.3629          0.3562          0.2530          0.2502 
#>        West Ham     Southampton       Liverpool         Chelsea       Blackburn 
#>          0.2046          0.1890          0.1594          0.0072          0.0000 
#>          Bolton   Coventry City    Leeds United   Middlesbrough Nott'ham Forest 
#>          0.0000          0.0000          0.0000          0.0000          0.0000 
#>             QPR  Sheffield Weds       Wimbledon    Derby County        Barnsley 
#>          0.0000          0.0000          0.0000          0.0000          0.0000 
#>    Charlton Ath   Bradford City    Ipswich Town          Fulham Birmingham City 
#>          0.0000          0.0000          0.0000          0.0000          0.0000 
#>      Portsmouth          Wolves  Wigan Athletic         Reading   Sheffield Utd 
#>          0.0000          0.0000          0.0000          0.0000          0.0000 
#>       Hull City         Burnley       Blackpool    Cardiff City        Brighton 
#>          0.0000          0.0000          0.0000          0.0000          0.0000 
#>    Huddersfield       Brentford      Stoke City         Everton    Swansea City 
#>          0.0000          0.0000         -0.0220         -0.0461         -0.0813 
#>       West Brom         Watford  Crystal Palace     Bournemouth      Sunderland 
#>         -0.1248         -0.1449         -0.1865         -0.1903         -0.2205 
#>   Newcastle Utd    Norwich City     Aston Villa 
#>         -0.2770         -0.3558         -0.6824
```

Now we load the function that is used to calculate the comparison
methods which are maximum likelihood and bias-reduced maximum
likelihood. Estimates with these methods are calculated with the
`bpolr.R` function of Kosmidis (2014). This function can be downloaded
from Ioannis Kosmidis’ [web page](https://www.ikosmidis.com/):

``` r
library("archive")
url_address <- "https://www.ikosmidis.com/files/rssb12025-sup-0002-supplementary-materials-ikosmidis.zip"
archive_extract(url_address, dir = "tmp")
source("tmp/supplementary/bpolr.R")
```

The next lines replicate the forecasting exercises of [Varin and Firth
(2024)](https://arxiv.org/abs/2406.09597). For each season, paired
comparison models are trained with the first $k$ matchweeks (with
$k=10, 15, 20, 25, 30$) and then validated by predicting the rest of the
tournament:

``` r
train_weeks <- seq(10, 30, by = 5)
n_weeks <- length(train_weeks)
methods <- c("PEB", "MLE", "BRMLE")
n_methods <- length(methods)
seasons <- unique(premier_data$season)
n_seasons <- length(seasons)
premier_skill_scores <- data.frame(season = rep(seasons, each = n_methods * n_weeks), week = rep(rep(train_weeks, each = n_methods, n_seasons)), skill_score = NA, method = rep(rep(methods, n_seasons), n_weeks))
premier_skill_scores$method <- factor(premier_skill_scores$method, levels = c("MLE", "BRMLE", "PEB"))

for (s in seasons) {
  print(s)
  for (w in train_weeks) {
    this_season <- premier_data %>% filter(season == s)
    this_season <- droplevels(this_season)
    train_set <- this_season %>% filter(week <= w)
    test_set <- this_season %>% filter(week > w)
    peb_fit <- peb(train_set, home_effect = TRUE, a = 40)
    mle_fit <- brmle(train_set, br = FALSE, maxit = 5000)
    br_fit <- brmle(train_set, br = TRUE, maxit = 5000)
    rows <- which(premier_skill_scores$season == s & 
                  premier_skill_scores$week == w)
    premier_skill_scores[rows, "skill_score"] <- 
      c(compute_skillscore(peb_fit, test_set, "probit", ref_probs),
        compute_skillscore(mle_fit, test_set, "logit", ref_probs),
        compute_skillscore(br_fit, test_set, "logit", ref_probs))
  }
}
#> [1] "1995-1996"
#> [1] "1996-1997"
#> [1] "1997-1998"
#> [1] "1998-1999"
#> [1] "1999-2000"
#> [1] "2000-2001"
#> [1] "2001-2002"
#> [1] "2002-2003"
#> [1] "2003-2004"
#> [1] "2004-2005"
#> [1] "2005-2006"
#> [1] "2006-2007"
#> [1] "2007-2008"
#> [1] "2008-2009"
#> [1] "2009-2010"
#> [1] "2010-2011"
#> [1] "2011-2012"
#> [1] "2012-2013"
#> [1] "2013-2014"
#> [1] "2014-2015"
#> [1] "2015-2016"
#> [1] "2016-2017"
#> [1] "2017-2018"
#> [1] "2018-2019"
#> [1] "2019-2020"
#> [1] "2020-2021"
#> [1] "2021-2022"
#> [1] "2022-2023"
```

The results of the forecasting exercises are summarized in Figure 2 of
[Varin and Firth (2024)](https://arxiv.org/abs/2406.09597) which reports
the boxplot summary of the distributions of the logarithmic skill scores
calculated after 10, 15, 20, 25 and 30 weeks of matches:

``` r
## select a set of colors that are safe for colorblind people:
my_colors <- c("#E69F00", "#D55E00", "#0072B2")
ref_color <- "#999999"

my_theme <- theme_minimal() +
  theme(axis.title.x = element_text(),
        axis.title.y = element_text(angle = 0))
    
boxplot_premier <- function() {
  week_names <- as_labeller(
      c(`10` = "Matchweek 10", 
        `15` = "Matchweek 15",
        `20` = "Matchweek 20", 
        `25` = "Matchweek 25",
        `30` = "Matchweek 30"))

  my_plot <- premier_skill_scores %>% 
    ggplot(aes(y = skill_score, x = method, fill = method)) +
    geom_boxplot(alpha = 0.8) +
    facet_wrap(~week, ncol = 5, labeller = week_names) + 
    geom_hline(aes(yintercept = 0.0), col = ref_color, alpha = 0.8) +
    scale_y_continuous(labels = signs_format(accuracy = .01)) +
    scale_fill_manual(values = my_colors) +
    labs(x = NULL, y = NULL, fill = NULL) + 
    my_theme +
    theme(axis.text.x = element_blank(), 
          legend.position = "bottom")
  ## 
  my_plot
}
boxplot_premier() 
```

<img src="README_files/figure-gfm/unnamed-chunk-18-1.png" width="100%" />

We now focus on the most recent Premier League season included in our
study and replicate Figure 3 of [Varin and Firth
(2024)](https://arxiv.org/abs/2406.09597):

``` r
skills_2223 <- premier_skill_scores %>%
  filter(season == "2022-2023")
    
skills_2223  %>% 
    ggplot(aes(y = skill_score, x = week, colour = method)) + 
    geom_point(aes(shape = method), alpha = 0.8, size = 3) +
    geom_line(alpha = 0.8, size = .8) + 
    labs(x = "Matchweeks used for training", y = NULL, title = NULL, colour = NULL, shape = NULL) +
    scale_color_manual(values = my_colors, guide = guide_legend(nrow = 1)) +
    scale_y_continuous(labels = signs_format(accurcy = 0.01)) +
    geom_hline(aes(yintercept = 0.0), col = ref_color, alpha = 0.5) +
    my_theme +
    theme(legend.position = "bottom", legend.text = element_text(size = 13))
```

<img src="README_files/figure-gfm/unnamed-chunk-19-1.png" width="100%" />

We conclude with Figure 4 from [Varin and Firth
(2024)](https://arxiv.org/abs/2406.09597). Drawing this figure is more
complicated than replicating previous figures due to the various
annotations included:

``` r
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## First step: extract the data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

season2223 <- premier_data %>% filter(season == "2022-2023")
season2223 <- droplevels(season2223)
## after 10 weeks:
train_w10 <- season2223 %>% filter(week <= 10)
peb_w10 <- peb(train_w10, home_effect = TRUE, a = 40)
mle_w10 <- brmle(train_w10, br = FALSE)
br_w10 <- brmle(train_w10, br = TRUE)
## at the end of the season:
peb_final <- peb(season2223, home_effect = TRUE, a = 40)
mle_final <- brmle(season2223, br = FALSE)
br_final <- brmle(season2223, br = TRUE)
## save results in a data.frame:
premier_2223 <- data.frame(week = c(rep(10, 20), rep(38, 20)),
                           team = factor(rep(names(peb_w10$mu), 2)),
                           PEB = c(peb_w10$mu, peb_final$mu),
                           MLE = c(mle_w10$mu, mle_final$mu),
                           BRMLE = c(br_w10$mu, br_final$mu))
premier_2223 <- pivot_longer(premier_2223, cols=3:5, names_to = "method", values_to = "log_score")
premier_2223$method <- factor(premier_2223$method, 
                              levels = c("MLE", "BRMLE", "PEB"))

## highlight Arsenal, Leicester City and Manchester City:
to_highlight <- which(premier_2223$team %in% c("Arsenal", "Leicester City", "Manchester City"))
premier_2223$highlight <- "no"
premier_2223$highlight[to_highlight] <- "yes"
premier_2223$highlight <- as.factor(premier_2223$highlight)
premier_2223$color <- as.character(premier_2223$method)
premier_2223$color[premier_2223$highlight == "no"] <- "no"
premier_2223$color <- factor(premier_2223$color, levels = c(levels(premier_2223$method), "no"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Second step: create the plot
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p <- premier_2223 %>%
  ggplot(aes(x = week, y = log_score, colour = color, group = team, alpha = color)) +
  geom_line(size = .8) + 
  geom_point(size = 1) + 
  scale_y_continuous(labels = "0", breaks = 0) +
  coord_cartesian(xlim=c(2, 47)) +
  my_theme + theme(axis.text.y  = element_text(colour = gray(.6)), 
                   axis.text.x  = element_blank(), 
                   axis.ticks.x = element_blank(),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.minor.y = element_blank()) +
  labs(y = NULL, x = NULL)+
  facet_wrap(~method) +
  theme(legend.position = "none") +
  scale_color_manual(values = c(my_colors, gray(.6)))  + 
  scale_alpha_manual(values = c(1, 1, 1, 0.3))

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Last step: adding annotations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## extract MLE, BRMLE and PEB estimates after 10 weeks:
initial_mle <- premier_2223$log_score[premier_2223$team %in% c("Arsenal", "Manchester City", "Leicester City") & premier_2223$method == "MLE" & premier_2223$week == 10]
initial_brmle <- premier_2223$log_score[premier_2223$team %in% c("Arsenal", "Manchester City", "Leicester City") & premier_2223$method == "BRMLE" & premier_2223$week == 10]
initial_peb <- premier_2223$log_score[premier_2223$team %in% c("Arsenal", "Manchester City", "Leicester City") & premier_2223$method == "PEB" & premier_2223$week == 10]

## extract MLE, BRMLE and PEB estimates at the end of the season:
final_mle <- premier_2223$log_score[premier_2223$team %in% c("Arsenal", "Manchester City", "Leicester City") & premier_2223$method == "MLE" & premier_2223$week == 38]
final_brmle <- premier_2223$log_score[premier_2223$team %in% c("Arsenal", "Manchester City", "Leicester City") & premier_2223$method == "BRMLE" & premier_2223$week == 38]
final_peb <- premier_2223$log_score[premier_2223$team %in% c("Arsenal", "Manchester City", "Leicester City") & premier_2223$method == "PEB" & premier_2223$week == 38]

## create a data.frame with the annotation info (some coordinates were chosen by trial and error to make the plot readable):
initial_x <- c(5, 5, 4.5) 
final_x <- c(43, 43, 43.5) 
w <- c(21, 15.5, 18, rep(initial_x, 3), rep(final_x, 3))
## add a small adjustment for PEB coordinates to avoid overlapping:
ls <- c(2.85, 2.1, -2.7, initial_mle, initial_brmle, initial_peb + c(.06, -.06, 0), final_mle, final_brmle, final_peb + c(-.02, .02, 0))
mt <- factor(c(rep("MLE", 3), rep(c(rep("MLE", 3), rep("BRMLE", 3), rep("PEB", 3)), 2)),
             levels = levels(premier_2223$method))

ann_text <- data.frame(week = w, log_score = ls, method = mt, color = mt, team = "Arsenal")
labels_ls <- c(initial_mle, initial_brmle, initial_peb,
               final_mle, final_brmle, final_peb)
label <- c("ARS", "MNC", "LEI", signs(labels_ls, accuracy = .01))
## this is the end
p + geom_text(data = ann_text, label = label)
```

<img src="README_files/figure-gfm/unnamed-chunk-20-1.png" width="100%" />

## References

Kosmidis I. (2014). [Improved estimation in cumulative link
models](https://academic.oup.com/jrsssb/article/76/1/169/7075945).
*Journal of the Royal Statistical Society Series B: Statistical
Methodology* **76**(1), 169–196.

R Core Team (2024). *R: A Language and Environment for Statistical
Computing*. R Foundation for Statistical Computing, Vienna, Austria.
<https://www.R-project.org/>.

Varin C. and Firth D. (2024). Ridge Regression for Paired Comparisons: A
Tractable New Approach, with Application to Premier League Football.
<https://arxiv.org/abs/2406.09597>
