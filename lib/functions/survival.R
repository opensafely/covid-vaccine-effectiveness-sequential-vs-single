library('tidyverse')
library('survival')
#library('flexsurv')

censor <- function(event_date, censor_date, na.censor=TRUE){
  # censors event_date to on or before censor_date
  # if na.censor = TRUE then returns NA if event_date>censor_date, otherwise returns min(event_date, censor_date)
  if (na.censor)
    dplyr::if_else(event_date>censor_date, as.Date(NA_character_), as.Date(event_date))
  else
    dplyr::if_else(event_date>censor_date, as.Date(censor_date), as.Date(event_date))
}

censor_indicator <- function(event_date, censor_date){
  # returns 0 if event_date is censored by censor_date, or if event_date is NA. Otherwise 1
  dplyr::if_else((event_date>censor_date) | is.na(event_date), FALSE, TRUE)
}

tte <- function(origin_date, event_date, censor_date, na.censor=FALSE){
  # returns time-to-event date or time to censor date, which is earlier

  if (na.censor)
    time <- dplyr::if_else(censor_date>=event_date, event_date-origin_date, as.Date(NA)-origin_date)
  else
    time <- pmin(event_date-origin_date, censor_date-origin_date, na.rm=TRUE)
  as.numeric(time)
}



round_tte <- function(time, width=7){
  # group follow-up time to be in periods of size `width`
  # eg, convert to weekly instead of dail with width=7
  # follow-up time of zero is always mapped to zero
  # then first period is mapped to `1`, second period is mapped to `2`, etc
  ceiling(time/width)
}

tidy_surv <-
  function(
    survfit,
    times = NULL,
    addtimezero=FALSE
  ) {

    # tidy post-fit survival dataset, with extra estimates than provided by broom::tidy.coxph

    mintime <- min(survfit$time)
    timezero <- min(0, mintime-1)


    if (is.null(times)) {
      output <-
        survfit %>%
        broom::tidy() %>%
        transmute(
          time,
          leadtime = lead(time),
          interval = leadtime - time,

          n.risk,
          n.event,
          n.censor,

          sumerand = n.event / ((n.risk - n.event) * n.risk),

          surv=cumprod(1 - n.event / n.risk),
          surv.ll = conf.low,
          surv.ul = conf.high,
          se.surv_greenwood = surv * sqrt(cumsum(sumerand)),

          # kaplan meier hazard estimates
          haz_km = n.event / (n.risk * interval), # =-(surv-lag(surv))/lag(surv)
          cml.haz_km = cumsum(haz_km), # =cumsum(haz_km)
          se.haz_km = haz_km * sqrt((n.risk - n.event) / (n.risk * n.event)),

          # actuarial hazard estimates
          haz_ac = n.event / ((n.risk - (n.censor / 2) - (n.event / 2)) * interval), # =(cml.haz-lag(cml.haz))/interval
          cml.haz_ac = -log(surv), #=cumsum(haz_ac)
          se.haz_ac = (haz_ac * sqrt(1 - (haz_ac * interval / 2)^2)) / sqrt(n.event),

          # log(-log()) scale

          llsurv = log(-log(surv)),
          se.llsurv = sqrt((1 / log(surv)^2) * cumsum(sumerand)),

      )
    }

    else {

      output <-
        survfit %>%
        broom::tidy() %>%
        complete(
          time = times,
          fill = list(n.event = 0, n.censor = 0)
        ) %>%
        fill(n.risk, .direction = c("up")) %>%
        transmute(
          time,
          leadtime = lead(time),
          interval = leadtime - time,

          n.risk,
          n.event,
          n.censor,

          sumerand = n.event / ((n.risk - n.event) * n.risk),

          surv=cumprod(1 - n.event / n.risk),
          surv.ll = conf.low,
          surv.ul = conf.high,

          se.surv_greenwood = surv * sqrt(cumsum(sumerand)),

          # kaplan meier hazard estimates
          haz_km = n.event / (n.risk * interval), # =-(surv-lag(surv))/lag(surv)
          cml.haz_km = cumsum(haz_km), # =cumsum(haz_km)
          se.haz_km = haz_km * sqrt((n.risk - n.event) / (n.risk * n.event)),

          # actuarial hazard estimates
          haz_ac = n.event / ((n.risk - (n.censor / 2) - (n.event / 2)) * interval), # =(cml.haz-lag(cml.haz))/interval
          cml.haz_ac = -log(surv), #=cumsum(haz_ac)
          se.haz_ac = (haz_ac * sqrt(1 - (haz_ac * interval / 2)^2)) / sqrt(n.event),

          # log(-log()) scale

          llsurv = log(-log(surv)),
          se.llsurv = sqrt((1 / log(surv)^2) * cumsum(sumerand)),
        )
    }

    if(addtimezero){
      output <- output %>%
        add_row(
          time = timezero,
          leadtime = mintime,
          interval = leadtime-time,
          sumerand=0,

          #estimate=1, std.error=0, conf.high=1, conf.low=1,

          surv=1,
          surv.ll=1,
          surv.ul=1,
          se.surv_greenwood=0,

          haz_km=0, se.haz_km=0, cml.haz_km=0,
          haz_ac=0, se.haz_ac=0, cml.haz_ac=0,
          .before=1
        )
    }

    return(output)
  }