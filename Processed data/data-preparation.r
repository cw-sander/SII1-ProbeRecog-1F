## ---------------------------
## Script name: data_preparation.r
##
## Purpose of script:
##  Creates a long data.frame with rts and ratings.
##  For each combination of five cutoff-criteria for slow outliers
##  (none, M + 2SD, 2500 ms, 2000 ms, 1500 ms), and three transformation
##  criteria (none, log, inverse) a column with reaction times is created
##  with the outliers being replaced by NAs.
##
## Author: Carsten Sander
##
## Date Created: 2022-02-18
##
## Copyright (c) Carsten Sander, 2022
## Email: carsten.sander@uni-hamburg.de
## ---------------------------
## Notes:
##  Before running this script, make sure your working directory
##  is set to a folder containing the /Processed data/ folder from
##  the github repository
##
##  R Version -- 4.1.2
## ---------------------------

# Load packages
library(tidyverse)

# Download raw data from Gitlab
paths <- read.csv("Processed data/raw-data-paths.csv", header = FALSE)[, 1]
pav <- lapply(paths, read_csv)
qua <- read.csv("https://gitlab.pavlovia.org/csander/sii1-proberecog-1f/raw/master/additional%20data/qualtrics.csv")[-c(1, 2), ] # nolint
stims <- read.csv("https://gitlab.pavlovia.org/csander/sii1-proberecog-1f/raw/master/additional%20data/stimuli.csv") # nolint

# ---------------------------

# Process pavlovia data
trials <- data.frame()
ratings <- data.frame()
subjects <- data.frame()

# Loop over all individual data.frames
for (i in seq(length(pav))) {
    # Select relevant data for analysis
    trials_i <- pav[[i]] %>%
        filter(
            trial_type == "test",
            probe_type %in% c("implied", "implied_other")
        ) %>%
        select(
            subject_id = ID, item_id, label = probe, probe_type,
            is_correct = c_probe_resp.corr, rt = c_probe_resp.rt
        ) %>%
        mutate(
            probe_type = recode(probe_type,
                "implied" = "im", "implied_other" = "io"))
    # Calculate individual cutoff (M + 2 * SD) based on correct responses
    rts <- trials_i$rt[trials_i$is_correct == 1]
    m2sd <- mean(rts) + 2 * sd(rts)
    # Define matrix of cutoff/transformation combinations
    corrections <- data.frame(
        cutoff = rep(c("none", "M2SD", "f25", "f20", "f15"), 3),
        cutoff_value = rep(c(Inf, m2sd, 2.5, 2.0, 1.5), 3),
        trans = c(rep("none", 5), rep("log", 5), rep("inv", 5)),
        trans_function = c(rep("{}", 5), rep("log({})", 5), rep("(1 / {}) * -1", 5))) # nolint
    # For each combination of cutoff criteria and transformations
    for (c in seq(NROW(corrections))) {
        # Get values for current correction
        var_name <- paste0("rt_", corrections$cutoff[c], "_", corrections$trans[c]) # nolint
        cutoff_value <- corrections$cutoff_value[c]
        trans_function <- corrections$trans_function[c]
        trans_function <- gsub("\\{\\}", "rts", trans_function)
        # Apply corrections to rts
        rts <- trials_i$rt
        rts[rts > cutoff_value] <- NA
        rts <- eval(parse(text = trans_function))
        # Append data to data.frame
        trials_i[, var_name] <- rts
    }
    # Build data.frame from individual datasets
    trials <- dplyr::bind_rows(trials, trials_i)

    # Get summary data on subjects (mean rt and rate of correct
    # responses to exclude those not meeting the pre-defined criteria)
    subjects[i, "ID"] <- trials_i$subject_id[1]
    subjects[i, "correct_rate"] <- sum(trials_i$is_correct) / nrow(trials_i)
    subjects[i, "mean_rt"] <- mean(trials_i$rt_M2SD_none[trials_i$is_correct == 1], na.rm = TRUE) # nolint

    # Get ratings of labels regarding valence
    # identification, and availability
    ratings_i <- pav[[i]] %>%
        filter(rating_block != "") %>%
        select(
            subject_id = ID,
            label = rated_label,
            response = rating_response,
            block = rating_block
        ) %>%
        pivot_wider(
            names_from = block,
            values_from = response,
            names_glue = "rating_{block}"
        )
    ratings <- dplyr::bind_rows(ratings, ratings_i)
}

# ---------------------------

# Calculate cutoff to exclude participants whose mean
# reaction time was slower than two SDs over the sample mean
cutoff <- mean(subjects$mean_rt) + 2 * sd(subjects$mean_rt)
subjects <- subjects %>%
    # Merge summary data with qualtrics data
    left_join(qua, by = "ID") %>%
    filter(part == 2) %>%
    # Recode survey items
    mutate(
        consent_given = ifelse(
            grepl("Ja", informed_consent.) |
            grepl("Nein", informed_consent_dc.),
            1, 0),
        compliance = as.numeric(recode(compliance.,
            "überhaupt nicht\1" = "1", "sehr gewissenhaft\n10" = "10")),
        age = as.numeric(str_extract(age, "\\d*")),
        gender = recode(gender,
            "anderes (z.B. nicht-binär)" = "other",
            "männlich" = "male", "weiblich" = "female"),
        pol_orientation = as.numeric(recode(pol_orientation,
            "links" = "1", "rechts" = "10")),
        pol_interest = as.numeric(recode(pol_interest.,
            "sehr stark" = "10", "überhaupt nicht" = "1")),
        pol_justice_freedom_1 = recode(pol_justice_freedom_1,
            "voll übereinstimmen" = 1, "weitgehend übereinstimmen" = 2,
            "weitgehend ablehnen" = 3, "voll und ganz ablehnen" = 4),
        pol_justice_freedom_2 = recode(pol_justice_freedom_2,
            "voll übereinstimmen" = 1, "weitgehend übereinstimmen" = 2,
            "weitgehend ablehnen" = 3, "voll und ganz ablehnen" = 4),
        pol_equal_treatment_1 = recode(pol_equal_treatment_1,
            "voll übereinstimmen" = 4, "weitgehend übereinstimmen" = 3,
            "weitgehend ablehnen" = 2, "voll und ganz ablehnen" = 1),
        pol_equal_treatment_2 = recode(pol_equal_treatment_2,
            "voll übereinstimmen" = 1, "weitgehend übereinstimmen" = 2,
            "weitgehend ablehnen" = 3, "voll und ganz ablehnen" = 4)
    ) %>%
    mutate(
        # Combine political satisfaction items into one scale
        pol_satisfaction = rowMeans(select(.,
            matches("pol_jus|pol_equ")), na.rm = TRUE)
    ) %>%
    # Exclusion according to the pre-registered exclusion criteria
    # criterion "a. who did not complete the probe recognition task" and
    # criterion "c. who self-report not being fluent in German"
    # are neccessarily not met due to the manner of data collection.
    # 170 --> 163 (7 excluded, read reasons below)
    filter(
        # b. who withdraw their consent to
        # data analysis after full debriefing
        consent_given == 1, # excludes 0
        # d. whose recognition performance
        # is < 60% in the probe recognition task
        correct_rate >= .6, # excludes 0
        # e. whose average response times are slower
        # than two standard deviations of the sample mean
        mean_rt < cutoff, # excludes 5
        # f. who self-report not having followed the instructions
        # conscientiously (5 or lower on a scale from 1 to 10) or
        # rate their own data to be unfit for analyses
        compliance > 5, # excludes 1
        data_quality. == "ja" # excludes 1
    ) %>%
    select(subject_id = ID, age, gender, pol_interest,
        pol_orientation, pol_satisfaction)

# ---------------------------

# Merge all data.frames
d_long <- trials %>%
    left_join(ratings) %>%
    filter(subject_id %in% subjects$subject_id) %>%
    left_join(subjects) %>%
    left_join(stims) %>%
    select(-rt)

# ---------------------------

# Export data
saveRDS(d_long, file = "Processed data/d-long.rds")
d_long <- readRDS("Processed data/d-long.rds")
