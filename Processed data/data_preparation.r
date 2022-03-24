library(tidyverse)
setwd("/Users/carsten/Documents/Eigene Experimente/1a Probe Recog 1F/Data Analysis") # nolint

# import files from folder
filenames_qua <- list.files("Qualtrics Data",
                            pattern = "*.csv",
                            full.names = TRUE)
filenames_pav <- list.files("Pavlovia Data",
                            pattern = "*.csv",
                            full.names = TRUE)
d_pav_raw <- lapply(filenames_pav, read.csv)
d_qua_raw <- lapply(filenames_qua, read.csv)[[1]][-c(1, 2), ]

# function
create_long_df <- function(d_pav_raw, d_qua_raw) {
    # Creates a long data.frame with rts and ratings
    #
    # For each combination of five cutoff-criteria for slow outliers
    # (none, M + 2SD, 2500 ms, 2000 ms, 1500 ms), and three transformation
    # criteria (none, log, inverse) a column with reaction times is created
    # with the outliers being replaced by NAs.
    #
    # Grouping vars: subject_id, item_id, label
    # Filter vars: is_correct
    # IVs: probe_type, rating_availability, rating_valence, rating_identity
    # DVs: rt_none_none etc.
    #
    # @param d_pav_raw A list of all individual Pavlovia data.frames
    # @param d_qua_raw A data.frame containing the Qualtrics data
    #
    # @return d Data.frame containing rts and ratings

    #### Pavlovia Data
    trials <- data.frame()
    ratings <- data.frame()
    subjects <- data.frame()

    # loop over all individual data.frames
    for (i in seq(length(d_pav_raw))) {
        ## GET TRIALS
        trials_i <- d_pav_raw[[i]] %>%
            # select test trials
            filter(
                trial_type == "test",
                probe_type %in% c("implied", "implied_other")
            ) %>%
            select(
                subject_id = ID,
                item_id,
                label = probe,
                probe_type,
                is_correct = c_probe_resp.corr,
                rt = c_probe_resp.rt
            ) %>%
            mutate(
                probe_type = recode(
                    probe_type,
                    "implied" = "im",
                    "implied_other" = "io"))

        # get reaction times
        rts <- trials_i$rt[trials_i$is_correct == 1]

        # define matrix of cutoff/transformation combinations
        m2sd <- mean(rts) + 2 * sd(rts)
        corrections <- data.frame(
            cutoff = rep(c("none", "M2SD", "f25", "f20", "f15"), 3),
            cutoff_value = rep(c(Inf, m2sd, 2.5, 2.0, 1.5), 3),
            trans = c(rep("none", 5), rep("log", 5), rep("inv", 5)),
            trans_function = c(rep("{}", 5), rep("log({})", 5), rep("(1 / {}) * -1", 5))) # nolint

        # for each combination of cutoff criteria and transformations
        for (c in seq(NROW(corrections))) {
            # get values for current correction
            var_name <- paste0("rt_", corrections$cutoff[c], "_", corrections$trans[c]) # nolint
            cutoff_value <- corrections$cutoff_value[c]
            trans_function <- corrections$trans_function[c]
            trans_function <- gsub("\\{\\}", "rts", trans_function)

            # apply corrections to rts
            rts <- trials_i$rt
            rts[rts > cutoff_value] <- NA
            rts <- eval(parse(text = trans_function))

            # append data to data.frame
            trials_i[, var_name] <- rts
        }
        trials <- dplyr::bind_rows(trials, trials_i)

        ## GET SUMMARY DATA ON PARTICIPANTS
        subjects[i, "ID"] <- trials_i$subject_id[1]
        subjects[i, "correct_rate"] <- sum(trials_i$is_correct) / nrow(trials_i)
        subjects[i, "mean_rt"] <- mean(trials_i$rt_M2SD_none[trials_i$is_correct == 1], na.rm = TRUE) # nolint

        ## GET RATINGS
        ratings_i <- d_pav_raw[[i]] %>%
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

    #### Merge with Qualtrics Data
    cutoff <- mean(subjects$mean_rt) + 2 * sd(subjects$mean_rt)
    # get subject_ids of subjects that are not excluded based on qualtrics data
    subjects <- subjects %>%
        left_join(d_qua_raw, by = "ID") %>%
        filter(part == 2) %>%
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
            # combine political satisfaction items into one scale
            pol_satisfaction = rowMeans(select(.,
                matches("pol_jus|pol_equ")), na.rm = TRUE)
        ) %>%
        # exclusion according to the pre-registered exclusion criteria
        # criterion "a. who did not complete the probe recognition task" and
        # criterion "c. who self-report not being fluent in German"
        # are neccessarily not met due to the manner of data collection
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

    d <- trials %>%
        left_join(ratings) %>%
        filter(subject_id %in% subjects$subject_id) %>%
        left_join(subjects)
    return(d)
}

# summarize and export data
d_long <- create_long_df(d_pav_raw, d_qua_raw)
saveRDS(d_long, file = "Outputs/d_long.rds")
