## ---------------------------
## Script name: analyses.r
##
## Purpose of script:
##  Code used to perform the analyses
##  reported in the paper.
##
## Author: Carsten Sander
##
## Date Created: 2022-02-23
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
library(effsize) # effsize_0.8.1
library(tidyverse) # tidyverse_1.3.1
library(ggpubr) # ggpubr_0.4.0
library(broom) # broom_0.7.12

# Read data
d_long <- readRDS("Processed data/d-long.rds")

# ---------------------------

# Prepare data for by-participant rt analysis
d_rt_long <- d_long %>%
    filter(is_correct == 1) %>%
    select(-starts_with("rating_"), -item_id, -label, -is_correct, -rt) %>%
    group_by(
        subject_id, probe_type, age, gender,
        pol_interest, pol_orientation, pol_satisfaction
    ) %>%
    summarize(across(everything(), mean, na.rm = TRUE)) %>%
    ungroup()
d_rt <- d_rt_long %>%
    pivot_wider(
        names_from = probe_type,
        values_from = starts_with("rt_"),
        names_sep = "_"
    ) %>%
    mutate(rt_M2SD_none_diff = rt_M2SD_none_im - rt_M2SD_none_io)




#### sample characteristics
# We collected data from a total of N = 170 participants.
# Following our pre-registered exclusion criteria, we excluded
# five participants whose average response times were slower
# than two standard deviations over the sample mean, one who
# self-reported not having followed the instructions conscientiously
# and one who rated their own data to be unfit for analysis.
# This resulted in a sample of N = 163 participants (84 female,
# 76 male, 3 other; average age M = 29.9 years, SD = 10.4, ranging
# from 18 to 69). Participants were recruited via the online platform
# Prolific (www.prolific.co) and received monetary compensation of
# 3.13 GBP for completing the 25-minute study.

# On average, participants reported to be rather left leaning
# (M = 3.9, SD = 1.7 on a scale from 1 = left to 10 = right), rather
# interested in politics (M = 6.6, SD = 2.2, on a scale ranging from
# 1 = not at all to 10 = very strongly), and moderately satisfied
# with the German political system (M = 2.3, SD = 0.5, on a scale
# ranging from 1 = satisfied to 4 = dissatisfied).
summarise_factor <- function(factor, percentage = FALSE) {
    out <- ""
    tab <- table(factor)
    props <- prop.table(tab)
    n <- length(tab)
    for (i in seq(n)) {
        if (percentage == FALSE) {
            out <- paste0(out, tab[i], " ", names(tab)[i])
        } else {
            out <- paste0(out, round(props[i] * 100), "% ", names(props)[i])
        }
        if (i < n) {
            out <- paste0(out, ", ")
        }
    }

    return(out)
}
summarise_normal <- function(var, range = TRUE) {
    out <- paste0(
        "M = ", round(mean(var), 1), ", ",
        "SD = ", round(sd(var), 1)
    )
    if (range == TRUE) {
        out <- paste0(out, ", ",
        "ranging from ", min(var), " to ", max(var))
    }
    return(out)
}
paste0("average age ", summarise_normal(d_rt$age))
paste0(summarise_factor(d_rt$gender))
paste0(summarise_normal(d_rt$pol_orientation, FALSE))
paste0(summarise_normal(d_rt$pol_interest, FALSE))
paste0(summarise_normal(d_rt$pol_satisfaction, FALSE))




#### pre-registered one-sided paired t-test
# differences don't seem severely non-normal upon visual inspection
hist(d_rt$rt_M2SD_none_diff)

# On average, participants responded significantly slower in the
# implied condition (M = 0.81, SD = 0.12) than in the implied other
# condition (M = 0.75, SD = 0.10), t(162) = 16.59, p < .001, dz = 1.30.
t.test(d_rt$rt_M2SD_none_im, d_rt$rt_M2SD_none_io,
       paired = TRUE, alternative = "greater")
boxplot(d_rt$rt_M2SD_none_diff)

paste0("d_z = ", round(mean(d_rt$rt_M2SD_none_im - d_rt$rt_M2SD_none_io) /
       sd(d_rt$rt_M2SD_none_im - d_rt$rt_M2SD_none_io), 2))
paste0("implied: M = ", round(mean(d_rt$rt_M2SD_none_im), 2),
       ", SD = ", round(sd(d_rt$rt_M2SD_none_im), 2))
paste0("implied other: M = ", round(mean(d_rt$rt_M2SD_none_io), 2),
       ", SD = ", round(sd(d_rt$rt_M2SD_none_io), 2))

# visualization
plot <- ggboxplot(d_rt_long, x = "probe_type", y = "rt_M2SD_none",
    palette = "jco") +
    stat_compare_means(paired = TRUE, method = "t.test",
    comparisons = list(c("im", "io")), label = "p.format") +
    rremove("legend") +
    scale_x_discrete(labels = c("implied", "implied other")) +
    theme(panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA))
plot <- ggpar(plot,
    xlab = "Probe type",
    ylab = "Mean reaction time")
plot
ggexport(plot, width = 4, height = 4,
    filename = "Outputs/H1 Boxplot.pdf")

#### multiverse analysis
# select all variables containing reaction time means
d_sub <- d_rt %>%
    select(starts_with("rt")) %>%
    pivot_longer(
        cols = everything(),
        names_to = c("del", "cutoff", "trans", "probe_type"),
        values_to = "rt",
        names_sep = "_"
    ) %>%
    select(-del)
# loop over all combinations of cutoffs and transformations
multi_v <- data.frame(
    cutoff = rep(c("none", "M2SD", "f25", "f20", "f15"), 3),
    trans = c(rep("none", 5), rep("log", 5), rep("inv", 5)))
for (c in seq(nrow(multi_v))) {
    # filter data using c-th cutoff and transformation
    d_sub_c <- d_sub %>%
        filter(cutoff == multi_v$cutoff[c] & trans == multi_v$trans[c])
    # calculate test and effect size
    t_test_c <- t.test(d_sub_c$rt[d_sub_c$probe_type == "im"],
                       d_sub_c$rt[d_sub_c$probe_type == "io"],
                       paired = TRUE, alternative = "greater")
    d_z_c <- mean(d_sub_c$rt[d_sub_c$probe_type == "im"] -
                  d_sub_c$rt[d_sub_c$probe_type == "io"]) /
                  sd(d_sub_c$rt[d_sub_c$probe_type == "im"] -
                  d_sub_c$rt[d_sub_c$probe_type == "io"])
    # save test statistics
    multi_v$t[c] <- t_test_c$statistic
    multi_v$df[c] <- t_test_c$parameter
    multi_v$p[c] <- t_test_c$p.value
    multi_v$d_z[c] <- d_z_c
}

# The multiverse analysis showed the effect to be robust
# against various cutoff criteria and transformations.
write.csv(multi_v, "Outputs/multiverse.csv")




#### exploratory error rate analysis
# data preparation
d_er_long <- d_long %>%
    filter(is_correct == 1) %>%
    group_by(subject_id, probe_type) %>%
    summarize(correct_responses = n()) %>%
    mutate(er = (24 - correct_responses) / 24) %>%
    select(-correct_responses)
d_er <- d_er_long %>%
    pivot_wider(names_from = probe_type, values_from = er)

# differences seem severely non-normal upon visual inspection
# therefore a wilcoxon signed-rank test is performed
hist(d_er$im - d_er$io)

# Error rates were significantly higher in the implied condition (Mdn = 0.04)
# than in the implied other condition (Mdn = 0), p < .001, r = -.44.
wcx <- wilcox.test(d_er$im, d_er$io,
                   paired = TRUE, correct = TRUE, alternative = "greater")
wcx
r_from_wilcox <- function(wilcox_model, n) {
    z <- qnorm(wilcox_model$p.value / 2)
    r <- z / sqrt(n)
    cat(wilcox_model$data.name, "Effect Size, r = ", r)
}
r_from_wilcox(wcx, nrow(d_er) * 2)
paste0("implied: Mdn = ", round(median(d$error_rate_im), 2))
paste0("implied other: Mdn = ", round(median(d$error_rate_io), 2))

# visualization
plot <- ggboxplot(d_er_long, x = "probe_type", y = "er",
    palette = "jco") +
    stat_compare_means(paired = TRUE, method = "wilcox.test",
    comparisons = list(c("im", "io")), label = "p.format") +
    rremove("legend") +
    scale_x_discrete(labels = c("implied", "implied other")) +
    theme(panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA))
plot <- ggpar(plot,
    xlab = "Probe type",
    ylab = "Error rate")
plot
ggexport(plot, width = 4, height = 4,
    filename = "Outputs/H1 Error Rate Boxplot.pdf")



#### correlation between error rates and reaction times
# effects in error rate and reaction times are uncorrelated
d_participant <- left_join(d_rt, d_er)
cor.test(d_participant$rt_M2SD_none_im - d_participant$rt_M2SD_none_io,
         d_participant$im - d_participant$io)




#### exploratory regression analysis (by participant)
d_rt[, "rt_M2SD_none_diff"] <- d_rt$rt_M2SD_none_im - d_rt$rt_M2SD_none_io
mod <- lm(rt_M2SD_none_diff ~ pol_interest + pol_orientation + pol_satisfaction, d_rt) # nolint
summary(mod)
par(mfrow = c(2, 2))
plot(mod, 4)
d_rt_mod1 <- augment(mod) %>%
    mutate(index = seq(nrow(d_rt))) %>%
ggplot(d_rt_mod1, aes(pol_orientation, rt_M2SD_none_diff)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = pol_orientation, yend = .fitted), color = "red", size = 0.3)
##




#### exploratory by-item analyses
stim_data <- read.csv("Stimulus Data/stimuli.csv")
# Prepare for by-item analysis aggregated for each item_id
# thus comparing rts for the two different labels of each
# item (one implied and one implied other)
d_item_by_id <- d_long %>%
    mutate(across(starts_with("rt_"), ~ifelse(is_correct == 0, NA, .))) %>%
    select(item_id, label, probe_type, matches("rt_|rating")) %>%
    group_by(item_id, probe_type, label) %>%
    summarize(across(matches("rt_|rating"), ~mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    pivot_wider(
        names_from = probe_type,
        values_from = matches("rt_|rating|label"),
        names_sep = "_") %>%
    mutate(
        diff_M2SD_none = rt_M2SD_none_im - rt_M2SD_none_io,
        diff_rating_availability = rating_availability_im - rating_availability_io, # nolint
        diff_rating_identity = rating_identity_im - rating_identity_io,
        diff_rating_valence = rating_valence_im - rating_valence_io) %>%
    select(-starts_with("rt_"), -starts_with("rating")) %>%
    left_join(stim_data)

apa.cor.table(select(d_item_by_id, starts_with("diff_"), sconsensus, cconsensus, label_score), filename = "Outputs/by_item_by_id.doc", show.conf.interval = F) # nolint

# Prepare for by-item analysis aggregated for each label
# thus comparing rts for the same label in two conditions
# (but presented after different statements)
d_item_by_label <- d_long %>%
    mutate(across(starts_with("rt_"), ~ifelse(is_correct == 0, NA, .))) %>%
    select(label, probe_type, matches("rt_|rating")) %>%
    group_by(label, probe_type) %>%
    summarize(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_wider(
        names_from = probe_type,
        values_from = starts_with("r"),
        names_sep = "_") %>%
    ungroup() %>%
    mutate(diff_M2SD_none = rt_M2SD_none_im - rt_M2SD_none_io) %>%
    select(
        label, diff_M2SD_none,
        rating_availability = rating_availability_im,
        rating_identity = rating_identity_im,
        rating_valence = rating_valence_im) %>%
    left_join(stim_data) %>%
    mutate(behavior_length = unlist(lapply(str_extract_all(.$behavior, "[[:alpha:]]+"), length))) # nolint

apa.cor.table(select(d_item_by_label, matches("diff_|rating|consensus|score|behavior_length")), filename = "Outputs/by_item_by_label.doc", show.conf.interval = F) # nolint
attach(d_item_by_label)
plot(rating_availability, diff_M2SD_none, main = "Scatterplot Example", xlab = "Availability", ylab = "SII Effect", pch = 19) # nolint
plot(rating_identity, diff_M2SD_none, main = "Scatterplot Example", xlab = "Identification", ylab = "SII Effect", pch = 19) # nolint
plot(rating_valence, diff_M2SD_none, main = "Scatterplot Example", xlab = "Valence", ylab = "SII Effect", pch = 19) # nolint
plot(behavior_length, diff_M2SD_none, main = "Scatterplot Example", xlab = "Behavior Length", ylab = "SII Effect", pch = 19) # nolint
detach(d_item_by_label)

#### exploratory analyses of ratings effects
d_rate <- d_long %>%
    mutate(across(starts_with("rt_"), ~ifelse(is_correct == 0, NA, .))) %>%
    select(subject_id, item_id, probe_type, matches("rt_|rating")) %>%
    group_by(subject_id, item_id) %>%
    pivot_wider(
        names_from = probe_type,
        values_from = matches("rt_|rating"),
        names_sep = "_") %>%
    ungroup() %>%
    select(
        subject_id, item_id, rt_M2SD_none_im, rt_M2SD_none_io,
        rating_availability_im, rating_availability_io,
        rating_identity_im, rating_identity_io,
        rating_valence_im, rating_valence_io)

d_rate_binned <- data.frame()
for (subject in unique(d_rate$subject_id)) {
    # select subject and shuffle rows
    d_sub <- d_rate[d_rate$subject_id == subject, ]
    d_sub <- d_sub[sample(seq(nrow(d_sub))), ]
    original_order <- order(d_sub$item_id)
    # get ascending order of ratings
    order_availability <- order(d_sub$rating_availability_im)
    order_identity <- order(d_sub$rating_identity_im)
    order_valence <- order(d_sub$rating_valence_im)
    # split into three equal-sized sets (low, med, high)
    bins <- as.factor(c(rep(1, 8), rep(2, 8), rep(3, 8)))
    levels(bins) <- c("low", "med", "high")
    d_sub[order_availability, "availability_bin"] <- bins
    d_sub[order_identity, "identity_bin"] <- bins
    d_sub[order_valence, "valence_bin"] <- bins
    # restore original order and save
    d_sub <- d_sub[original_order, ]
    d_rate_binned <- dplyr::bind_rows(d_rate_binned, d_sub)
}

d_identity <- d_rate_binned %>%
    select(matches("subject|rt_|identity_bin")) %>%
    group_by(subject_id, identity_bin) %>%
    summarize(across(everything(), mean, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(rt_diff = rt_M2SD_none_im - rt_M2SD_none_io)
boxplot((d_identity$rt_diff) ~ d_identity$identity_bin)

d_availability <- d_rate_binned %>%
    select(matches("subject|rt_|availability_bin")) %>%
    group_by(subject_id, availability_bin) %>%
    summarize(across(everything(), mean, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(rt_diff = rt_M2SD_none_im - rt_M2SD_none_io)
boxplot((d_availability$rt_diff) ~ d_availability$availability_bin)

d_valence <- d_rate_binned %>%
    select(matches("subject|rt_|valence_bin")) %>%
    group_by(subject_id, valence_bin) %>%
    summarize(across(everything(), mean, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(rt_diff = rt_M2SD_none_im - rt_M2SD_none_io)
boxplot((d_valence$rt_diff) ~ d_valence$valence_bin)
