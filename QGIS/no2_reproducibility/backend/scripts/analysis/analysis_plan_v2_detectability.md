# Analysis plan v2: detectability of changes in weather-normalised TROPOMI NO2 (Algeria)

Written: 2026-11-20. Authors: M. Gana et al.
Status: written BEFORE any of the analyses below was run. Save this file next to analysis_plan_lockdown_normalisation.md and do not edit it afterwards. Later changes go in a dated "Deviations" section at the end.

## 0. What is already known (so this plan does not pretend to be blind)

All results of draft v1 are known: elasticities, the 2020 containment estimates (Table 8), the placebo windows (Table 9), gradient boosting (Table 10), trends (Table 11), the Friday cycle (Table 12), and the Level-2 replication with controls. The Friday analysis was exploratory and is NOT in plan v1. Plan v1 was written after the elasticity results and before the containment estimates. Nothing below changes how those existing estimates were made.

## 1. Central claim to be tested

After weather normalisation, a one-off change in the NO2 column in an 83-day window can be detected at a study site only if it is larger than a minimum detectable effect (MDE). The 2020 containment change is below that MDE at most sites. A recurring weekly modulation (Friday) is detectable at all sites where it exceeds a few percent.

## 2. Frozen choices (not to be changed in any analysis below)

- Sites: the six study sites and the two existing controls. 10 km disc, 50-100 km ring.
- Data: Level-3 daily columns as in draft v1; qa and coverage rules as in v1 (at least 30% valid pixels in disc and ring).
- Model: ln(NO2) on window dummy, linear trend, month, weekday, Ramadan, ln wind10, ln BLH, T2m; Newey-West errors, lag 7 (the "weather-adjusted" model of v1). The unadjusted model is reported alongside.
- Window length 83 days. Primary containment window 23 March to 13 June 2020.
- No change of model, window, radius or counterfactual method after seeing any result of this plan.

## 3. Analysis A: synthetic injection (the core of the paper)

1. Windows: all 83-day placebo windows used in v1 (start every 7 days inside sample A, outside the true window, at least 15 usable days). Also report the subset of non-overlapping windows separately.
2. For each site (study sites and controls) and each window, multiply the observed NO2 on the window days by (1 + d), with d in {0, -5%, -10%, -15%, -20%, -30%}. Fit the v1 weather-adjusted and unadjusted models and the gradient-boosting counterfactual exactly as in v1 (the true containment dummy stays in the model).
3. Report per site, method and d: mean recovered change, bias (recovered minus d), share of windows whose 95% interval excludes zero (power; at d = 0 this is the false-positive rate), and for gradient boosting the share outside its bootstrap interval.
4. MDE = smallest |d| on the grid with power >= 0.80 and bias within 20% of d. If power never reaches 0.80, report "above 30%".
5. Because windows overlap, state the effective number of windows and give results for the non-overlapping subset too.
6. Interpretation rule, fixed now: the 2020 change at a site is called "detected" only if its adjusted 95% interval excludes zero AND its absolute value exceeds that site's MDE. Otherwise it is called "not separable from ordinary variability".

## 4. Analysis B: seasonal placebo windows

1. Windows: 23 March to 13 June in each of 2021, 2022, 2023, 2024 and 2025 (years without containment).
2. Model: the frozen model fitted within each year on 1 January to 31 August of that year, with a window dummy; processor-version step dummies for any step date (as listed in v1 section 2.10) that falls inside that period.
3. The 2020 window is estimated with the same within-year model (1 January to 31 August 2020) as an additional estimate. The primary 2020 estimate in Table 8 is unchanged.
4. Report the five seasonal placebo estimates per site and where the 2020 value ranks among the six. Only five springs exist; say so. No p-value is claimed.

## 5. Analysis C: Friday cycle as a pre-registered prediction

Predictions written now, before opening any new data:

- Sign: negative Friday effect (weather + background model) at every urban or industrial site.
- Ranking of magnitude, from the source-mix reading of draft v1: Algiers and Oran largest; Constantine and Sétif next; Arzew smaller; Hassi Messaoud smallest.
- Confirmatory holdout 1: the six study sites over 1 January 2026 to the latest complete month available, same model. Prediction: negative sign with 95% interval excluding zero at Algiers, Oran, Constantine and Sétif; Arzew and Hassi Messaoud smaller in magnitude than those four (their intervals may include zero, as the sample is about 40 Fridays).
- Confirmatory holdout 2 (needs a site list fixed BEFORE downloading): additional Algerian cities chosen by a rule (population above a threshold, at least 100 km from the six sites), each classified as traffic-dominated or industrial from static layers only, before any NO2 is opened.
- Contrast (to be checked first): Morocco and Tunisia, where the weekend is Saturday-Sunday (verify against official sources before use). Prediction: no Friday dip beyond the regional background.
- Success rule: signs as predicted at the listed sites. Rank agreement is reported descriptively (Spearman), with no threshold, because n is small. Failure is reported as failure.

## 6. Analysis D: more background controls (limited purpose)

Purpose: estimate the false-positive rate of the gradient-boosting counterfactual and the regional Friday dip. Not used to build the null distribution for polluted sites (that is Analysis A).

- Rule: v1 eligibility rules (static layers only), up to 10 sites, at least 100 km from each study site and from each other, one per 0.25 degree cell, lowest night light first. Fixed and archived before any NO2 or wind data for them are opened.
- Report the number of sites and an effective number allowing for spatial correlation (correlation of their daily residuals).

## 7. Outputs

Tables: MDE per site and method (A); seasonal placebo estimates (B); Friday predictions versus outcomes (C); control false-positive rates (D). Scripts and CSVs are named detect_*.py / detect_*.csv.

## 8. Explicitly exploratory (to be labelled as such in the paper)

Anything not listed above, including wind-sector transport for Fridays (separate paper, separate plan), wind-rotation and EDGAR comparisons.

## 9. Open items for the author before running

- Choose the city list and population threshold for Analysis C (holdout 2), and confirm the Morocco/Tunisia weekend from official sources.
- Confirm Level-3 data for 2026 are available to download.
- Save the SHA-256 of this file (given in the chat) with the date.

## Deviations (add dated entries below; do not edit the text above)

(none yet)

### 2026-11-20 (entered after running Analyses A and B; plan text above is unchanged)

1. Analysis A, gradient boosting: no time-blocked cross-validation and no block bootstrap (computing cost). Detection threshold = 2.5th percentile of the gradient-boosting shifts of the same site at d = 0 (empirical null). The true containment window was excluded from the training days of every placebo fit.
2. Analysis A: two-sided false-positive rate at d = 0 was 3% to 25% (nominal 5%): Algiers 21%, Oran 12%, Arzew 3%, Constantine 4%, Setif 18%, Hassi Messaoud 14%, both controls 25%. Power values are therefore optimistic. MDE was still computed by the rule in section 3 (power as share of intervals below zero, at least 0.80; bias within 20% of d), on a coarse grid.
3. Analysis B: months September to December are absent from the January-August data, so their dummies were dropped. Any column linearly dependent on earlier ones was dropped (greedy rank check): the version step of 1 July 2021 equals the July dummy, and the step of 17 July 2022 was dropped or kept according to rank.
4. Analysis B: with one year of data and month dummies, the window dummy is weakly identified (window covers parts of months 3 to 6). The 2020 estimate of this model differs from Table 8 of draft v1 (different baseline); it is not used as the primary 2020 estimate.

### 2026-11-20 (protocol for Analyses C and D, fixed BEFORE any NO2 or wind data for these sites were opened)

**Analysis D.** select_extra_controls.py was run on the static candidate table (Control_candidates.csv, 578 candidates) with the v1 eligibility rule plus the two v1 amendments (latitude >= 34.5 N; 100 km <= distance to the nearest study site <= 150 km), 100 km minimum separation between all controls (including the two existing ones), lowest night light first. 23 candidates were eligible and only 4 new sites satisfy the separation: ctl_x01 (35.75 N, 2.25 E), ctl_x02 (35.75 N, 3.50 E), ctl_x03 (36.00 N, 8.00 E), ctl_x04 (35.25 N, 7.00 E). The plan said "up to 10"; four is what the frozen rule yields. Total controls: six. The effective number of controls will be reported with their spatial correlation.

**Analysis C, holdout 1.** Six study sites, 1 January to 14 November 2026 (ERA5 availability in Earth Engine to be checked; if ERA5 ends earlier, the holdout ends there and this is reported).

**Analysis C, holdout 2.** Cities chosen by the rule: wilaya capitals, agglomeration of about 150,000 or more (approximate figures from memory, to be checked against the census before downloading), at least 100 km from each of the six study sites (distances computed from approximate city-centre coordinates). Selected: Annaba (36.90 N, 7.77 E), Tlemcen (34.88, -1.32), Biskra (34.85, 5.73), Djelfa (34.67, 3.26), Chlef (36.17, 1.33), Tebessa (35.40, 8.12), Tiaret (35.37, 1.32). Prior classification, written now: Annaba = mixed with heavy industry (steel, fertiliser; predicted smaller Friday dip); the other six = traffic and residential (predicted dip as in Algiers, Oran, Constantine, Setif). Prediction: negative Friday effect (weather + background model) with a 95% interval excluding zero at the six traffic-and-residential cities; Annaba smaller in magnitude than the median of those six. Period 7 August 2019 to 14 November 2026. Morocco and Tunisia contrast: not selected yet (weekend to be verified from official sources first).


## Deviations and results, added 21 November 2026 (the text above is unchanged)

- Analysis D, gradient-boosting counterfactual: its false-positive rate at the four added sites (3.1%) is 2.5% by construction, because the threshold is the 2.5th percentile of the same site's shifts at d = 0. It is not an estimate of the method's error and is not reported as one.
- Ramadan 2026 was entered as 18 February to 19 March 2026 (approximate dates) in the Friday model.
- Analysis C, holdout 1 (run 21 November 2026 with detect_friday_holdout.py, model of friday_cycle.py unchanged): the prediction held. Ring-controlled Friday effect, 2026: Algiers -33.4% (-41.1 to -24.7), Oran -16.2% (-25.9 to -5.3), Constantine -14.0% (-20.8 to -6.5), Setif -15.4% (-22.3 to -8.0); Arzew -10.8% (-21.0 to +0.7) and Hassi Messaoud -4.8% (-10.2 to +0.9) smaller in magnitude.
- Analysis C, holdout 2 (run 21 November 2026 with detect_friday_cities.py, same model): the first part of the prediction held (negative effect with interval excluding zero at all six traffic-and-residential cities: Tlemcen -19.4%, Tiaret -10.6%, Djelfa -9.9%, Tebessa -9.7%, Biskra -9.5%, Chlef -9.4%). The second part FAILED: Annaba was -12.1% (-15.1 to -9.0), larger in magnitude than the median of the six (-9.8%).
- The Morocco and Tunisia contrast was not carried out.
- The city classification and the populations were from general knowledge and were not checked against census data.

## Deviations added 21 November 2026 (calibration of the inference; the text above is unchanged)

- Newey-West lag: `newey_west` (analyze_no2_ventilation.py, line 565) counts the lag in retained observations, after days without valid data are removed, so a lag of 7 spans 7 to about 10 calendar days. It was used unchanged in Table 8, the placebo windows, the injection and the Friday models.
- Calibration (detect_calibrate.py, results in detect_calibration_*.csv), added after the results: (a) the mean Newey-West standard error was 0.66 to 1.11 times the SD of the 64 placebo shifts at the six study sites (0.54 to 1.11 at the twelve sites); (b) calendar-day bandwidths of 28 and 90 days made the false-positive rate worse (up to 39% and 52%), not better; (c) with a threshold at the 2.5th percentile of the same site's shifts at d = 0, the adjusted MDE was unchanged at the six study sites (30, 30, 20, 20, 30, 20 for Algiers, Oran, Arzew, Constantine, Setif, Hassi Messaoud) and rose at the background sites; (d) intervals of the 2020 estimate +/- 1.96 placebo SD included zero at all six study sites, including Hassi Messaoud (-24.1% to +4.1%). The placebo windows overlap, so only about seven are independent; a block bootstrap was not run.
- The rule "a 2020 change is detected only if its interval excludes zero and its absolute value exceeds the MDE" (section 2.14 of the paper) is dropped. The MDE is a property of the test, not a threshold for one estimate.
- Exact false-positive rates (two-sided, adjusted, 64 windows) are 21.875% (Algiers), 12.5% (Oran), 3.1% (Arzew), 4.7% (Constantine), 18.75% (Setif), 14.1% (Hassi Messaoud); the summary file rounded them down.
