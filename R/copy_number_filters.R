map_short_segment_id_to_long_segment_id <- function(dt) {
    dt[, .(segmentID, longSegmentID = rleid(.SD[[2]]), call = .SD[[2]])]
}

#' @export
long_segment_filter <- function(samples_list, calls_table) {
    new_calls <- copy(calls_table)
    for (samplename_ in names(samples_list)) {
        loginfo("Long segment filter: %s", samplename_)
        sample_calls <- calls_table[, .SD, .SDcols = c("segmentID", paste0(samplename_, ".totalCN"))]
        mapping <- map_short_segment_id_to_long_segment_id(sample_calls)

        sample_data <- samples_list[[samplename_]]
        updates <- unique(sample_data[mapping, , on = "segmentID"][, .(segmentID, call=as.integer(round(mean(total_cn)))), by = longSegmentID])
        new_calls[updates, (paste0(samplename_, ".totalCN")) := call, on = "segmentID"]
    }
    polymorphic <- new_calls[, apply(as.matrix(.SD), 1, function(row) !all(row[1] == row)), .SDcols = grep("^s.+\\.totalCN$", colnames(new_calls))]
    new_calls[, isPolymorphic := polymorphic]
    new_calls
}

#' @export
make_breakpoint_test_lookup_table <- function(samples_list, compute_pval = FALSE) {
    breakpoint_test_lookup_table <- list()
    for (samplename_ in names(samples_list)) {
        loginfo("Breakpoint testing: %s", samplename_)
        dt <- samples_list[[samplename_]]
        result <- dt[, .(mn = mean(total_cn)),
                     by = segmentID][,
                                     .(segmentID,
                                       samplename = samplename_,
                                       effect_size = abs(shift(mn) - mn))][!is.na(effect_size)]

        if (compute_pval) {
            nseg <- dt[, max(segmentID)]
            pvals <-
                sapply(2:nseg, function(i) {
                    wilcox.test(dt[segmentID == i, total_cn], dt[segmentID == i-1, total_cn])$p.value
                })
            result[, p_value := pvals]
        } else {
            result[, p_value := 0]
        }

        breakpoint_test_lookup_table <- listappend(breakpoint_test_lookup_table, result)
    }
    breakpoint_test_lookup_table <- rbindlist(breakpoint_test_lookup_table)
    setcolorder(breakpoint_test_lookup_table, c("p_value", "effect_size", "samplename", "segmentID"))
    breakpoint_test_lookup_table
}


# Filters used in processing total copy number
apply_adj_filter_until_stable <- function(samples_list, samplename_, calls_dt, breakpoint_lookup, threshold) {
    # Define a class
    failuresList <- setRefClass("failuresList",
                                fields=list(iteration = "numeric", indices = "numeric"),
                                methods=list(
                                    equals = function(v) {
                                        return (isTRUE(all.equal(v, indices)))
                                    }
                                ))

    nrounds <- 1
    stabilised <- FALSE
    post_adj <- copy(calls_dt)

    prev_failures <- list()
    prev_failures <- listappend(prev_failures, failuresList(iteration = nrounds, indices = get_failures(samplename_, calls_dt, breakpoint_lookup)))

    while(isFALSE(stabilised)) {
        post_adj <- adjacent_segment_filter_iteration(samplename_, copy(post_adj), samples_list[[samplename_]], breakpoint_lookup, threshold)
        nrounds <- nrounds + 1
        current_failures <- get_failures(samplename_, post_adj, breakpoint_lookup)

        if (length(current_failures) == 0) {
            return (list(calls=post_adj, oscillators=list()))
        }

        # Check for stability - the failures repeat a previous configuration, or there are no failures
        for (configuration in prev_failures) {
            if (configuration$equals(current_failures)) {
                matching_iter <- configuration$iteration
                stabilised = TRUE
                logdebug(paste("Match at", matching_iter))
            }
        }

        prev_failures <- listappend(prev_failures, failuresList(iteration = nrounds, indices = current_failures))
    }

    oscillating_rounds <- seq(matching_iter, nrounds-1)
    oscillators <- units(sort(unique(Reduce(c, lapply(prev_failures[oscillating_rounds], function(item) item$indices)))))

    return(list(calls=post_adj, oscillators=oscillators))
}

get_failures <- function(samplename_, calls_dt, breakpoint_lookup, min_effect_size = 0.6) {
    changepoint_segments <- get_changepoint_segids(calls_dt, paste0(samplename_, ".totalCN"))

    # 2) Test changepoints and collect failures (Note: this can be done once on every breakpoint and cached, right?)
    tests <- breakpoint_lookup[samplename==samplename_ & segmentID %in% changepoint_segments]
    failures <- tests[p_value >= 0.05 | effect_size < min_effect_size, segmentID]
    return (failures)
}

#' @export
get_changepoint_segids <- function(calls_dt, samplename) {
    if (!endsWith(samplename, ".totalCN")) {
        samplename <- paste0(samplename, ".totalCN")
    }
    dt <- calls_dt[, .SD, .SDcols = c("segmentID", "nsnp", samplename)]
    setnames(dt, old = samplename, new = "call")
    dt[call != shift(call), segmentID]
}

#' @export
get_breakpoint_positions <- function(segment_ids, breakpoint_table) {
    breakpoint_table[segmentID %in% segment_ids, start]
}


adjacent_segment_filter_iteration <- function(samplename_, calls_dt, sample_dt, breakpoint_test_table,
                                              min_effect_size = 0.6) {
    # 0) Preamble - make copies of data that should stay constant
    #    and subset that data to the current sample
    filtered_calls <- copy(calls_dt[, .SD, .SDcols = c("segmentID", "nsnp", paste0(samplename_, ".totalCN"))])
    sample_data <- copy(sample_dt[, .(segmentID, Index, total_cn)])
    setkey(filtered_calls, segmentID)
    setkey(sample_data, segmentID)
    sample_data <- sample_data[filtered_calls]

    # 1) Select changepoints
    changepoint_segments <- get_changepoint_segids(filtered_calls, paste0(samplename_, ".totalCN"))

    # 2) Test changepoints and collect failures (Note: this can be done once on every breakpoint and cached, right?)
    tests <- breakpoint_test_table[samplename==samplename_ & segmentID %in% changepoint_segments]
    failures <- tests[p_value >= 0.05 | effect_size < min_effect_size, segmentID]

    if (length(failures) == 0) {
        return (filtered_calls)
    }
    suggested_updates <- list()
    for (unit in units(failures)) {
        suggested_integer <- suggest_new_call(unit, samplename_, sample_data, breakpoints_table, filtered_calls)
        suggestions <- data.table(segmentID = c(min(unit) - 1, unit),
                                  integer_method = suggested_integer)
        suggested_updates <- listappend(suggested_updates, suggestions)
    }
    suggested_updates <- rbindlist(suggested_updates)

    # 3) Apply the changes
    filtered_calls[suggested_updates, (paste0(samplename_, ".totalCN")) := integer_method]
}

#' Break a vector into contiguous units
#' @example
#' units(c(655, 657, 667, 672, 684, 686, 687, 694, 709, 710))
#' # equiv. to list(655, 657, 667, 672, 684, c(686, 687), c(709, 710))
units <- function(v) {
    N <- length(v)
    if (N < 1) stop("Vector is empty")
    out <- list()
    current_unit <- v[1]
    if (N > 1) {
        for (i in 2:N) {
            elem <- v[i]
            if (elem - current_unit[length(current_unit)] == 1) {
                current_unit <- append(current_unit, elem)
            } else {
                out <- listappend(out, current_unit)
                current_unit <- elem
            }
        }
    }
    out <- listappend(out, current_unit)
    out
}

suggest_new_call <- function(breakpoint_id, samplename, sample_data, breakpoints_table, current_filtered_calls, search_range = 10000, method="integer") {
    if (isFALSE(method %in% c("integer", "local_distribution"))) {
        stop("method should be one of 'integer' or 'local_distribution")
    }

    segments_to_merge <- c(min(breakpoint_id) - 1, breakpoint_id)
    current_cns <- current_filtered_calls[segmentID %in% segments_to_merge, .SD[[1]], .SDcols = paste0(samplename, ".totalCN")]
    if (all(current_cns[1] == current_cns)) {
        return (current_cns[1])
    }

    data <- sample_data[segmentID %in% segments_to_merge, total_cn]

    if (method == "local_distribution") {
        breakpoint_position <- get_breakpoint_positions(breakpoint_id, breakpoints_table)
        cns_to_check <- seq(min(current_cns), max(current_cns), 1)
        local_sample_data <- sample_data[Index > (breakpoint_position - search_range) & Index < (breakpoint_position + search_range)]
        local_means <- local_sample_data[, .(mean_total_cn=mean(total_cn)), by = totalCN][order(totalCN)]
        local_means[, dist := abs(mean_total_cn - mean(data))]
        logdebug(local_means)
        return (local_means[dist == min(dist), totalCN])
    } else {
        return (as.integer(round(mean(data),0)))
    }
}

resolve_oscillators <- function(samplename_, sample_dt, calls_dt, breakpoint_test_lookup_table, oscillators) {
    .resolve_duo <- function(context, approved) {
        if (context$cn_left != context$cn_segment & abs(context$cn_left_mean - context$cn_segment_mean) > 0.6) {
            updated_calls[segmentID >= unit[1] & segmentID < unit[length(unit)], (paste0(samplename_, ".totalCN")) := context$cn_segment]
            tmp <- data.table(segmentID=unit, approved=FALSE)
            tmp[segmentID==context$start, approved := TRUE]
            approved <- listappend(approved, tmp)
        } else {
            updated_calls[segmentID >= unit[1] & segmentID < unit[length(unit)], (paste0(samplename_, ".totalCN")) := context$cn_left]
            approved <- listappend(approved, data.table(segmentID=unit, approved=FALSE))
        }
        return (approved)
    }

    .resolve_trio <- function(context, approved) {
        # This function closes over variables sample_dt, unit
        lr_diff <- abs(context$cn_left_mean - context$cn_right_mean)

        if (context$cn_left != context$cn_right & lr_diff > 0.6) {
            # Find a breakpoint within the current unit
            scores <- sapply(unit, function(bk) {
                score_merge_with_left <-
                    sample_dt[segmentID >= unit[1] & segmentID < bk, sum(abs(total_cn - context$cn_left))]
                score_merge_with_right <-
                    sample_dt[segmentID >= bk & segmentID < unit[length(unit)], sum(abs(total_cn - context$cn_right))]
                score <- score_merge_with_left + score_merge_with_right
                logdebug(sprintf("bk=%d, score_l=%.3f, score_r=%.3f, score=%.3f\n",
                                 bk, score_merge_with_left, score_merge_with_right, score))
                score
            })

            position_of_best_score <- which.min(scores)
            best_breakpoint <- unit[position_of_best_score]

            # Everything to the left of the best breakpoint is assigned the same CN as the left region,
            # and everything to the right gets the same CN as the right region
            updated_calls[segmentID >= unit[1] & segmentID < best_breakpoint, (paste0(samplename_, ".totalCN")) := context$cn_left]
            updated_calls[segmentID >= best_breakpoint & segmentID < unit[length(unit)], (paste0(samplename_, ".totalCN")) := context$cn_right]

            # Add best breakpoint to approved list
            tmp <- data.table(segmentID=unit, approved=FALSE)
            tmp[segmentID==best_breakpoint, approved := TRUE]
            approved <- listappend(approved, tmp)
        } else {
            if (context$cn_left != context$cn_right) {
                logdebug(sprintf("WARNING: mismatched left and right integer copy numbers for %d-%d and %d-%d",
                                 context$left, context$start, context$end, context$right))
            }
            updated_calls[segmentID >= context$left & segmentID < context$right, (paste0(samplename_, ".totalCN")) := context$cn_overall]

            # No breakpoints are approved
            approved <- listappend(approved, data.table(segmentID=unit, approved=FALSE))
        }
        return (approved)
    }

    changepoints <- c(calls_dt[, min(segmentID)],
                      intersect(get_changepoint_segids(calls_dt, paste0(samplename_, ".totalCN")),
                                get_significant_breakpoints(breakpoint_test_lookup_table, samplename_, 0.6)),
                      calls_dt[, max(segmentID)])

    updated_calls <- copy(calls_dt)

    approved <- list(data.table(segmentID=changepoints, approved=TRUE))

    for (unit in oscillators) {
        ctxt <- get_unit_context(unit, changepoints, sample_dt)
        if (is.null(ctxt$cn_right_mean) | is.nan(ctxt$cn_right_mean)) {
            approved <- .resolve_duo(ctxt, approved)
        } else {
            approved <- .resolve_trio(ctxt, approved)
        }
    }

    approved <- rbindlist(approved)
    approved <- rbind(approved, data.table(segmentID=setdiff(calls_dt$segmentID, approved$segmentID), approved=FALSE))[order(segmentID)]

    list(calls=updated_calls, approved_breakpoints=approved)
}

get_significant_breakpoints <- function(lookup_table, samplename_, effect_size_threshold, p_value_threshold = 0.05) {
    lookup_table[samplename == samplename_ & effect_size > effect_size_threshold & p_value < p_value_threshold, segmentID]
}

get_unit_context <- function(unit, sigbks, sample_dt) {
    start_ <- unit[1]
    end_ <- ifelse(length(unit) == 1, start_+1, unit[length(unit)])
    left_ <- max(sigbks[sigbks < start_])

    if (end_ >= max(sigbks)) {
        right_ <- end_
    } else {
        right_ <- min(sigbks[sigbks >= end_])
    }

    cn_left_mean_ <- sample_dt[(segmentID >= left_) & (segmentID < start_), mean(total_cn)]
    cn_segment_mean_ <- sample_dt[(segmentID >= start_) & (segmentID < end_), mean(total_cn)]
    cn_right_mean_ <- sample_dt[(segmentID >= end_) & (segmentID < right_),  mean(total_cn)]
    cn_overall_mean_ <- sample_dt[(segmentID >= left_) & (segmentID < right_),  mean(total_cn)]

    list(start=start_, end=end_, left=left_, right=right_,
         cn_left_mean=cn_left_mean_,
         cn_segment_mean=cn_segment_mean_,
         cn_right_mean=cn_right_mean_,
         cn_overall_mean=cn_overall_mean_,
         cn_left=as.integer(round(cn_left_mean_)),
         cn_segment=as.integer(round(cn_segment_mean_)),
         cn_right=as.integer(round(cn_right_mean_)),
         cn_overall=as.integer(round(cn_overall_mean_)))
}



#' Runs the adjacent segments filter on each sample in the samples list.
#' @param samples_list List of data.tables for each sample. Data tables should
#' already be marked with segment IDs of the initial segmentation.
#' @param comp_table Table comparing the segment means for each sample in the set.
#' Should contain entries for each sample in the samples list. Can be created using
#' `make_comparison_table(samples_list)`, and can be further filtered with
#' `pairwise_filter(comp_table)`.
#' @param breakpoint_lookup_table Table that allows for quick check of breakpoint
#' significance level. Create using `make_breakpoint_test_lookup_table(samples_list)`
#' @param threshold for effect size in breakpoint assessment
#' @export
adjacent_segments_filter <- function(samples_list, comp_table, breakpoint_lookup_table, threshold = 0.6) {
    new_calls <- copy(comp_table)
    for (samplename_ in names(samples_list)) {
        loginfo(samplename_)
        postadj <- apply_adj_filter_until_stable(samples_list, samplename_, comp_table, breakpoint_lookup_table, threshold)
        resolve_osc <- resolve_oscillators(samplename_, samples_list[[samplename_]], postadj$calls, breakpoint_lookup_table, postadj$oscillators)
        new_calls[, eval(paste0(samplename_, ".totalCN")) := resolve_osc$calls[, .SD, .SDcols = paste0(samplename_, ".totalCN")]]
    }
    polymorphic <- new_calls[, apply(as.matrix(.SD), 1, function(row) !all(row[1] == row)), .SDcols = grep("^s.+\\.totalCN$", colnames(new_calls))]
    new_calls[, isPolymorphic := polymorphic]
    new_calls
}

#' @export
make_comparison_table <- function(samples_list) {
    ploidy=2.0
    calc_total_copynumber <- cnpipe::calc_total_copynumber
    comp <- samples_list[[1]][, .(nsnp=unique(nsnp)), by = segmentID]

    for (samplename_ in names(samples_list)) {
        loginfo("Making comparison table: analysing %s", samplename_)
        comp[, paste0(samplename_, ".totalCN") := samples_list[[samplename_]][, unique(totalCN), by = segmentID][, V1]]
        comp[, paste0(samplename_, ".meanCN") := samples_list[[samplename_]][, mean(total_cn), by = segmentID][, V1]]
    }

    loginfo("Making comparison table: comparing all samples")
    find_most_frequent <- function(v) {
        t = table(v)
        as.integer(names(t)[which.max(t)])
    }

    most_frequent_state <- apply(comp[, .SD, .SDcols = grep("s.+\\.totalCN$", colnames(comp))], 1, find_most_frequent)
    comp[, mode.totalCN := most_frequent_state]

    selected_state_means <- sapply(1:nrow(comp), function(i) {
        mfs <- comp[i, mode.totalCN]
        index <- which(unlist(comp[i, .SD, .SDcols = grep("s.+\\.totalCN$", colnames(comp))]) == mfs)
        means.index <- sub("totalCN$", "meanCN", names(index))
        means <- unlist(comp[i, .SD, .SDcols = means.index])
        mean(means)
    })

    comp[, mean.of.mode.state := selected_state_means]

    setcolorder(comp, c("segmentID", "nsnp", "mode.totalCN", "mean.of.mode.state"))
}


update_copynumber_calls_multiple_comparison_strategy <- function(comp_table, threshold, majority) {
    stopifnot(majority > 0.0 & majority <= 1.0)
    stopifnot(threshold > 0.0)
    comp_table_copy <- copy(comp_table)
    update_row <- function(row_index) {
        new_calls <- multiple_comparison_check(comp_table_copy, row_index, threshold, majority)
        as.data.table(t(new_calls))
    }

    new_calls <- list()
    for (i in 1:nrow(comp_table_copy)) {
        if (i %% 100 == 0) {
            loginfo("Pairwise filter processing segment ID: %d",
                    comp_table_copy[i, segmentID])
        }
        new_calls[[i]] <- update_row(row_index=i)
    }
    new_calls <- rbindlist(new_calls)
    comp_table_copy[, colnames(new_calls) := new_calls]
    comp_table_copy[, isPolymorphic := apply(as.matrix(.SD[, 3:38]), 1, function(row) !all(row[1] == row))]
    comp_table_copy
}

#' @param copy_number_table = table of data with a row for each segment. A row contains the integer
#' copy number state of each sample, the mean of the per-SNP copy number estimates for each sample,
#' and the most frequent copy number state.
#' @param segment_id - this function only looks at one segment of the table
#' @param threshold - compared distances have to be greater than this value to pass
#' @param majority - the proportion of passes has to be greater than this value to count
#' @return integer_counts - vector of integer copy number estimates for each sample, possibly updated.
#'
#' Procedure:
#' Extract a row from the table of copy number data.
#' If there are multiple different copy number estimates among the samples, do this:
#' * for each distinct pair of copy number states, compare the state with the fewest samples
#'   (minority state) against the majority state. (i.e. measure the distance between each sample
#'   in the minority set against each sample in the majority set).
#'
#' * Count the fraction of sample comparisons that exceed the distance threshold, for each
#'   sample that's in the minority set.
#' * If this fraction is too small (fails the check), then the minority state is changed to
#'   the majority state, for this sample.
#'
#' If there is no majority (equal numbers of samples in each of the two copy number states)
#' we don't know which samples should be changed. So in this case we measure the fraction
#' of sample comparisons that exceed the distance threshold for both the minority and majority sets.
#' We also compute the distance among samples within the two sets, so that samples have the opportunity
#' to remain unchanged. Otherwise a likely result is that the two sets will simply swap states.
multiple_comparison_check <- function(copy_number_table, segment_id, threshold, majority) {
    # Some steps of preparation before doing comparisons...

    # 1st step - extract our target row from the table
    table_row <- copy_number_table[segmentID == segment_id]

    # 2nd step - identify the columns that hold the copy number means and integer values,
    # then extract this data into 'cn_means' and 'integer_calls'. Also identify the unique calls.
    mean_columns <- grep("^s.+\\.meanCN$", colnames(table_row))
    integer_columns <- grep("^s.+\\.totalCN$", colnames(table_row))
    cn_means <- unlist(table_row[, .SD, .SDcols = mean_columns])
    integer_calls <- unlist(table_row[, .SD, .SDcols = integer_columns])
    unique_calls <- sort(unique(integer_calls))

    # If there's only one state present in the segment then no correction is possible,
    # so let's just return the original calls unchanged
    if (length(unique_calls) <= 1) {
        return(integer_calls)
    }

    # 3rd step - compute all the pairwise distances. Easiest to do this in one go,
    # and pull out the distances for different copy number states later.
    distances <- as.matrix(dist(cn_means))

    # Preparation over. Now loop over the possible combinations of copy number states.
    # This should treat the case where there are onnly 2 unique states correctly.

    # Results go in this list.
    results <- list()

    # Loop over unique states in pairs - comparing state_i to state_j (e.g. 1+2, 1+3, 1+4, 2+3, 2+4, etc...)
    for (i in 1:(length(unique_calls)-1)) {
        state_i <- unique_calls[i]
        number_of_state_i <- length(integer_calls[integer_calls == state_i])

        for (j in (i+1):length(unique_calls)) {
            state_j <- unique_calls[j]
            number_of_state_j <- length(integer_calls[integer_calls == state_j])

            # Find the majority state or if there is no clear majority
            majority_state <- ifelse(number_of_state_i > number_of_state_j, state_i, state_j)
            minority_state <- ifelse(number_of_state_i > number_of_state_j, state_j, state_i)
            no_clear_majority <- number_of_state_i == number_of_state_j

            # Pull out the distances for state_i and state_j. Put the minority state in the rows.
            # If there is no clear majority, put state_i in the rows, state_j in the columns.
            # Meaning: if the minority state is CN=1, and the majority state is CN=2, then the rows
            # of 'dists' will be the samples with CN=1, and the columns will be the samples with CN=2.
            if (!no_clear_majority) {
                dists <- distances[which(integer_calls == minority_state),
                                   which(integer_calls == majority_state),
                                   drop = FALSE] # This is needed to force R to keep this as a matrix, even if it only has 1 row or column
            } else {
                dists <- distances[which(integer_calls == state_i),
                                   which(integer_calls == state_j),
                                   drop = FALSE]
            }

            # Convert the distances into TRUE/FALSE - TRUE meaning that the distance is above the theshold value.
            passing_threshold <- dists > threshold

            if (!no_clear_majority) {
                # Count the proportion of comparisons that pass the threshold, for the samples in the rows (which is the
                # minority sample unless there is no clear majority). Because we only want to change the copy number estimate
                # of the minority sample to the majority state, we only need to count along the rows.
                proportion_passing_threshold <- rowMeans(passing_threshold)
                failing_samples <- which(proportion_passing_threshold < majority)

                # Store the details of any failing samples. Don't make any copy number updates yet, save them for after the loop.
                # Possibly there will be some ambiguities to resolve. To help in this, also record the average distance between
                # mean copy number in the failing sample and each sample with the candidate copy number.
                if (length(failing_samples) > 0) {
                    result <- data.table(samplename =   names(failing_samples),
                                         from =         minority_state,
                                         to =           majority_state,
                                         avg_distance = rowMeans(dists[failing_samples, , drop=FALSE]))
                    results <- listappend(results, result)
                }
            } else  {
                # If there is no clear majority, we can't restrict ourselves to only considering changing the samples with
                # the minority copy number state. We need to look in the column directions as well.
                proportion_passing_threshold_row <- rowMeans(passing_threshold)
                failing_samples_row <- which(proportion_passing_threshold_row < majority)
                if (length(failing_samples_row) > 0) {
                    result_row <- data.table(samplename =   names(failing_samples_row),
                                             from =         state_i,
                                             to =           state_j,
                                             avg_distance = rowMeans(dists[failing_samples_row, , drop=FALSE]))
                    results <- listappend(results, result_row)
                }

                proportion_passing_threshold_col <- colMeans(passing_threshold)
                failing_samples_col <- which(proportion_passing_threshold_col < majority)
                if (length(failing_samples_col) > 0) {
                    result_col <- data.table(samplename =   names(failing_samples_col),
                                             from =         state_j,
                                             to =           state_i,
                                             avg_distance = colMeans(dists[, failing_samples_col, drop=FALSE]))
                    results <- listappend(results, result_col)
                }

                # And finally also allow for samples to stay in the original state, by looking at the average
                # distance between themselves and all the other samples in their copy number state
                dists_i <- distances[which(integer_calls == state_i),
                                     which(integer_calls == state_i),
                                     drop=FALSE]
                avg_i <- apply(dists_i, 1, function(row) {
                    mean(row[row>0])
                })

                result_i <- data.table(samplename =   names(avg_i),
                                       from =         state_i,
                                       to =           state_i,
                                       avg_distance = avg_i)
                results <- listappend(results, result_i)

                dists_j <- distances[which(integer_calls == state_j),
                                     which(integer_calls == state_j),
                                     drop=FALSE]
                avg_j <- apply(dists_j, 1, function(row) {
                    mean(row[row>0])
                })

                result_j <- data.table(samplename =   names(avg_j),
                                       from =         state_j,
                                       to =           state_j,
                                       avg_distance = avg_j)
                results <- listappend(results, result_j)
            }
        } # END OF LOOP OVER STATE j
    } # END OF LOOP OVER STATE i
    results <- rbindlist(results)

    # if results is empty, because no samples were flagged for change, return the original calls
    if (nrow(results) == 0) {
        return(integer_calls)
    }

    # To resolve any ambiguities, e.g. a sample that was flagged for a new copy number state in multiple comparisons,
    # (e.g. a sample has CN=3, but failed the majority vote compared to CN=2 and CN=4), simply choose the state with
    # the smallest avg_distance.
    results <- results[, .SD[which.min(avg_distance)], by = samplename]
    samples_to_change <- sub("\\.meanCN$", ".totalCN", results$samplename)
    stopifnot(isTRUE(all.equal(unname(integer_calls[samples_to_change]), results$from))) # assert that these samples currently have the 'from' state we expect them to have!

    integer_calls[samples_to_change] <- results$to
    return(integer_calls)
}

#' Runs the pairwise filter on the samples given in comp_table. This filter compares
#' the copy number call for each sample to the call in the majority of samples. If
#' the copy number is different to the majority, then it must have a mean total copy
#' number that differs from the majority mean by at least `threshold`.
#' See documentation for `multiple_comparison_check` for details.
#' @param comp_table Table comparing the segment means for each sample in the set.
#' Can be created using `make_comparison_table(samples_list)`
#' @param threshold
#' @param majority
#' @export
pairwise_filter <- function(comp_table, threshold = 0.8, majority = 0.8) {
    calls <- update_copynumber_calls_multiple_comparison_strategy(comp_table, 0.8, 0.8)
    polymorphic <- calls[, apply(as.matrix(.SD), 1, function(row) !all(row[1] == row)), .SDcols = grep("^s.+\\.totalCN$", colnames(calls))]
    calls[, isPolymorphic := polymorphic]
}

get_genotyping_breakpoints_for_sample <- function(comp_table, all_breakpoints, samplename_, breakpoint_lookup) {
    extra_breakpoints_in_sample <- get_significant_breakpoints(breakpoint_lookup, samplename_, 0.9)
    candidates <- sort(unique(c(all_breakpoints, extra_breakpoints_in_sample)))

    # Remove any existing changepoints from consideration
    existing_changepoints <- get_changepoint_segids(comp_table, samplename_)
    candidates <- setdiff(candidates, existing_changepoints)

    # Only look at breakpoints with large enough step size
    large_enough_step_size <- get_significant_breakpoints(breakpoint_lookup, samplename_, 0.6)
    candidates <- intersect(candidates, large_enough_step_size)

    return (sort(unique(candidates)))
}

get_all_breakpoints <- function(comp_table, samplenames) {
    all_bks <- vector("integer")
    for (samplename_ in samplenames) {
        breakpoint_segids <- get_changepoint_segids(comp_table, samplename_)
        all_bks <- union(all_bks, breakpoint_segids)
    }
    return (sort(unique(all_bks)))
}

recall_genotyping_segments <- function(samplename_, test_segments, sample_dt, calls_dt) {
    updated_calls <- copy(calls_dt)
    updates <- sample_dt[segmentID %in% test_segments, .(new_call=as.integer(round(mean(total_cn)))), by = segmentID]
    setkey(updates, segmentID)
    updated_calls[updates, (paste0(samplename_, ".totalCN")) := new_call, on = "segmentID"]
    return (updated_calls)
}

#' @export
genotype_copynumber_calls <- function(samples_list, comp_table, breakpoint_lookup_table) {
    all_candidates <- get_all_breakpoints(comp_table, names(samples_list))
    new_calls <- copy(comp_table)
    for (samplename_ in names(samples_list)) {
        candidates <- get_genotyping_breakpoints_for_sample(comp_table, all_candidates, samplename_, breakpoint_lookup_table)
        new_calls <- recall_genotyping_segments(samplename_, candidates, samples_list[[samplename_]], comp_table)
    }
    polymorphic <- new_calls[, apply(as.matrix(.SD), 1, function(row) !all(row[1] == row)), .SDcols = grep("^s.+\\.totalCN$", colnames(new_calls))]
    new_calls[, isPolymorphic := polymorphic]
    new_calls
}
