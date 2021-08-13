map_short_segment_id_to_long_segment_id <- function(dt) {
    dt[, .(segmentID, longSegmentID = rleid(.SD[[2]]), call = .SD[[2]])]
}

#' @export
#' @importFrom "progress" progress_bar
long_segment_filter <- function(samples_list, calls_table) {
    new_calls <- copy(calls_table)
    pb <- progress_bar$new(format = "Long segment filter [:bar] :what :percent eta: :eta",
                          total = length(samples_list),
                          clear = FALSE,
                          width = 80)
    longest_label <- max(nchar(names(samples_list)))
    for (samplename_ in names(samples_list)) {
        pb$tick(tokens=list(what=paste(c(samplename_, rep(' ', longest_label - nchar(samplename_))), collapse="")))
        #loginfo("Long segment filter: %s", samplename_)
        sample_calls <- calls_table[, .SD, .SDcols = c("segmentID", paste0(samplename_, ".totalCN"))]
        mapping <- map_short_segment_id_to_long_segment_id(sample_calls)

        sample_data <- samples_list[[samplename_]]
        updates <- unique(sample_data[mapping, , on = "segmentID"][, .(segmentID, call=max(0, as.integer(round(median(total_cn))))), by = longSegmentID])
        new_calls[updates, (paste0(samplename_, ".totalCN")) := call, on = "segmentID"]
    }
    polymorphic <- new_calls[, apply(as.matrix(.SD), 1, function(row) !all(row[1] == row)), .SDcols = grep("^s.+\\.totalCN$", colnames(new_calls))]
    new_calls[, isPolymorphic := polymorphic]
    new_calls
}

#' Makes a lookup table to speed up the adjacent segments filter.
#' @export
make_breakpoint_test_lookup_table <- function(samples_list, compute_pval = FALSE) {
    breakpoint_test_lookup_table <- list()
    pb <- progress_bar$new(format = "Making lookup table [:bar] :what :percent eta: :eta",
                           total = length(samples_list),
                           clear = FALSE,
                           width = 80)
    longest_label <- max(nchar(names(samples_list)))
    for (samplename_ in names(samples_list)) {
        pb$tick(tokens=list(what=paste(c(samplename_, rep(' ', longest_label - nchar(samplename_))), collapse="")))
        #loginfo("Breakpoint testing: %s", samplename_)
        dt <- samples_list[[samplename_]]
        result <- dt[, .(med = median(total_cn)),
                     by = segmentID][,
                                     .(segmentID,
                                       samplename = samplename_,
                                       effect_size = abs(shift(med) - med))][!is.na(effect_size)]

        # Precompute the median total cn of pairs of adjacent segments, to quickly see the copy number
        # that would occur if they were merged. This speeds up the adjacent segments filter.
        dt[, mergeOddID := segmentID - (segmentID %% 2)+1]
        dt[, mergeEvenID := segmentID + (segmentID %% 2)]
        mergeOdds <- dt[, .(merged_median = median(total_cn)), by = mergeOddID]
        mergeEvens <- dt[, .(merged_median = median(total_cn)), by = mergeEvenID]
        dt[, c("mergeOddID", "mergeEvenID") := NULL]
        setnames(mergeOdds, old = "mergeOddID", new = "segmentID")
        setnames(mergeEvens, old = "mergeEvenID", new = "segmentID")
        merge_data <- rbind(mergeOdds, mergeEvens)[order(segmentID)][segmentID > 1]
        result[merge_data, merged_median := i.merged_median, on = "segmentID"]

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
apply_adj_filter_until_stable <- function(samples_list, samplename_, calls_dt, breakpoint_lookup, threshold, method = "clonal") {
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
        nrounds <- nrounds + 1
        post_adj <- adjacent_segment_filter_iteration(samplename_, copy(post_adj), samples_list[[samplename_]], breakpoint_lookup, threshold, method)
        current_failures <- get_failures(samplename_, post_adj, breakpoint_lookup)
        logdebug("Round %d, %d Failures: %s", nrounds, length(current_failures), current_failures)

        if (length(current_failures) == 0) {
            return (list(calls=post_adj, oscillators=list()))
        }

        # Check for stability - the failures repeat a previous configuration, or there are no failures
        for (configuration in prev_failures) {
            if (configuration$equals(current_failures)) {
                matching_iter <- configuration$iteration
                stabilised = TRUE
                logdebug("Adjacent segments filter has stabilised - matching configurations found at iterations %d and %d", matching_iter, nrounds)
            }
        }

        prev_failures <- listappend(prev_failures, failuresList(iteration = nrounds, indices = current_failures))
    }

    # Collect oscillators
    period <- nrounds - matching_iter + 1
    if (nrounds - matching_iter < 0) { # Identify the oscillators
        oscillating_rounds <- seq(matching_iter+1, nrounds)
        breakpoints_in_any_oscillating_round <- Reduce(union, lapply(prev_failures[oscillating_rounds], function(item) item$indices))
        breakpoints_in_common <- Reduce(intersect, lapply(prev_failures[oscillating_rounds], function(item) item$indices))
        oscillators <- units(sort(unique(setdiff(breakpoints_in_any_oscillating_round, breakpoints_in_common))))
        logdebug("Found %d oscillators - %s", length(oscillators), oscillators)
    } else { # No oscillators to find
        logdebug("Found 0 oscillators")
        oscillators <- vector("integer", 0L)
    }

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
                                              min_effect_size = 0.6, method = "clonal") {
    # 0) Preamble - make copies of data that should stay constant
    #    and subset that data to the current sample
    filtered_calls <- copy(calls_dt[, .SD, .SDcols = c("segmentID", "nsnp", paste0(samplename_, ".totalCN"))])
    sample_data <- copy(sample_dt[, .(segmentID, Index, total_cn)])
    setkey(filtered_calls, segmentID)
    setkey(sample_data, segmentID)
    sample_data <- sample_data[filtered_calls]

    # 1) Select changepoints
    changepoint_segments <- get_changepoint_segids(filtered_calls, paste0(samplename_, ".totalCN"))

    # 2) Test changepoints and collect failures (breakpoints where the step change is below the threshold)
    tests <- breakpoint_test_table[samplename==samplename_ & segmentID %in% changepoint_segments]
    failures <- tests[p_value >= 0.05 | effect_size < min_effect_size, segmentID]

    if (length(failures) == 0) {
        return (filtered_calls)
    }
    suggested_updates <- list()

    suggest_new_call <- function(bk, sample_lookup_table, method="clonal") {
        merged_median <- sample_lookup_table[segmentID==bk, merged_median]
        if (method == "clonal") {
            suggested_call <- round(merged_median)
        } else if (method == "subclonal") {
            suggested_call <- round(2*merged_median)/2
        } else { stop("Unrecognised method") }
        suggested_call
    }

    sample_lookup_table <- breakpoint_test_table[samplename==samplename_]
    for (bk in failures) {
        #suggested_call <- suggest_new_call(bk, samplename_, sample_data, breakpoint_test_table, filtered_calls, method)
        suggested_call <- suggest_new_call(bk, sample_lookup_table, method)
        suggestions <- data.table(segmentID = bk,
                                  suggested_new_call = suggested_call)
        suggested_updates <- listappend(suggested_updates, suggestions)
    }
    suggested_updates <- rbindlist(suggested_updates)

    # 3) Apply the changes
    filtered_calls[suggested_updates, (paste0(samplename_, ".totalCN")) := suggested_new_call]
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

# suggest_new_call <- function(breakpoint_id, samplename, sample_data, breakpoints_table, current_filtered_calls, method="clonal", search_range = 10000) {
#     if (isFALSE(method %in% c("clonal", "local_distribution", "subclonal"))) {
#         stop("method should be one of 'integer', 'halfinteger' or 'local_distribution'")
#     }
#
#     segments_to_merge <- c(min(breakpoint_id) - 1, breakpoint_id)
#     current_cns <- current_filtered_calls[segmentID %in% segments_to_merge, .SD[[1]], .SDcols = paste0(samplename, ".totalCN")]
#     if (all(current_cns[1] == current_cns)) {
#         return (current_cns[1])
#     }
#
#     data <- sample_data[segmentID %in% segments_to_merge, total_cn]
#
#     if (method == "local_distribution") {
#         breakpoint_position <- get_breakpoint_positions(breakpoint_id, breakpoints_table)
#         cns_to_check <- seq(min(current_cns), max(current_cns), 1)
#         local_sample_data <- sample_data[Index > (breakpoint_position - search_range) & Index < (breakpoint_position + search_range)]
#         local_medians <- local_sample_data[, .(median_total_cn=median(total_cn)), by = totalCN][order(totalCN)]
#         local_medians[, dist := abs(median_total_cn - median(data))]
#         logdebug(local_medians)
#         return (local_medians[dist == min(dist), totalCN])
#     } else if (method == "clonal") {
#         return (as.integer(round(median(data), 0)))
#     } else {
#         return (round(median(data*2), 0) / 2)
#     }
# }

resolve_oscillators <- function(samplename_, sample_dt, calls_dt, breakpoint_test_lookup_table, oscillators, threshold = 0.6) {
    .resolve_duo <- function(context, approved) {
        if (context$cn_left != context$cn_segment & abs(context$cn_left_median - context$cn_segment_median) > threshold) {
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
        lr_diff <- abs(context$cn_left_median - context$cn_right_median)

        if (context$cn_left != context$cn_right & lr_diff > threshold) {
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
        logdebug("Resolving oscillators: current unit = %s", unit)
        ctxt <- get_unit_context(unit, changepoints, sample_dt)
        if (is.null(ctxt$cn_right_median) | is.nan(ctxt$cn_right_median)) {
            logdebug("Resolving as a duo: %d-%d | %d-%d", ctxt$left, ctxt$start, ctxt$start, ctxt$end)
            approved <- .resolve_duo(ctxt, approved)
        } else {
            logdebug("Resolving as a trio: %d-%d | %d-%d | %d-%d", ctxt$left, ctxt$start, ctxt$start, ctxt$end, ctxt$end, ctxt$right)
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
        right_ <- min(sigbks[sigbks > end_])
    }

    cn_left_median_ <- sample_dt[(segmentID >= left_) & (segmentID < start_), median(total_cn)]
    cn_segment_median_ <- sample_dt[(segmentID >= start_) & (segmentID < end_), median(total_cn)]
    cn_right_median_ <- sample_dt[(segmentID >= end_) & (segmentID < right_),  median(total_cn)]
    cn_overall_median_ <- sample_dt[(segmentID >= left_) & (segmentID < right_),  median(total_cn)]

    list(start=start_, end=end_, left=left_, right=right_,
         cn_left_median=cn_left_median_,
         cn_segment_median=cn_segment_median_,
         cn_right_median=cn_right_median_,
         cn_overall_median=cn_overall_median_,
         cn_left=as.integer(round(cn_left_median_)),
         cn_segment=as.integer(round(cn_segment_median_)),
         cn_right=as.integer(round(cn_right_median_)),
         cn_overall=as.integer(round(cn_overall_median_)))
}



#' Runs the adjacent segments filter on each sample in the samples list.
#' @param samples_list List of data.tables for each sample. Data tables should
#' already be marked with segment IDs of the initial segmentation.
#' @param comp_table Table comparing the segment medians for each sample in the set.
#' Should contain entries for each sample in the samples list. Can be created using
#' `make_comparison_table(samples_list)`, and can be further filtered with
#' `pairwise_filter(comp_table)`.
#' @param breakpoint_lookup_table Table that allows for quick check of breakpoint
#' significance level. Create using `make_breakpoint_test_lookup_table(samples_list)`
#' @param threshold for effect size in breakpoint assessment
#' @export
adjacent_segments_filter <- function(samples_list, comp_table, breakpoint_lookup_table, samplenames = NULL, threshold = 0.6, method = "clonal") {
    new_calls <- copy(comp_table)
    if (is.null(samplenames)) {
        samplenames <- names(samples_list)
    } else {
        samplenames <- intersect(names(samples_list), samplenames)
    }
    pb <- progress_bar$new(format = "Adjacent segments filter [:bar] :what :percent eta: :eta",
                           total = length(samplenames),
                           clear = FALSE,
                           width = 80)
    longest_label <- max(nchar(samplenames))
    for (samplename_ in samplenames) {
        pb$tick(tokens=list(what=paste(c(samplename_, rep(' ', longest_label - nchar(samplename_))), collapse="")))
        # loginfo("Running adjacent segments filter on %s [%d/%d, %.2f%%]", samplename_, iter, N, 100*iter/N)
        postadj <- apply_adj_filter_until_stable(samples_list, samplename_, comp_table, breakpoint_lookup_table, threshold, method)
        resolve_osc <- resolve_oscillators(samplename_, samples_list[[samplename_]], postadj$calls, breakpoint_lookup_table, postadj$oscillators)
        new_calls[, eval(paste0(samplename_, ".totalCN")) := resolve_osc$calls[, .SD, .SDcols = paste0(samplename_, ".totalCN")]]
    }
    polymorphic <- new_calls[, apply(as.matrix(.SD), 1, function(row) !all(row[1] == row)), .SDcols = grep("^s.+\\.totalCN$", colnames(new_calls))]
    new_calls[, isPolymorphic := polymorphic]
    new_calls
}

#' @export
make_comparison_table <- function(samples_list) {
    comp <- samples_list[[1]][, .(nsnp=unique(nsnp)), by = segmentID]
    pb <- progress_bar$new(format = "Making comparison table [:bar] :what :percent eta: :eta",
                           total = length(samples_list) + 1,
                           clear = FALSE,
                           width = 80)
    longest_label <- max(nchar(names(samples_list)))
    for (samplename_ in names(samples_list)) {
        pb$tick(tokens=list(what=paste(c(samplename_, rep(' ', longest_label - nchar(samplename_))), collapse="")))
        #loginfo("Making comparison table: analysing %s", samplename_)
        comp[, paste0(samplename_, ".totalCN") := samples_list[[samplename_]][, unique(totalCN), by = segmentID][, V1]]
        comp[, paste0(samplename_, ".medianCN") := samples_list[[samplename_]][, median(total_cn), by = segmentID][, V1]]
    }

    pb$tick(tokens=list(what="Final comparison"))
    #loginfo("Making comparison table: comparing all samples")
    find_most_frequent <- function(v) {
        t = table(v)
        as.integer(names(t)[which.max(t)])
    }

    most_frequent_state <- apply(comp[, .SD, .SDcols = grep("s.+\\.totalCN$", colnames(comp))], 1, find_most_frequent)
    comp[, mode.totalCN := most_frequent_state]

    selected_state_medians <- sapply(1:nrow(comp), function(i) {
        mfs <- comp[i, mode.totalCN]
        index <- which(unlist(comp[i, .SD, .SDcols = grep("s.+\\.totalCN$", colnames(comp))]) == mfs)
        medians.index <- sub("totalCN$", "medianCN", names(index))
        medians <- unlist(comp[i, .SD, .SDcols = medians.index])
        median(medians)
    })

    comp[, median.of.mode.state := selected_state_medians]

    setcolorder(comp, c("segmentID", "nsnp", "mode.totalCN", "median.of.mode.state"))
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
    pb <- progress_bar$new(format = "Pairwise filter [:bar] :percent eta: :eta",
                           total = nrow(comp_table_copy),
                           clear = FALSE,
                           width = 80)
    for (i in 1:nrow(comp_table_copy)) {
        # if (i %% 100 == 0) {
        #     loginfo("Pairwise filter processing segment ID: %d",
        #             comp_table_copy[i, segmentID])
        # }
        pb$tick()
        new_calls[[i]] <- update_row(row_index=i)
    }
    new_calls <- rbindlist(new_calls)
    comp_table_copy[, colnames(new_calls) := new_calls]
    comp_table_copy[, isPolymorphic := apply(as.matrix(.SD[, 3:38]), 1, function(row) !all(row[1] == row))]
    comp_table_copy
}

#' @param copy_number_table = table of data with a row for each segment. A row contains the integer
#' copy number state of each sample, the median of the per-SNP copy number estimates for each sample,
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

    # 2nd step - identify the columns that hold the copy number medians and integer values,
    # then extract this data into 'cn_medians' and 'integer_calls'. Also identify the unique calls.
    median_columns <- grep("^s.+\\.medianCN$", colnames(table_row))
    integer_columns <- grep("^s.+\\.totalCN$", colnames(table_row))
    cn_medians <- unlist(table_row[, .SD, .SDcols = median_columns])
    integer_calls <- unlist(table_row[, .SD, .SDcols = integer_columns])
    unique_calls <- sort(unique(integer_calls))

    # If there's only one state present in the segment then no correction is possible,
    # so let's just return the original calls unchanged
    if (length(unique_calls) <= 1) {
        return(integer_calls)
    }

    # 3rd step - compute all the pairwise distances. Easiest to do this in one go,
    # and pull out the distances for different copy number states later.
    distances <- as.matrix(dist(cn_medians))

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
                # median copy number in the failing sample and each sample with the candidate copy number.
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
    samples_to_change <- sub("\\.medianCN$", ".totalCN", results$samplename)
    stopifnot(isTRUE(all.equal(unname(integer_calls[samples_to_change]), results$from))) # assert that these samples currently have the 'from' state we expect them to have!

    integer_calls[samples_to_change] <- results$to
    return(integer_calls)
}

#' Runs the pairwise filter on the samples given in comp_table. This filter compares
#' the copy number call for each sample to the call in the majority of samples. If
#' the copy number is different to the majority, then it must have a median total copy
#' number that differs from the majority median by at least `threshold`.
#' See documentation for `multiple_comparison_check` for details.
#' @param comp_table Table comparing the segment medians for each sample in the set.
#' Can be created using `make_comparison_table(samples_list)`
#' @param threshold float
#' @param majority float
#' @export
pairwise_filter <- function(comp_table, threshold = 0.8, majority = 0.8) {
    calls <- update_copynumber_calls_multiple_comparison_strategy(comp_table, threshold, majority)
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

#' @export
get_all_breakpoints <- function(comp_table, samplenames = NULL) {
    if (is.null(samplenames)) {
        samplenames <- sub("\\.totalCN$", "", grep("\\.totalCN", colnames(comp_table), value = TRUE))
        samplenames <- samplenames[samplenames != "mode"]
    }
    all_bks <- vector("integer")
    for (samplename_ in samplenames) {
        breakpoint_segids <- get_changepoint_segids(comp_table, samplename_)
        all_bks <- union(all_bks, breakpoint_segids)
    }
    return (sort(unique(all_bks)))
}

recall_genotyping_segments <- function(samplename_, test_segments, sample_dt, calls_dt) {
    updated_calls <- copy(calls_dt)
    updates <- sample_dt[segmentID %in% test_segments, .(new_call=as.integer(round(median(total_cn)))), by = segmentID]
    setkey(updates, segmentID)
    updated_calls[updates, (paste0(samplename_, ".totalCN")) := new_call, on = "segmentID"]
    return (updated_calls)
}

recall_genotyping_segments <- function(samplename_, test_segments, sample_dt, calls_dt, allow_half_integer=FALSE) {
    updated_calls <- copy(calls_dt)
    if (allow_half_integer) {
        updates <- sample_dt[segmentID %in% test_segments, .(new_call=0.5*round(median(2*total_cn))), by = segmentID]
    } else {
        updates <- sample_dt[segmentID %in% test_segments, .(new_call=as.integer(round(median(total_cn)))), by = segmentID]
    }
    setkey(updates, segmentID)
    updated_calls[updates, (paste0(samplename_, ".totalCN")) := new_call, on = "segmentID"]
    return (updated_calls)
}

#' @export
genotype_copynumber_calls <- function(samples_list, comp_table, breakpoint_lookup_table) {
    all_candidates <- get_all_breakpoints(comp_table, names(samples_list))
    new_calls <- copy(comp_table)
    samplenames <- names(samples_list)
    pb <- progress_bar$new(format = "Genotyping breakpoints [:bar] :what :percent eta: :eta",
                           total = length(samplenames),
                           clear = FALSE,
                           width = 80)
    longest_label <- max(nchar(samplenames))
    for (samplename_ in samplenames) {
        pb$tick(tokens=list(what=paste(c(samplename_, rep(' ', longest_label - nchar(samplename_))), collapse="")))
        #loginfo("Running genotyping on %s [%d/%d, %.2f%%]", samplename_, iter, N, 100*iter/N)

        candidates <- get_genotyping_breakpoints_for_sample(comp_table, all_candidates, samplename_, breakpoint_lookup_table)
        new_calls <- recall_genotyping_segments(samplename_, candidates, samples_list[[samplename_]], comp_table)
    }
    polymorphic <- new_calls[, apply(as.matrix(.SD), 1, function(row) !all(row[1] == row)), .SDcols = grep("^s.+\\.totalCN$", colnames(new_calls))]
    new_calls[, isPolymorphic := polymorphic]
    new_calls
}

#' @export
unify_medians_filter <- function(samples_list, calls_table, method = "clonal") {
    if (!(method %in% c("clonal", "subclonal"))) {
        stop("Method must be one of 'clonal' or 'subclonal'")
    }
    segid_map <- get_unified_segment_ids(samples_list, calls_table)

    new_calls <- copy(calls_table)
    samplenames <- names(samples_list)
    pb <- progress_bar$new(format = "Unify medians filter [:bar] :what :percent eta: :eta",
                           total = length(samplenames),
                           clear = FALSE,
                           width = 80)
    longest_label <- max(nchar(samplenames))
    for (samplename_ in samplenames) {
        pb$tick(tokens=list(what=paste(c(samplename_, rep(' ', longest_label - nchar(samplename_))), collapse="")))
        #loginfo("Running unify medians filter (%s) on %s [%d/%d, %.2f%%]", method, samplename_, iter, N, 100*iter/N)
        dt <- samples_list[[samplename_]]
        if (method == "clonal") {
            median_call <- dt[segid_map, , on = "segmentID"][, .(median_call=max(0, as.integer(round(median(total_cn))))), by = ulsID][segid_map, median_call, on = "ulsID"]
        } else {
            median_call <- dt[segid_map, , on = "segmentID"][, .(median_call=max(0, as.integer(0.5*round(median(2*total_cn))))), by = ulsID][segid_map, median_call, on = "ulsID"]
        }
        new_calls[, (paste0(samplename_, ".totalCN")) := median_call]
    }

    polymorphic <- new_calls[, apply(as.matrix(.SD), 1, function(row) !all(row[1] == row)), .SDcols = grep("^s.+\\.totalCN$", colnames(new_calls))]
    new_calls[, isPolymorphic := polymorphic]
    new_calls
}

get_unified_segment_ids <- function(samples_list, calls_table) {
    all_breakpoints <- get_all_breakpoints(calls_table, names(samples_list))

    tmp <- data.table(segmentID=1:max(calls_table$segmentID), is_break = FALSE)
    tmp[segmentID %in% all_breakpoints, is_break := TRUE]
    tmp[, ulsID := cumsum(is_break) + 1]
    tmp
}


do_segment_calculation <- function(sample_dt, selecter, region_bounds = NULL, use_tetraploid_copynumber = FALSE) {
    if (!is.null(region_bounds)) {
        if (!is.numeric(region_bounds)) {
            stop("Region bounds must be numeric")
        } else if (length(region_bounds) != 2) {
            stop ("Region bounds must be length 2")
        }

        minStart <- region_bounds[1]
        maxEnd <- region_bounds[2]
        sample_dt <- sample_dt[START >= minStart & END <= maxEnd]
    }

    setkey(sample_dt, CHROM, START, END)
    if (use_tetraploid_copynumber) {
        ol <- foverlaps(sample_dt, selecter)[, .(CHROM, START=i.START, END=i.END, total_cn=total_cn_tetraploid, selecterID, segmentID)]
    } else {
        ol <- foverlaps(sample_dt, selecter)[, .(CHROM, START=i.START, END=i.END, total_cn=total_cn_diploid, selecterID, segmentID)]

    }
    ol[, .(START=min(START), END=max(END), minSegmentID=min(segmentID), maxSegmentID=max(segmentID), medianCN=median(total_cn), clonalCN=round(median(total_cn)), subclonalCN=0.5*round(median(2*total_cn))), by=selecterID]
}

#' @export
get_changepoints_from_breakpoint_list <- function(segmentation_dt, sample_dt, breakpoints, method = "clonal", region_bounds = NULL, use_tetraploid_copynumber = FALSE) {
    # Activate only the breakpoints in the breakpoints list
    # ACTIVATION
    # 1) Divide the data into segments that run between each activated breakpoint
    # 2) For each segment:
    #   a) Calculate the median total copy number
    #   b) Assign either a clonal or subclonal copy number call based on the median
    # 3) Retrieve only activated breakpoints that also induce a copy number change
    # RESOLUTION
    # 4) Divide the data into segments that run between each activated changepoint
    # 5) For each segment:
    #   a) Calculate the median total copy number
    #   b) Assign either a clonal or subclonal copy number call based on the median

    if (!(method %in% c("clonal", "subclonal"))) {
        stop("Method should be one of 'clonal' or 'subclonal'")
    }

    # 1) & 2)
    segment_selecter <- make_selecter(segmentation_dt, breakpoints)
    activated_segment_stats <- do_segment_calculation(sample_dt, segment_selecter, region_bounds, use_tetraploid_copynumber = use_tetraploid_copynumber)

    # 3)
    if (method == "clonal") {
        changepoints <- activated_segment_stats[clonalCN != shift(clonalCN, type = "lag"), c(1, minSegmentID)]
    } else {
        changepoints <- activated_segment_stats[subclonalCN != shift(subclonalCN, type = "lag"), c(1, minSegmentID)]
    }

    # 4) & 5)
    segment_selecter <- make_selecter(segmentation_dt, changepoints)
    changepoint_segment_stats <- do_segment_calculation(sample_dt, segment_selecter, region_bounds, use_tetraploid_copynumber = use_tetraploid_copynumber)

    return (list(active_breakpoints = breakpoints,
                 changepoints = changepoints,
                 active_breakpoint_stats = activated_segment_stats,
                 changepoint_stats = changepoint_segment_stats,
                 method = method,
                 region = region_bounds))
}

#' Updates the copynumber calls in calls_dt for the subset of samples in samplenames,
#' using the provided set of breakpoints
update_calls <- function(calls_dt, samples_list, samplenames, breakpoints, segmentation_dt, use_tetraploid_copynumber = FALSE, method = "subclonal") {
    if (!(method %in% c("clonal", "subclonal"))) {
        stop("Method should be one of 'clonal' or 'subclonal'")
    }

    updated_calls <- copy(calls_dt)

    selecter <- make_selecter(segmentation_dt, breakpoints)
    pb <- progress_bar$new(format = "Updating calls [:bar] :what :percent eta: :eta",
                           total = length(samplenames),
                           clear = FALSE,
                           width = 80)
    longest_label <- max(nchar(samplenames))
    for (samplename_ in samplenames) {
        pb$tick(tokens=list(what=paste(c(samplename_, rep(' ', longest_label - nchar(samplename_))), collapse="")))
        new_calls <- make_calls_for_sample(samples_list[[samplename_]], selecter, use_tetraploid_copynumber, method)
        updated_calls[new_calls, (paste0(samplename_, ".totalCN")) := i.call, on = "segmentID"]
    }

    return (updated_calls)
}

median_shift_approved_breakpoints <- function(samplenames, segmentation_dt, samples_list, input_breakpoints, use_tetraploid_copynumber = FALSE, median_threshold = 0.3, method = "subclonal") {
    if (!(method %in% c("clonal", "subclonal"))) {
        stop("Method should be one of 'clonal' or 'subclonal'")
    }

    # For each sample, get the size of the median shift at each changepoint
    pb <- progress_bar$new(format = "Median shift filter: scanning samples [:bar] :what :percent eta: :eta",
                           total = length(samplenames),
                           clear = FALSE,
                           width = 80)
    longest_label <- max(nchar(samplenames))
    bks <- vector("integer")
    sample_bks <- list()
    for (samplename_ in samplenames) {
        changepoints <- get_changepoints_from_breakpoint_list(segmentation_dt,
                                                              samples_list[[samplename_]],
                                                              input_breakpoints,
                                                              method,
                                                              region_bounds = NULL,
                                                              use_tetraploid_copynumber = use_tetraploid_copynumber)
        changepoints$changepoint_stats[, medianShift := abs(medianCN - shift(medianCN, type = "lag"))]
        approved <- changepoints$changepoint_stats[is.na(medianShift) | medianShift > median_threshold, minSegmentID]
        sample_bks[[samplename_]] <- approved
        bks <- c(bks, approved)
        pb$tick(tokens=list(what=paste(c(samplename_, rep(' ', longest_label - nchar(samplename_))), collapse="")))
    }
    approved_breaks <- unique(sort(bks))
    return (list(all=approved_breaks, by_sample=sample_bks))
}

#' @export
median_shift_filter <- function(samplenames, segmentation_dt, samples_list, calls_to_filter, calls_to_update, use_tetraploid_copynumber, median_threshold = 0.3, method = "subclonal") {
    if (!(method %in% c("clonal", "subclonal"))) {
        stop("Method should be one of 'clonal' or 'subclonal'")
    }

    # Approve the breakpoints using the median filter
    approved_tet_bks <- median_shift_approved_breakpoints(samplenames, segmentation_dt, samples_list, get_all_breakpoints(calls_to_filter), use_tetraploid_copynumber = use_tetraploid_copynumber, median_threshold, method)
    calls_median_filtered <- copy(calls_to_update)

    pb <- progress_bar$new(format = "Median shift filter: updating calls [:bar] :what :percent eta: :eta",
                           total = length(samplenames),
                           clear = FALSE,
                           width = 80)
    longest_label <- max(nchar(samplenames))

    for (samplename_ in samplenames) {
        calls_median_filtered <- update_calls(calls_median_filtered, samples_list, samplename_, approved_tet_bks$by_sample[[samplename_]], segmentation_dt, use_tetraploid_copynumber = use_tetraploid_copynumber, method)
        pb$tick(tokens=list(what=paste(c(samplename_, rep(' ', longest_label - nchar(samplename_))), collapse="")))
    }
    return (calls_median_filtered)
}

#' 1) Take a subset of samples
#' 2) Scan subset for potential new subclonal breakpoints
#' 3) Update the calls table for the subset (possibly augmented subset) using the new breakpoints
#' 4) Median-Shift filter the calls
#' 5) Adjacent-Segments filter the calls (optional)
#' @export
subclonal_search_pipeline <- function(samples_to_search, samples_to_update, input_calls, lookup_table, samples_list, segmentation_dt, effect_size_threshold = 0.3, update_method = "subclonal", do_adjacent_segments_filter = FALSE) {
    existing_breakpoints <- get_all_breakpoints(input_calls)
    new_breakpoints <- get_subclonal_breakpoints(input_calls,
                                                 lookup_table,
                                                 samples_to_search=samples_to_search,
                                                 samples_to_exclude=NULL,
                                                 effect_size_threshold)
    union_of_breakpoints <- sort(union(existing_breakpoints, new_breakpoints))

    augmented_calls <- update_calls(input_calls,
                                    samples_list,
                                    samples_to_update,
                                    union_of_breakpoints,
                                    segmentation_dt,
                                    use_tetraploid_copynumber = FALSE,
                                    method = update_method)

    filtered_calls <- median_shift_filter(samples_to_update,
                                          segmentation_dt,
                                          samples_list,
                                          augmented_calls,
                                          input_calls,
                                          use_tetraploid_copynumber = FALSE,
                                          median_threshold = effect_size_threshold,
                                          method = update_method)

    if (do_adjacent_segments_filter) {
        filtered_calls <- adjacent_segments_filter(samples_list,
                                                   filtered_calls,
                                                   lookup_table,
                                                   samplenames = samples_to_update,
                                                   threshold = effect_size_threshold,
                                                   method = update_method)
    }

    return (filtered_calls)
}

#' Uses TOST (two one-sided tests) to check if the copy number data within two adjacent segments
#' differs by at least `tol`, at a significance level of `alpha`.
#' If this is true, then the breakpoint will be retained in the "approved_breaks" list, otherwise
#' it is removed from the segmentation, and the segmentation is rebuilt.
tost_prefilter <- function(dta, seg, tol=0.3, alpha=0.05) {
    # Returns TRUE if the distributions a and b are equivalent within tolerance `tol`
    tost.equivalent <- function(a, b, tol = 0.5, alpha = 0.05) {
        test1 <- t.test(a, b, mu = tol, var.equal=TRUE, alternative = "less")$p.value
        test2 <- t.test(a, b, mu = -tol, var.equal=TRUE, alternative = "greater")$p.value
        test1 < alpha & test2 < alpha
    }

    dta[, bk := rep(seq(1, seg$nIntervals), times = seg$Lengde)]
    dta[, segment_median := median(total_cn), by = bk]

    nseg <- dta[, max(bk)]

    if (nseg > 1) {
        tost_results <- sapply(1:(nseg-1), function(i) {
            result <- tost.equivalent(dta[segment_id==i, total_cn],
                                      dta[segment_id==(i+1), total_cn],
                                      tol, alpha)
            return (result)
        })
        approved_breaks <- (2:nseg)[!tost_results]
    } else {
        approved_breaks <- 1
    }

    dta[is_breakpoint_start, tmp := 1, on = c("START")]
    dta[, filtered_segment_id := cumsum(tmp) + 1]
    dta[, tmpI := .I]

    result <- dta[, .(Lengde=.N, sta=min(tmpI), mean=mean(total_cn)), by = filtered_segment_id]
    result <- list(Lengde = result[, Lengde],
                   sta = result[, sta],
                   mean = result[, mean],
                   nIntervals = dta[, max(filtered_segment_id)])

    dta[, tmp := NULL]
    dta[, tmpI := NULL]
    dta[, unfiltered_segment_id := bk]
    dta[, bk := NULL]
    dta[, segment_median := NULL]

    return (result)
}
