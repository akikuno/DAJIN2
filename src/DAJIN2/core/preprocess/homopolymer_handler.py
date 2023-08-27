from __future__ import annotations

import re
from collections import defaultdict
import numpy as np
import scipy
from statsmodels.nonparametric.smoothers_lowess import lowess as sm_lowess


def get_repeat_regions(sequence: str, loci: set[int]) -> list[tuple[int, int]]:
    """
    Find homopolymers in the sequence but discard them that
    are adjacent to candidate mutation loci because they are
    likely to be covered by the real mutations
    """
    pattern = r"A{4,}|C{4,}|G{4,}|T{4,}|N{4,}"
    repeat_regions = []
    for start, end in (match.span() for match in re.finditer(pattern, sequence)):
        if not (start - 1 in loci and end + 1 in loci):
            repeat_regions.append((start, end))
    return repeat_regions


def get_counts_homopolymer(
    indels_sample_mut, indels_control_mut, repeat_regions
) -> tuple[dict[int, list[float]], list[tuple[int, int, np.array]]]:
    # Initialize default dictionaries to hold counts
    mutation_counts = defaultdict(list)
    mutation_counts_regions = []
    # Iterate through each repeat region
    for start, end in repeat_regions:
        # Calculate mutations for each sample and control
        sample_mutations = np.array(indels_sample_mut[start:end])
        control_mutations = np.array(indels_control_mut[start:end])
        # Total mutations is the sum of sample and control
        total_mutations = sample_mutations + control_mutations
        # If the total number of mutations is greater in the last position, reverse the order
        if total_mutations[0] > total_mutations[-1]:
            total_mutations = total_mutations[::-1]
            start, end = end, start
        # Apply a log transformation to the total mutations
        total_mutations_log = np.log(total_mutations)
        # Append the start, end, and total log mutations to the region count
        mutation_counts_regions.append((start, end, total_mutations_log))
        # Append the log mutation count to each position
        for position, value in enumerate(total_mutations_log):
            mutation_counts[position].append(value)
    return mutation_counts, mutation_counts_regions


def _smooth_data(input_x, input_y, input_xgrid):
    # Sample 50 data points from the input x and y
    sampled_indices = np.random.choice(len(input_x), 50, replace=True)
    sampled_y = input_y[sampled_indices]
    sampled_x = input_x[sampled_indices]
    # Apply lowess smoothing to the sampled data
    smoothed_y = sm_lowess(sampled_y, sampled_x, frac=1.0 / 5.0, it=5, return_sorted=False)
    # Interpolate the smoothed data onto the input grid
    interpolated_y_grid = scipy.interpolate.interp1d(sampled_x, smoothed_y, fill_value="extrapolate")(input_xgrid)
    return interpolated_y_grid


def return_thresholds(mutation_counts) -> list(float):
    # Initialize empty lists to hold x and y data
    mutation_positions = []
    mutation_counts_log = []
    # Populate the x and y data with the mutation counts and positions
    for position, counts in mutation_counts.items():
        position_values = [position] * len(counts)
        mutation_positions.extend(position_values)
        mutation_counts_log.extend(counts)
    # Convert x and y data to numpy arrays
    mutation_positions = np.array(mutation_positions)
    mutation_counts_log = np.array(mutation_counts_log)
    # Define a grid of x values from the minimum to the maximum position
    x_grid = np.linspace(mutation_positions.min(), mutation_positions.max(), mutation_positions.max() + 1)
    # Smooth the y data K times and stack the results
    num_smoothings = 100
    smoothed_data = np.stack(
        [_smooth_data(mutation_positions, mutation_counts_log, x_grid) for _ in range(num_smoothings)]
    ).T
    # Calculate the mean and standard error of the smoothed data
    mean_smoothed_data = np.nanmean(smoothed_data, axis=1)
    stderr_smoothed_data = scipy.stats.sem(smoothed_data, axis=1)
    stderr_smoothed_data = np.nanstd(smoothed_data, axis=1, ddof=0)
    # Define the thresholds
    thresholds = mean_smoothed_data + 1.90 * stderr_smoothed_data
    return thresholds


def get_errors_in_homopolyer(mutation_counts_regions, thresholds) -> set(int):
    # Initialize a set to hold the locations of mutations
    sequence_error_loci = set()
    # Iterate through each region and its associated mutation counts
    for start, end, log_mutations in mutation_counts_regions:
        # +-1 because 0-index is not considered as homopolymer
        if start > end:
            region = set(range(end, start - 1))
        else:
            region = set(range(start + 1, end))
        # Initialize a list to hold the indices of mutations
        mutation_indices = []
        # Iterate through each mutation count
        for index, log_mutation in enumerate(log_mutations):
            # Skip the first count (since it has no previous count to compare to)
            if index == 0:
                continue
            # If the count exceeds the threshold, record the index
            if log_mutation > thresholds[index]:
                mutation_indices.append(index)
        # Add the absolute location of each mutation to the set of mutation locations
        mutations = set()
        for mutation_index in mutation_indices:
            if start > end:
                mutations.add(start - 1 - mutation_index)
            else:
                mutations.add(end + mutation_index)
        sequence_error_loci |= region - mutations
    return sequence_error_loci


###########################################################
# main
###########################################################


def extract_errors(sequence, indels_sample, indels_control, candidate_loci: dict[set]) -> dict[str, set(int)]:
    errors_in_homopolymer = dict()
    for mut in ["+", "-", "*"]:
        repeat_regions = get_repeat_regions(sequence, candidate_loci[mut])
        if repeat_regions == []:
            errors_in_homopolymer[mut] = set()
            continue
        mutation_counts, mutation_counts_regions = get_counts_homopolymer(
            indels_sample[mut], indels_control[mut], repeat_regions
        )
        thresholds = return_thresholds(mutation_counts)
        errors_in_homopolymer[mut] = get_errors_in_homopolyer(mutation_counts_regions, thresholds)
    return errors_in_homopolymer
