from __future__ import annotations

import re
from collections import defaultdict

import scipy
import numpy as np
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


def get_directions(indels_mut: np.ndarray, repeat_regions: list[tuple]) -> dict[tuple, bool]:
    directions = dict()
    for start, end in repeat_regions:
        error_counts = indels_mut[start:end]
        # If there are many errors at the beginning, direction is False.
        if error_counts[0] <= error_counts[-1]:
            directions[(start, end)] = True
        else:
            directions[(start, end)] = False
    return directions


def get_counts_homopolymer(
    indels_mut: np.ndarray, directions: dict[tuple, bool]
) -> dict[int, list[float], dict[tuple[int, int], np.ndarray]]:
    # Initialize default dictionaries to hold counts
    error_counts = defaultdict(list)
    error_counts_regions = dict()
    # Iterate through each repeat region
    for (start, end), direction in directions.items():
        indels_counts = indels_mut[start:end]
        # If there are many errors at the beginning, reverse the order
        if direction is False:
            indels_counts = indels_counts[::-1]
            start, end = end, start
        # Append the start, end, and total log errors to the region count
        error_counts_regions[(start, end)] = indels_counts
        # Append the log errors count to each position
        for position, value in enumerate(indels_counts):
            error_counts[position].append(value)
    return dict(error_counts), error_counts_regions


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


def return_thresholds(error_counts: dict[int, list[float]]) -> list[float]:
    # Initialize empty lists to hold x and y data
    error_positions = []
    error_counts = []
    # Populate the x and y data with the mutation counts and positions
    for position, counts in error_counts.items():
        position_values = [position] * len(counts)
        error_positions.extend(position_values)
        error_counts.extend(counts)
    # Convert x and y data to numpy arrays
    error_positions = np.array(error_positions)
    error_counts = np.array(error_counts)
    # Define a grid of x values from the minimum to the maximum position
    x_grid = np.linspace(error_positions.min(), error_positions.max(), error_positions.max() + 1)
    # Smooth the y data K times and stack the results
    num_smoothings = 100
    smoothed_data = np.stack([_smooth_data(error_positions, error_counts, x_grid) for _ in range(num_smoothings)]).T
    # Calculate the mean and standard error of the smoothed data
    mean_smoothed_data = np.nanmean(smoothed_data, axis=1)
    stderr_smoothed_data = scipy.stats.sem(smoothed_data, axis=1)
    stderr_smoothed_data = np.nanstd(smoothed_data, axis=1, ddof=0)
    # Define the thresholds
    thresholds = mean_smoothed_data + 1.90 * stderr_smoothed_data
    return thresholds


def get_errors_in_homopolyer(
    error_counts_regions: dict[tuple[int, int], np.ndarray], thresholds: list[float]
) -> set[int]:
    # Initialize a set to hold the locations of errors
    sequence_error_loci = set()
    # Iterate through each region and its associated error counts
    for (start, end), counts in error_counts_regions.items():
        # +-1 because 0-index is not considered as homopolymer
        if start > end:
            candidate_error_region = set(range(end, start - 1))
        else:
            candidate_error_region = set(range(start + 1, end))
        # Initialize a list to hold the indices of mutations
        mutation_indices = []
        # Iterate through each mutation count
        for index, count in enumerate(counts):
            # Skip the first count (since it has no previous count to compare to)
            if index == 0:
                continue
            # If the count exceeds the threshold, record the index
            if count > thresholds[index]:
                mutation_indices.append(index)
        # Add the absolute location of each mutation to the set of mutation locations
        mutations = set()
        for mutation_index in mutation_indices:
            if start > end:
                mutations.add(start - 1 - mutation_index)
            else:
                mutations.add(end + mutation_index)
        sequence_error_loci |= candidate_error_region - mutations
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
        directions = get_directions(indels_control[mut], repeat_regions)
        error_counts_control, _ = get_counts_homopolymer(indels_control[mut], directions)
        _, error_counts_regions_sample = get_counts_homopolymer(indels_sample[mut], directions)
        thresholds = return_thresholds(error_counts_control)
        errors_in_homopolymer[mut] = get_errors_in_homopolyer(error_counts_regions_sample, thresholds)
    return errors_in_homopolymer
