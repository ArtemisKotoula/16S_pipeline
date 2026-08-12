#!/usr/bin python

########################################################

import argparse
from scipy.stats import binom


def calculate_min_reads(frequency, confidence, min_taxon_reads):
    """
    Calculate the minimum total number of reads required so that
    P(X >= min_taxon_reads) >= confidence, where
    X ~ Binomial(N, frequency).
    """

    # Start at the point where the expected number of taxon reads
    # reaches min_taxon_reads
    n = int(min_taxon_reads / frequency)

    while True:
        probability = binom.sf(
            min_taxon_reads - 1,
            n,
            frequency
        )

        if probability >= confidence:
            return n

        n += 1


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Calculate the minimum sequencing depth required to "
            "detect a taxon at a specified frequency with a given "
            "confidence and minimum number of reads."
        )
    )

    parser.add_argument(
        "--frequency",
        type=float,
        required=True,
        help="Expected taxon frequency, e.g. 0.01 for 1%%"
    )

    parser.add_argument(
        "--confidence",
        type=float,
        required=True,
        help="Required detection probability, e.g. 0.95 for 95%%"
    )

    parser.add_argument(
        "--min_taxon_reads",
        type=int,
        required=True,
        help="Minimum number of reads required for the taxon"
    )

    args = parser.parse_args()

    # Validate arguments
    if not 0 < args.frequency < 1:
        parser.error("--frequency must be between 0 and 1")

    if not 0 < args.confidence < 1:
        parser.error("--confidence must be between 0 and 1")

    if args.min_taxon_reads < 1:
        parser.error("--min_taxon_reads must be >= 1")

    # Calculate cutoff
    min_reads_required = calculate_min_reads(
        frequency=args.frequency,
        confidence=args.confidence,
        min_taxon_reads=args.min_taxon_reads
    )

    print(min_reads_required)


if __name__ == "__main__":
    main()