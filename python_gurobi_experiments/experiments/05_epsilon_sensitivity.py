"""Enumerate 4,683 weak orders and reproduce the paper's epsilon table."""

from _bootstrap import parser, pipeline


def main() -> None:
    args = parser("Exhaustive small-instance epsilon sensitivity experiment").parse_args()
    pipeline(args.config).run_sensitivity()


if __name__ == "__main__":
    main()
