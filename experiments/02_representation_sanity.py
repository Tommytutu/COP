"""Run Proposition 1's n=3,...,9 computational feasibility sanity check."""

from _bootstrap import parser, pipeline


def main() -> None:
    args = parser("Representation-theorem computational sanity check").parse_args()
    pipeline(args.config).run_sanity()


if __name__ == "__main__":
    main()
