"""Run the paper's three protected-judgment house-buying EVRIM cases."""

from _bootstrap import parser, pipeline


def main() -> None:
    args = parser("Saaty house-buying EVRIM cases A/B/C").parse_args()
    pipeline(args.config).run_house()


if __name__ == "__main__":
    main()
