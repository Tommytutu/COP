"""MNVDM runtime distributions and certification-rate experiment."""

from _bootstrap import parser, pipeline


def main() -> None:
    command = parser("MNVDM runtime and certification experiment")
    command.add_argument("--limit", type=int, help="Optional number of dataset rows for a quick check")
    args = command.parse_args()
    experiment = pipeline(args.config)
    experiment.run_priority(limit=args.limit)
    experiment.run_summarize()


if __name__ == "__main__":
    main()
