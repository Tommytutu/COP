"""Run EVRIM and EVRIM+T on high-noise and cyclic PCMs."""

from _bootstrap import parser, pipeline


def main() -> None:
    command = parser("Difficult-PCM EVRIM protected-judgment experiment")
    command.add_argument("--limit", type=int, help="Optional number of selected instances")
    args = command.parse_args()
    pipeline(args.config).run_evrim(limit=args.limit)


if __name__ == "__main__":
    main()
