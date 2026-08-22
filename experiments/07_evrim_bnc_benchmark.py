"""Compare direct convex MIQCP with lazy-OA Branch-and-Check EVRIM."""

from _bootstrap import parser, pipeline


def main() -> None:
    args = parser("Equal-settings EVRIM direct-versus-OA/B&C benchmark").parse_args()
    pipeline(args.config).run_bnc()


if __name__ == "__main__":
    main()
