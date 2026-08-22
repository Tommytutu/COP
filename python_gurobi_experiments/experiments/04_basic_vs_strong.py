"""Compare MNVLLSM-basic and MNVLLSM-strong Stage-1 formulations."""

from _bootstrap import parser, pipeline


def main() -> None:
    command = parser("Equal-settings basic-versus-strong formulation experiment")
    command.add_argument("--limit", type=int, help="Optional number of selected instances")
    args = command.parse_args()
    pipeline(args.config).run_formulation(limit=args.limit)


if __name__ == "__main__":
    main()
