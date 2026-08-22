"""Build all manuscript CSV summaries and the representability figure."""

from _bootstrap import parser, pipeline


def main() -> None:
    args = parser("Build manuscript tables and figures from raw results").parse_args()
    pipeline(args.config).run_summarize()


if __name__ == "__main__":
    main()
