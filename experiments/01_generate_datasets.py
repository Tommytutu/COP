"""Generate the low/moderate/high/cyclic PCM dataset used by the paper."""

from _bootstrap import parser, pipeline


def main() -> None:
    args = parser("Generate all synthetic PCMs and latent weights").parse_args()
    pipeline(args.config).run_generate()


if __name__ == "__main__":
    main()
