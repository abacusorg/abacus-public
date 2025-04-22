#!/usr/bin/env python3

from Abacus import abacus


def main():
    parser = abacus.default_parser()
    parser.add_argument(
        'parfn',
        help='The parameter file, relative to CONFIG_DIR',
        nargs='?',
        default='abacus.par2',
    )
    args = parser.parse_args()
    args = vars(args)

    args['config_dir'] = '.'

    retcode = abacus.run(**args)

    exit(retcode)


if __name__ == '__main__':
    main()
