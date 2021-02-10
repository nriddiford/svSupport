#!/usr/bin/env python
from __future__ import division
import sys

from parseConfig import parse_config
from utils import make_dirs, cleanup
from getArgs import get_args
from worker import worker


def main():
    options, args = get_args()
    make_dirs(options.out_dir)

    if options.config:
        cleanup(options.out_dir)
        parse_config(options)
        sys.exit()

    elif options.test:
        print()
        print("* Running in test mode...")
        print()

        options.region = '3L:9892365-9894889'
        options.out_dir = 'test/test_out'
        options.in_file = 'test/data/test.bam'
        options.debug = True
        options.guess = True

    if options.in_file and options.region:
        try:
            worker(options)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return


if __name__ == "__main__":
    sys.exit(main())
