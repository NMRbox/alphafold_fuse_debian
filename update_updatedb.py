#!/usr/bin/env python3
import argparse
import fileinput
import logging
"""Remove mountpoint from updatedb scan"""

_logger = logging.getLogger(' ')
logging.basicConfig()
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('mountpoint', help="fuse system mountpoint")
parser.add_argument('-l', '--loglevel', default='WARN', help="Python logging level")
parser.add_argument('-e', '--update-config', default='/etc/updatedb.conf')

args = parser.parse_args()
_logger.setLevel(getattr(logging, args.loglevel))

pathline = None
with fileinput.input(args.update_config, inplace=True) as f:
    for line in f:
        sline = line.strip('\n')
        _logger.debug(sline)
        if sline.startswith('PRUNEPATHS'):
            pathline = sline
            if not args.mountpoint in sline:
                parts = sline.split('"', maxsplit=2)
                if len(parts) == 3:
                    newprune = f'{parts[0]}"{parts[1] + " " + args.mountpoint}"{parts[2]}'
                    print(newprune)
                    _logger.info(f"Updated {newprune}")
                else:
                    _logger.warning(f"{sline} did not split into 3 parts {parts}")
                    print(sline)
            else:
                print(sline)
            continue

        print(sline)
    if sline is None:
        raise ValueError(f"PRUNEPATHS not found in {args.update_config}")
