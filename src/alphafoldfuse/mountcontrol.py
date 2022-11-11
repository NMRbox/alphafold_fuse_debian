#!/usr/bin/env python3
import argparse
import logging
import os.path
import subprocess
import sys
import time
from typing import Optional

import psutil
"""Manage mount point helper"""

_logger = logging.getLogger('mountcontrol')
_HELP = """umount: attempt to unmount system
query: display mounted state
forceumount: force unmount by killing processes using mount"""


class MountInfo:
    """Mount information"""

    def __init__(self, line):
        parts = line.strip('\n').split(maxsplit=3)
        self.name = parts[0]
        self.mountpoint = parts[1]
        self.type = parts[2]
        self.options = parts[3]

    @property
    def description(self):
        return f"{self.mountpoint} ({self.name}) is {self.type} mount with options {self.options}"


class MountControl:

    def __init__(self, path):
        if not os.path.isabs(path):
            raise ValueError(f"{path} must be absolute")
        if os.path.isfile(path):
            raise ValueError(f"{path} is a file, not a directory")
        self.mountpoint = path

    def _findmount(self) -> Optional[MountInfo]:
        with open('/proc/mounts') as f:
            mounts = [MountInfo(line) for line in f]
            mp = list(filter(lambda m: m.mountpoint == self.mountpoint, mounts))
            if len(mp) == 1:
                return mp[0]
            if len(mp) >= 1:
                raise ValueError(f"duplicate {mp}")
        return None

    def _umount(self) -> Optional[str]:
        """Attempt an umount. Return error if it fails"""
        if self._findmount() is None:
            return None
        try:
            cmd = ('/usr/bin/umount', self.mountpoint)
            subprocess.run(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True, check=True)
        except subprocess.CalledProcessError as cpe:
            return cpe.stderr.strip('\n')

    @staticmethod
    def _summarize(p: psutil.Process) -> str:
        """Pretty print process"""
        return f"{p.pid} {p.name()} {p.cwd()} {p.open_files()}"

    def _kill_users(self):
        """Find and kill users of mountpoint"""
        using = []
        for p in psutil.process_iter(attrs=('name', 'cmdline', 'open_files', 'cwd')):
            try:
                if p.cwd().startswith(self.mountpoint):
                    using.append(p)
                    continue
                for f in p.open_files():
                    _logger.debug(self._summarize(p))
                    if f.path.startswith(self.mountpoint):
                        using.append(p)
            except:
                _logger.error(f"{p.pid} query")
        if not using:
            return
        for u in using:
            _logger.info(f"User {self._summarize(u)}")
            u.terminate()
        time.sleep(0.5)
        remaining = list(filter(lambda u: u.is_running(), using))
        if remaining:
            _logger.info(f"Waiting for process exit")
            time.sleep(30)
            for r in remaining:
                r.kill()

    def query(self):
        """Check mountpoint and print description"""
        if not os.path.exists(self.mountpoint):
            print(f"{self.mountpoint} does not exists")
        if os.path.isdir(self.mountpoint):
            if os.path.ismount(self.mountpoint):
                if (mi := self._findmount()) is not None:
                    print(mi.description)
                else:
                    raise ValueError(f"Can't find mount info for {self.mountpoint}")
            else:
                print(f"{self.mountpoint} is an unmounted directory")
        else:
            print(f"{self.mountpoint} is not a directory")

    def umount(self):
        """Try to unmount directory gracefully"""
        if (err := self._umount()) is None:
            return
        print(err, file=sys.stderr)
        # noinspection PyProtectedMember
        os._exit(1)

    def forceunmount(self):
        """Forcefully unmount mountpoint by killing running processes"""
        if (err := self._umount()) is None:
            return
        if 'target is busy' in err:
            self._kill_users()
            if self._umount() is None:
                return
            _logger.error(f"Unable to umount {self.mountpoint}")
            os._exit(1)
        else:
            _logger.warning(err)


def main():
    logging.basicConfig()
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, epilog=_HELP)
    parser.add_argument('action', choices=('umount', 'query', 'forceunmount'))
    parser.add_argument('mountpoint')
    parser.add_argument('-l', '--loglevel', default='WARN', help="Python logging level")

    args = parser.parse_args()
    _logger.setLevel(getattr(logging, args.loglevel))
    mc = MountControl(args.mountpoint)
    getattr(mc, args.action)()


if __name__ == "__main__":
    main()
