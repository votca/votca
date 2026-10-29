#!/usr/bin/env python3
#
# Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""Checks for the 'key value unit' result files of grw, kmctof, kmcmultiple
and kmclifetime (integration tests).

  check_results.py reference RESULT REFERENCE [--rtol R]
      every key of REFERENCE must be in RESULT with the same value
      (relative tolerance; strings exactly)
  check_results.py velocity VELOCITIES RESULT [--rtol R]
      the velocities per segment written by vaverage (nm/s) add up to the
      exact single-carrier drift velocity of the kmcmultiple RESULT file
      within R relative to its length (the occupations come from KMC)
  check_results.py agree RESULT KEY ERRORKEY EXACTFILE EXACTKEY [--sigma S]
                         [--max-rel M]
      |RESULT[KEY] - EXACTFILE[EXACTKEY]| < S * RESULT[ERRORKEY], and the
      error is at most M times the exact value (so that the test means
      something)
"""
import argparse
import sys


def read(name):
    values = {}
    with open(name) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            try:
                values[parts[0]] = float(parts[1])
            except ValueError:
                values[parts[0]] = parts[1]
    return values


def reference(args):
    result = read(args.result)
    ok = True
    for key, ref in read(args.reference).items():
        if key not in result:
            print(f"missing key {key}")
            ok = False
            continue
        value = result[key]
        if isinstance(ref, str) or isinstance(value, str):
            good = value == ref
        else:
            good = abs(value - ref) <= args.rtol * max(abs(ref), abs(value)) \
                or abs(value - ref) <= args.atol
        print(f"{key}: {value} reference {ref} {'ok' if good else 'DIFFERS'}")
        ok = ok and good
    return ok


def agree(args):
    result = read(args.result)
    exact = read(args.exactfile)[args.exactkey]
    value = result[args.key]
    error = result[args.errorkey]
    z = (value - exact) / error if error > 0 else float('inf')
    print(f"{args.key} = {value} +- {error}, {args.exactkey} = {exact}, "
          f"{z:.2f} standard errors")
    ok = error > 0 and abs(z) < args.sigma
    if error > args.max_rel * abs(exact):
        print(f"standard error larger than {args.max_rel} of the value")
        ok = False
    return ok


def velocity(args):
    total = [0.0, 0.0, 0.0]
    with open(args.velocities) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            for a in range(3):
                total[a] += float(parts[5 + a]) * 1e-9  # m/s
    result = read(args.result)
    exact = [result['exact_velocity_' + a] for a in 'xyz']
    diff = sum((t - e)**2 for t, e in zip(total, exact))**0.5
    norm = sum(e**2 for e in exact)**0.5
    print(f"sum of the segment velocities {total} m/s, exact {exact} m/s, "
          f"relative difference {diff / norm:.3f}")
    return diff < args.rtol * norm


parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
sub = parser.add_subparsers(dest='command', required=True)
p = sub.add_parser('reference')
p.add_argument('result')
p.add_argument('reference')
p.add_argument('--rtol', type=float, default=1e-6)
p.add_argument('--atol', type=float, default=0.0)
p = sub.add_parser('velocity')
p.add_argument('velocities')
p.add_argument('result')
p.add_argument('--rtol', type=float, default=0.05)
p = sub.add_parser('agree')
p.add_argument('result')
p.add_argument('key')
p.add_argument('errorkey')
p.add_argument('exactfile')
p.add_argument('exactkey')
p.add_argument('--sigma', type=float, default=4.5)
p.add_argument('--max-rel', type=float, default=0.05)
args = parser.parse_args()
commands = {'reference': reference, 'agree': agree, 'velocity': velocity}
ok = commands[args.command](args)
sys.exit(0 if ok else 1)
