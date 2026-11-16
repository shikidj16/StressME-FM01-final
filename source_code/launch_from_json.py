"""Launch run_stressme.py with args from a JSON file (avoids bash quoting issues
when args contain # or other shell-significant characters)."""
import json
import os
import sys

if len(sys.argv) != 2:
    print('Usage: python launch_from_json.py <args.json>', file=sys.stderr)
    sys.exit(2)

# Use utf-8-sig to tolerate a leading BOM (PowerShell writes one by default)
import io
with io.open(sys.argv[1], encoding='utf-8-sig') as f:
    args = json.load(f)

if not isinstance(args, list):
    print('args.json must be a JSON array of strings', file=sys.stderr)
    sys.exit(2)

# Replace sys.argv with [script_path, ...args] and exec the real script.
script_path = os.path.join(os.path.dirname(__file__), 'run_stressme.py')
sys.argv = [script_path] + [str(a) for a in args]

with open(script_path) as f:
    code = f.read()

# Use exec with a proper __name__ context so the script's `if __name__ == "__main__":`
# guards don't change behavior (run_stressme.py is procedural, no guard needed)
exec(compile(code, script_path, 'exec'), {'__name__': '__main__', '__file__': script_path})
