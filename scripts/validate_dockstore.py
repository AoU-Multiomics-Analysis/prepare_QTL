#!/usr/bin/env python3
"""Validate the repository's Dockstore descriptor registration."""

from pathlib import Path
import re
import sys

root = Path(__file__).resolve().parents[1]
dockstore = root / ".dockstore.yml"
text = dockstore.read_text(encoding="utf-8")
paths = re.findall(r"^\s*primaryDescriptorPath:\s*/(.+\.wdl)\s*$", text, re.MULTILINE)
names = re.findall(r"^\s*name:\s*(\S+)\s*$", text, re.MULTILINE)

if not paths:
    sys.exit("No WDL primaryDescriptorPath entries found in .dockstore.yml")
if len(paths) != len(names):
    sys.exit("Each Dockstore workflow must have exactly one name and primaryDescriptorPath")
if len(paths) != len(set(paths)):
    sys.exit("Duplicate Dockstore primaryDescriptorPath entries found")
if len(names) != len(set(names)):
    sys.exit("Duplicate Dockstore workflow names found")

missing = [path for path in paths if not (root / path).is_file()]
if missing:
    sys.exit("Dockstore descriptor file(s) missing: " + ", ".join(missing))

print("Validated", len(paths), "Dockstore WDL descriptor(s)")
print("\n".join(paths))
