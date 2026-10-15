"""The fully resolved condition is copied into each independent run snapshot."""
import json
from pathlib import Path

CONFIG = json.loads(Path(__file__).with_name('condition.json').read_text())
PROTOCOL = CONFIG['channel_protocol']
