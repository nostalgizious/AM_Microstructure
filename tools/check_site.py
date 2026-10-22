#!/usr/bin/env python3
"""Check every built HTML local link/asset and HTML fragment without network access."""
from html.parser import HTMLParser
from pathlib import Path
import hashlib
import json
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1] / "website/dist"


class Page(HTMLParser):
    def __init__(self, path):
        super().__init__()
        self.links, self.ids = [], set()
        self.feed(path.read_text())

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            self.ids.add(attrs["id"])
        for key in ("src", "href"):
            if attrs.get(key):
                self.links.append(attrs[key])


def main():
    pages = {path: Page(path) for path in ROOT.rglob("*.html")}
    if not pages:
        raise SystemExit("No built HTML found; run tools/build_site.py first.")
    errors, count = [], 0
    for path, page in pages.items():
        for link in page.links:
            url = urlsplit(link)
            if url.scheme or url.netloc:
                continue
            target = (path.parent / unquote(url.path)).resolve() if url.path else path
            if target.is_dir():
                target /= "index.html"
            count += 1
            if not target.is_relative_to(ROOT) or not target.is_file():
                errors.append(f"{path.relative_to(ROOT)}: missing/outside target {link}")
            elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
                errors.append(f"{path.relative_to(ROOT)}: missing anchor {link}")
    # Playback references these images dynamically, outside static HTML links.
    timeline = json.loads((ROOT / "assets/demo/evolution.json").read_text())
    for frame in timeline["frames"]:
        image = ROOT / "assets/demo/evolution" / f"frame-{frame['step']:03d}.png"
        count += 1
        if not image.is_file():
            errors.append(f"Missing simulation frame: {image.name}")
    manifest = json.loads((ROOT / "assets/demo/manifest.json").read_text())
    for entry in manifest["files"]:
        asset = ROOT / "assets/demo" / entry["file"]
        if not asset.is_file() or hashlib.sha256(asset.read_bytes()).hexdigest() != entry["sha256"]:
            errors.append(f"Figure hash mismatch: {entry['file']}")
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"Checked {len(pages)} HTML pages and {count} local links/assets; all resolve.")


if __name__ == "__main__":
    main()
