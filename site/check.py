"""Check the internal links of the generated site under a Pages project prefix."""
import json
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import quote, unquote, urljoin, urlsplit

ROOT = Path(__file__).resolve().parent.parent
DIST = ROOT / "dist"
BASE = "https://preview.invalid/project/"


class Page(HTMLParser):
    def __init__(self, file):
        super().__init__()
        self.ids = set()
        self.links = []
        self.feed(file.read_text(encoding="utf-8"))
        self.close()

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if values.get("id"):
            self.ids.add(values["id"])
        for key in ("href", "src"):
            if values.get(key):
                self.links.append(values[key])


def check():
    files = [DIST / "index.html", DIST / "guide/index.html"]
    for topic in json.loads((ROOT / "site/topics.json").read_text(encoding="utf-8")):
        directory = DIST / topic["id"]
        files.append(directory / "index.html")
        files.extend((directory / "articles").rglob("*.html"))
    pages = {file: Page(file) for file in files}
    problems = []
    checked = 0
    for file, page in pages.items():
        name = file.relative_to(DIST).as_posix()
        for link in page.links:
            url = urlsplit(urljoin(BASE + quote(name), link))
            if url.scheme not in ("http", "https") or url.netloc != urlsplit(BASE).netloc:
                continue
            prefix = urlsplit(BASE).path
            if not url.path.startswith(prefix):
                problems.append(f"{name}: {link} escapes the Pages project prefix")
                continue
            target = DIST / unquote(url.path[len(prefix):])
            if target.is_dir():
                target /= "index.html"
            if not target.is_file():
                problems.append(f"{name}: missing target for {link}")
            elif url.fragment and target in pages and unquote(url.fragment) not in pages[target].ids:
                problems.append(f"{name}: missing anchor for {link}")
            checked += 1
    if problems:
        raise SystemExit("\n".join(problems))
    print(f"Verified {checked} internal links on {len(pages)} pages.")


if __name__ == "__main__":
    check()
