#!/usr/bin/env python3
"""Static checks for the documentation site in `site/`.

The site is hand-written HTML with no build step, so nothing else would catch a
dangling link, a class with no rule behind it, or a page that still carries its
placeholder callout. This script is that safety net, and it runs in CI.

It also enforces two prose rules the site is held to: no em dashes, and none of
the register words that make technical writing read as machine-generated. Those
are style rules, not correctness rules, but they are the kind of thing that
erodes one commit at a time unless something checks.

Exit status is 0 when every check passes and 1 otherwise, with every failure
printed rather than only the first, so one run tells you everything to fix.

Standard library only, to match the rest of the repository's tooling.
"""

from __future__ import annotations

import html.parser
import re
import sys
from pathlib import Path

SITE = Path(__file__).resolve().parent.parent / "site"

# U+2014 EM DASH and U+2013 EN DASH. The en dash is included because it is the
# usual substitute once the em dash is banned, and a spaced en dash reads the
# same way. A hyphen in a compound word is fine and is not matched here.
DASHES = {"—": "em dash", "–": "en dash"}

# Register words. Matched on word boundaries and case-insensitively. Kept to
# terms that have a plain alternative in every context this site uses, so a hit
# is always fixable rather than something to argue about.
BANNED_WORDS = [
    "delve", "leverage", "robust", "seamless", "powerful", "comprehensive",
    "unlock", "elevate", "harness", "empower", "streamline", "effortlessly",
    "cutting-edge", "game-changer", "utilize", "myriad", "plethora",
]
# Terms that contain a banned word but are proper names, not register. The
# estimator really is called KING-robust, in the 2010 paper and in 46 places in
# this repository's own documentation. Banning the word outright made a page
# drop the estimator's name, which trades an accuracy loss for a style win.
# These are removed from the text before the banned-word scan, so bare "robust"
# as an adjective is still caught.
ALLOWED_COMPOUNDS = ["KING-robust"]

BANNED_PHRASES = [
    "at its core", "it's worth noting", "it is worth noting", "in today's world",
    "dive in", "a wide range of", "when it comes to", "the world of",
    "needless to say", "rest assured",
]

# Markers that a page is still holding scaffold content. The class is the
# reliable one: the scaffold wraps every placeholder note in `note--draft`, so a
# page that still carries it has not been written yet, whatever its prose says.
# Matching on the visible word "Draft" alone would miss a page whose placeholder
# copy happens not to use it, and would fire on a page that legitimately
# discusses a draft.
PLACEHOLDER_CLASS = "note--draft"
PLACEHOLDER_MARKERS = ["Lorem ipsum", "placeholder text", "TODO:", "FIXME:"]

VOID_ELEMENTS = {
    "area", "base", "br", "col", "embed", "hr", "img", "input",
    "link", "meta", "source", "track", "wbr",
}


class PageParser(html.parser.HTMLParser):
    """Collects ids, classes, link targets and an element stack for balance."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.ids: set[str] = set()
        self.classes: set[str] = set()
        self.links: list[tuple[str, str]] = []  # (attribute, value)
        self.stack: list[tuple[str, int]] = []
        self.unbalanced: list[str] = []

    def handle_starttag(self, tag, attrs):
        d = dict(attrs)
        if d.get("id"):
            self.ids.add(d["id"])
        for cls in (d.get("class") or "").split():
            self.classes.add(cls)
        for attr in ("href", "src"):
            if d.get(attr):
                self.links.append((attr, d[attr]))
        if tag not in VOID_ELEMENTS:
            self.stack.append((tag, self.getpos()[0]))

    def handle_endtag(self, tag):
        if tag in VOID_ELEMENTS:
            return
        if not self.stack:
            self.unbalanced.append(f"stray </{tag}> at line {self.getpos()[0]}")
            return
        open_tag, line = self.stack.pop()
        if open_tag != tag:
            self.unbalanced.append(
                f"<{open_tag}> opened at line {line} closed by </{tag}> "
                f"at line {self.getpos()[0]}"
            )


def visible_text(source: str) -> str:
    """The page text with <script>, <style> and tags removed.

    Prose rules apply to what a reader sees. A banned word inside a CSS class
    name or a code sample is not a style problem, and matching it would make the
    check something people route around.
    """
    source = re.sub(r"<script\b.*?</script>", " ", source, flags=re.S | re.I)
    source = re.sub(r"<style\b.*?</style>", " ", source, flags=re.S | re.I)
    source = re.sub(r"<code\b.*?</code>", " ", source, flags=re.S | re.I)
    source = re.sub(r"<pre\b.*?</pre>", " ", source, flags=re.S | re.I)
    return re.sub(r"<[^>]+>", " ", source)


def line_of(source: str, index: int) -> int:
    return source.count("\n", 0, index) + 1


def check_page(path: Path, css_classes: set[str], failures: list[str]) -> None:
    rel = path.relative_to(SITE.parent)
    source = path.read_text(encoding="utf-8")
    text = visible_text(source)

    for char, name in DASHES.items():
        for m in re.finditer(re.escape(char), source):
            failures.append(f"{rel}:{line_of(source, m.start())}: {name} is not allowed")

    scan = text
    for allowed in ALLOWED_COMPOUNDS:
        scan = re.sub(re.escape(allowed), " ", scan, flags=re.I)

    for word in BANNED_WORDS:
        for m in re.finditer(rf"\b{re.escape(word)}\b", scan, re.I):
            failures.append(f"{rel}: banned word {word!r} (near {scan[m.start():m.start()+60].strip()!r})")

    for phrase in BANNED_PHRASES:
        for m in re.finditer(re.escape(phrase), text, re.I):
            failures.append(f"{rel}: banned phrase {phrase!r}")

    if PLACEHOLDER_CLASS in source:
        failures.append(f"{rel}: still carries a {PLACEHOLDER_CLASS} note; the page is scaffold, not content")
    for marker in PLACEHOLDER_MARKERS:
        if marker in source:
            failures.append(f"{rel}: still contains placeholder marker {marker!r}")

    parser = PageParser()
    parser.feed(source)
    for problem in parser.unbalanced:
        failures.append(f"{rel}: unbalanced markup: {problem}")
    for tag, line in parser.stack:
        failures.append(f"{rel}:{line}: <{tag}> is never closed")

    for attr, value in parser.links:
        if value.startswith(("http://", "https://", "mailto:", "data:")):
            continue
        if value.startswith("#"):
            if value[1:] and value[1:] not in parser.ids:
                failures.append(f"{rel}: {attr}={value!r} points at an id that does not exist")
            continue
        target_path, _, fragment = value.partition("#")
        target = (path.parent / target_path).resolve()
        if not target.exists():
            failures.append(f"{rel}: {attr}={value!r} points at a missing file")
        elif fragment and target.suffix == ".html":
            other = PageParser()
            other.feed(target.read_text(encoding="utf-8"))
            if fragment not in other.ids:
                failures.append(f"{rel}: {attr}={value!r} points at an id missing from {target_path}")

    for cls in sorted(parser.classes - css_classes):
        failures.append(f"{rel}: class {cls!r} is used but has no rule in style.css")


def main() -> int:
    if not SITE.is_dir():
        print(f"no site directory at {SITE}", file=sys.stderr)
        return 1

    css_path = SITE / "assets" / "style.css"
    if not css_path.exists():
        print(f"missing {css_path}", file=sys.stderr)
        return 1

    # Every class named anywhere in a selector. Deliberately loose: the point is
    # to catch a class that exists nowhere in the stylesheet, not to model CSS
    # specificity.
    css_classes = set(re.findall(r"\.([A-Za-z_][-\w]*)", css_path.read_text(encoding="utf-8")))

    pages = sorted(SITE.glob("*.html"))
    if not pages:
        print(f"no pages found in {SITE}", file=sys.stderr)
        return 1

    failures: list[str] = []
    for page in pages:
        check_page(page, css_classes, failures)

    if failures:
        for f in failures:
            print(f)
        print(f"\n{len(failures)} problem(s) across {len(pages)} page(s)")
        return 1

    print(f"site checks pass: {len(pages)} pages, {len(css_classes)} classes in style.css")
    return 0


if __name__ == "__main__":
    sys.exit(main())
