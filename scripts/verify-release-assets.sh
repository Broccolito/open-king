#!/usr/bin/env bash
# Verify a published release's assets are what the documentation promises.
#
# The release workflow builds and uploads unsigned archives. The macOS pair is
# then signed and notarized by hand (scripts/notarize-macos-release.sh), which
# means the promise the README and the site make about macOS is kept by a human
# step rather than by CI. This script is the guard on that step: run it after
# publishing and it fails loudly if the macOS assets went out unsigned.
#
# Usage:  scripts/verify-release-assets.sh [tag]     (default: the latest tag)
#
# Requires: gh, unzip, shasum. codesign checks run on macOS only; elsewhere the
# signature assertions are skipped with a notice rather than passing silently.

set -euo pipefail

# gh needs a repository context, and every command below runs against a temp
# directory that is not a checkout, so name the repo explicitly rather than
# relying on the working directory.
REPO_SLUG="${OPEN_KING_REPO:-Broccolito/open-king}"
TAG="${1:-$(gh release view --repo "$REPO_SLUG" --json tagName --jq .tagName)}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

EXPECTED=(
  open-king-macos-arm64.zip
  open-king-macos-x86_64.zip
  open-king-linux-x86_64.zip
  open-king-windows-x86_64.zip
  SHA256SUMS.txt
)
MAC_ASSETS=(open-king-macos-arm64.zip open-king-macos-x86_64.zip)
TEAM_ID="F3YYBXAFJ8"

fail() { printf '  FAIL  %s\n' "$*"; FAILED=$((FAILED + 1)); }
pass() { printf '  ok    %s\n' "$*"; }
FAILED=0

echo "Verifying release $TAG"
echo

echo "Asset set"
ACTUAL="$(gh release view "$TAG" --repo "$REPO_SLUG" --json assets --jq '.assets[].name' | sort)"
for a in "${EXPECTED[@]}"; do
  grep -qx "$a" <<<"$ACTUAL" && pass "$a present" || fail "$a missing"
done
# Anything extra is usually a stale asset from a previous packaging scheme.
while read -r a; do
  [ -z "$a" ] && continue
  printf '%s\n' "${EXPECTED[@]}" | grep -qx "$a" || fail "unexpected asset: $a"
done <<<"$ACTUAL"
echo

echo "Download and checksums"
if gh release download "$TAG" --repo "$REPO_SLUG" --dir "$WORK" --clobber >/dev/null 2>&1; then
  pass "downloaded every asset"
else
  fail "could not download the assets for $TAG"
fi
if ( cd "$WORK" && shasum -a 256 -c SHA256SUMS.txt >/dev/null 2>&1 ); then
  pass "every checksum in SHA256SUMS.txt matches"
else
  fail "checksum mismatch, or SHA256SUMS.txt does not cover every archive"
fi
echo

echo "Archive contents"
# The documented promise is one executable per archive and nothing else, so a
# reader can unzip and run without tidying up afterwards.
for z in open-king-macos-arm64.zip open-king-macos-x86_64.zip open-king-linux-x86_64.zip; do
  n="$(unzip -Z1 "$WORK/$z" | wc -l | tr -d ' ')"
  entry="$(unzip -Z1 "$WORK/$z" | head -1)"
  [ "$n" = "1" ] && [ "$entry" = "open-king" ] \
    && pass "$z holds exactly open-king" \
    || fail "$z holds $n entr(y/ies), first is '$entry'"
done
n="$(unzip -Z1 "$WORK/open-king-windows-x86_64.zip" | wc -l | tr -d ' ')"
entry="$(unzip -Z1 "$WORK/open-king-windows-x86_64.zip" | head -1)"
[ "$n" = "1" ] && [ "$entry" = "open-king.exe" ] \
  && pass "open-king-windows-x86_64.zip holds exactly open-king.exe" \
  || fail "windows zip holds $n entr(y/ies), first is '$entry'"
echo

echo "macOS signing and notarization"
if ! command -v codesign >/dev/null 2>&1; then
  echo "  note  codesign is unavailable on this host, so the signature assertions"
  echo "        were NOT run. Re-run this script on macOS before announcing."
else
  for z in "${MAC_ASSETS[@]}"; do
    d="$WORK/x_${z%.zip}"; mkdir -p "$d"
    ( cd "$d" && unzip -qo "$WORK/$z" )
    info="$(codesign -dv --verbose=2 "$d/open-king" 2>&1 || true)"
    grep -q "Authority=Developer ID Application" <<<"$info" \
      && pass "$z carries a Developer ID signature" \
      || fail "$z is NOT Developer ID signed (this is the hand step being skipped)"
    grep -q "flags=0x10000(runtime)" <<<"$info" \
      && pass "$z has the hardened runtime" \
      || fail "$z lacks the hardened runtime, so Apple would refuse to notarize it"
    grep -q "TeamIdentifier=$TEAM_ID" <<<"$info" \
      && pass "$z is signed by team $TEAM_ID" \
      || fail "$z has the wrong or no TeamIdentifier"
    grep -q "^Timestamp=" <<<"$info" \
      && pass "$z has a secure timestamp" \
      || fail "$z has no secure timestamp"
  done
fi
echo

if [ "$FAILED" -gt 0 ]; then
  echo "$FAILED check(s) failed."
  exit 1
fi
echo "All checks passed."
