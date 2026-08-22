#!/usr/bin/env bash
# Sign and notarize the macOS assets of a published release, then replace them.
#
# WHY THIS IS A SCRIPT AND NOT A CI STEP. The release workflow runs on GitHub's
# runners and holds no Apple credential: `grep -rn 'secrets\.' .github/workflows/`
# returns nothing. Adding the Developer ID certificate and an app-specific
# password to repository secrets is a decision with its own blast radius, so for
# now the macOS assets are signed on a trusted machine after CI has built them.
#
# That leaves a gap a person can walk into: tag a release, let CI publish four
# archives, and the macOS pair goes out unsigned while the README and the site
# say it is notarized. scripts/verify-release-assets.sh is the guard on exactly
# that, and this script is the step it checks for. Run this, then run that.
#
# The binaries this signs are the ones CI built from the tagged commit. It does
# not rebuild them locally, because then the shipped bytes would not be the ones
# the release pipeline produced.
#
# Usage:
#   APPLE_ID=... APPLE_APP_PASSWORD=... scripts/notarize-macos-release.sh v0.1.0
#
# Environment:
#   APPLE_ID             Apple ID that owns the app-specific password
#   APPLE_APP_PASSWORD   app-specific password, never a real account password
#   APPLE_TEAM_ID        defaults to F3YYBXAFJ8
#   SIGN_IDENTITY        defaults to the UCSF Developer ID Application identity
#
# Requires macOS with Xcode command line tools, plus gh and zip.

set -euo pipefail

TAG="${1:?usage: notarize-macos-release.sh <tag>}"
REPO_SLUG="${OPEN_KING_REPO:-Broccolito/open-king}"
TEAM_ID="${APPLE_TEAM_ID:-F3YYBXAFJ8}"
IDENTITY="${SIGN_IDENTITY:-Developer ID Application: University of California at San Francisco ($TEAM_ID)}"
: "${APPLE_ID:?set APPLE_ID}"
: "${APPLE_APP_PASSWORD:?set APPLE_APP_PASSWORD (an app-specific password)}"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

ASSETS=(open-king-macos-arm64 open-king-macos-x86_64)

echo "Signing and notarizing the macOS assets of $TAG"
echo

for a in "${ASSETS[@]}"; do
  echo "== $a"
  d="$WORK/$a"; mkdir -p "$d"
  gh release download "$TAG" --repo "$REPO_SLUG" --pattern "$a.zip" --dir "$WORK" --clobber
  ( cd "$d" && unzip -qo "../$a.zip" )

  # Hardened runtime and a secure timestamp are both required: Apple rejects a
  # submission without them, and the rejection arrives minutes later rather than
  # at signing time.
  codesign --force --options runtime --timestamp --sign "$IDENTITY" "$d/open-king"
  # Capture first, then match. Piping straight into `grep -q` looks equivalent and is
  # not: `grep -q` exits at the first match, which is line 4 of codesign's 14, so
  # codesign takes SIGPIPE and exits 141 -- and `set -o pipefail` above turns that into
  # the pipeline's status. The guard then fires on a signature that is perfectly good.
  # That is not hypothetical: it is why v0.1.1 shipped unsigned.
  siginfo="$(codesign -dv --verbose=2 "$d/open-king" 2>&1)"
  grep -qE 'flags=0x10000\(runtime\)' <<<"$siginfo" \
    || { echo "hardened runtime not set"; echo "$siginfo"; exit 1; }
  grep -q 'Timestamp=' <<<"$siginfo" \
    || { echo "no secure timestamp"; echo "$siginfo"; exit 1; }
  grep -q "TeamIdentifier=$TEAM_ID" <<<"$siginfo" \
    || { echo "wrong or missing TeamIdentifier"; echo "$siginfo"; exit 1; }

  # Repackage flat, so the archive holds one file with no enclosing directory
  # and no AppleDouble sidecars.
  rm -f "$WORK/$a.zip"
  ( cd "$d" && zip -q -X "$WORK/$a.zip" open-king )

  xcrun notarytool submit "$WORK/$a.zip" \
    --apple-id "$APPLE_ID" --team-id "$TEAM_ID" --password "$APPLE_APP_PASSWORD" \
    --wait --timeout 20m
  echo
done

# A bare Mach-O cannot be stapled: Apple staples tickets to app bundles, disk
# images and installer packages, not to loose executables. The ticket lives on
# Apple's side and Gatekeeper resolves it on first run, so there is deliberately
# no `xcrun stapler` call here.

echo "== checksums"
gh release download "$TAG" --repo "$REPO_SLUG" --pattern 'open-king-linux-x86_64.zip' --dir "$WORK" --clobber
gh release download "$TAG" --repo "$REPO_SLUG" --pattern 'open-king-windows-x86_64.zip' --dir "$WORK" --clobber
( cd "$WORK" && shasum -a 256 open-king-*.zip > SHA256SUMS.txt && cat SHA256SUMS.txt )
echo

echo "== uploading"
gh release upload "$TAG" --repo "$REPO_SLUG" --clobber \
  "$WORK/open-king-macos-arm64.zip" \
  "$WORK/open-king-macos-x86_64.zip" \
  "$WORK/SHA256SUMS.txt"
echo

echo "Now verify what actually shipped:"
echo "  scripts/verify-release-assets.sh $TAG"
