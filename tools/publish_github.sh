#!/usr/bin/env bash
# Create the GitHub repository and publish the combined research site and docs.
set -euo pipefail

project_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"
git rev-parse --verify HEAD >/dev/null
if [[ -n "$(git status --porcelain)" ]]; then
    echo "Commit the project changes before publishing." >&2
    exit 1
fi
if [[ "$(git branch --show-current)" != main ]]; then
    echo "Switch to the main branch before publishing." >&2
    exit 1
fi

gh auth status
if git remote get-url origin >/dev/null 2>&1; then
    repository="$(gh repo view --json nameWithOwner --jq .nameWithOwner)"
    if [[ "${repository##*/}" != AM_Microstructure ]]; then
        echo "The existing origin does not point to AM_Microstructure." >&2
        exit 1
    fi
else
    owner="$(gh api user --jq .login)"
    repository="$owner/AM_Microstructure"
    gh repo create "$repository" --public --source=. --remote=origin \
        --description "Cellular automaton modelling of grain growth in PBF-LB/M — Ferdinand Schnell's thesis work"
fi

git push --set-upstream origin main
if gh api "repos/$repository/pages" --silent 2>/dev/null; then
    gh api --method PUT "repos/$repository/pages" -f build_type=workflow --silent
else
    gh api --method POST "repos/$repository/pages" -f build_type=workflow --silent
fi
pages_url="$(gh api "repos/$repository/pages" --jq .html_url)"
gh repo edit "$repository" --homepage "$pages_url"
# Run after enabling Pages, including when the first push preceded that setting.
gh workflow run pages.yml --repo "$repository" --ref main
printf 'Pages deployment started. Website: %s\n' "$pages_url"
printf 'Deployment status: https://github.com/%s/actions/workflows/pages.yml\n' "$repository"
