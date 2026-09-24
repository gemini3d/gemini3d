#!/usr/bin/env bash
set -euo pipefail
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
root="${here}/../build/local"
arguments=("$@")
for ((i=0; i<${#arguments[@]}; i++)); do
  case "${arguments[i]}" in
    --root)
      ((i+=1))
      [[ $i -lt ${#arguments[@]} && -n "${arguments[i]}" ]] || { echo "--root requires a directory" >&2; exit 1; }
      root="${arguments[i]}"
      ;;
    --root=*) root="${arguments[i]#--root=}" ;;
  esac
done
case "$root" in
  "~") root="$HOME" ;;
  "~/"*) root="$HOME/${root:2}" ;;
  "~"*) echo "Expand the --root home directory before invoking the installer." >&2; exit 1 ;;
esac
[[ -n "$root" ]] || { echo "--root requires a directory" >&2; exit 1; }

for argument in "$@"; do
  if [[ "$argument" == "--system-deps" ]]; then
    # Package managers can fail after changing native libraries.
    rm -f -- "${root}/environment.json"
    bash "${here}/install-system-deps.sh"
    break
  fi
done

python="${PYTHON:-python3}"
if [[ "$(uname -s)" == Darwin ]] && command -v brew >/dev/null; then
  brew_prefix="$(brew --prefix)"
  if [[ -z "${PYTHON:-}" && -x "${brew_prefix}/opt/python@3.12/bin/python3.12" ]]; then
    python="${brew_prefix}/opt/python@3.12/bin/python3.12"
  fi
  export PATH="${brew_prefix}/bin:${PATH}"
  export CMAKE_PREFIX_PATH="${brew_prefix}${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}"
fi
exec "$python" "${here}/local_environment.py" install "$@"
