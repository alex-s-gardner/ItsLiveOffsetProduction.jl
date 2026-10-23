#!/usr/bin/env bash
# Make every library in a bundle find its dependencies inside the bundle.
#
#   app/link_libraries.sh <bundle>
#
# A bundled library finds a dependency in another artifact by the dependency's SONAME (Linux) or install name
# (macOS), which the Julia process registers when it loads the dependency first. A library whose dependency
# is requested under a different name searches only its own directory (RUNPATH `$ORIGIN`, LC_RPATH
# `@loader_path`) and then the system: libproj and libgdal ask for `libsqlite3.so.0` / `libsqlite3.0.dylib`,
# while SQLite_jll ≥ 3.51 names its library `libsqlite3.so.3.<version>` / `libsqlite3.dylib`. Such a library
# fails to load on a host without that system library, or loads the host's copy. For each such dependency
# this links the file of that name from the artifact that ships it into the requesting library's directory,
# and fails if no artifact ships it. C library and loader dependencies on Linux come from the host.
set -euo pipefail
out=$1
arts="$out/share/julia/artifacts"
case "$(uname -s)" in
    Linux)
        libs() { find "$out" -type f -name '*.so*'; }
        names() { readelf -d "$1" 2>/dev/null | sed -n 's/.*(SONAME).*\[\(.*\)\]/\1/p'; }
        deps() { readelf -d "$1" 2>/dev/null | sed -n 's/.*(NEEDED).*\[\(.*\)\]/\1/p'; }
        host='^(ld-linux.*|lib(c|dl|m|pthread|resolv|rt|util))\.so\.[0-9]+$' ;;
    Darwin)
        libs() { find "$out" -type f -name '*.dylib'; }
        names() { otool -D "$1" 2>/dev/null | sed -n '2s#.*/##p'; }
        deps() { otool -L "$1" 2>/dev/null | sed -n 's#^[[:space:]]*@rpath/\([^ ]*\) .*#\1#p'; }
        host='^$' ;;
    *) echo "unsupported OS $(uname -s)" >&2; exit 1 ;;
esac
known=$(libs | while IFS= read -r f; do names "$f"; done | sort -u)
missing=
while IFS= read -r lib; do
    dir=$(dirname "$lib")
    for need in $(deps "$lib"); do
        grep -qxF "$need" <<< "$known" || [[ $need =~ $host ]] || [ -e "$dir/$need" ] && continue
        # An artifact's own links stay inside its directory; one leading out of it was made here.
        provider=$(find "$arts" -mindepth 3 -maxdepth 3 -path '*/lib/*' -name "$need" ! -lname '../*' | head -n 1)
        if [ -n "$provider" ] && [[ $dir == "$arts"/* ]]; then
            # Relative, so the bundle stays relocatable.
            up=$(sed 's#[^/]*#..#g' <<< "${dir#"$arts"/}")
            ln -s "$up/${provider#"$arts"/}" "$dir/$need"
            echo "linked ${dir#"$out"/}/$need -> $up/${provider#"$arts"/}"
        else
            missing+="  ${lib#"$out"/} needs $need"$'\n'
        fi
    done
done < <(libs)
[ -z "$missing" ] || { printf 'no bundled library provides:\n%s' "$missing" >&2; exit 1; }
