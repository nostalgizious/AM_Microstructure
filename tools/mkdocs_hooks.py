"""Include repository guides and resolve their links in the MkDocs source tree."""
from pathlib import Path
import posixpath
import re
from urllib.parse import unquote, urlsplit, urlunsplit

from mkdocs.structure.files import File, InclusionLevel
from build_site import ROOT, source_files


def source_uri(path):
    rel = path.relative_to(ROOT)
    if rel.parts[0] == "docs":
        return Path(*rel.parts[1:]).as_posix()
    uri = "source/" + rel.as_posix()
    if path.name == "README.md":
        uri = str(Path(uri).with_name("overview.md"))
    return uri + ".txt" if path.suffix == ".html" else uri


def on_config(config):
    # Pages handles missing routes at the site root, above this /docs subtree.
    config.theme.static_templates.discard("404.html")
    return config


def rewrite_links(content, source, uri, targets):
    def replace(match):
        url = urlsplit(match[2])
        if url.scheme or url.netloc or not url.path:
            return match[0]
        target = (source.parent / unquote(url.path)).resolve()
        if target in targets:
            target_uri = targets[target]
        elif target.is_dir() and target.is_relative_to(ROOT):
            target_uri = "source/" + target.relative_to(ROOT).as_posix() + "/index.md"
        else:
            raise ValueError(f"Missing packaged target in {source.relative_to(ROOT)}: {url.path}")
        href = posixpath.relpath(target_uri, posixpath.dirname(uri) or ".")
        return match[1] + urlunsplit(("", "", href, url.query, url.fragment)) + match[3]
    return re.sub(r'(!?\[[^\]]*\]\()([^\s)]+)(\))', replace, content)


def on_files(files, config):
    paths = source_files()
    targets = {path.resolve(): source_uri(path) for path in paths}
    for path in paths:
        uri = targets[path.resolve()]
        if path.suffix == ".md":
            content = rewrite_links(path.read_text(), path, uri, targets)
            file = File.generated(config, uri, content=content)
        elif path.is_relative_to(ROOT / "docs"):
            continue
        else:
            file = File.generated(config, uri, abs_src_path=str(path))
        if existing := files.get_file_from_path(uri):
            files.remove(existing)
        if uri.startswith("source/"):
            file.inclusion = InclusionLevel.NOT_IN_NAV
        files.append(file)

    # Browsable source directories for links to examples and coupling routines.
    children = {}
    for file in list(files):
        if not file.src_uri.startswith("source/"):
            continue
        parts = file.src_uri.split("/")
        for depth in range(1, len(parts)):
            parent = "/".join(parts[:depth])
            name = parts[depth]
            target = name if depth == len(parts)-1 else name + "/index.md"
            children.setdefault(parent, {})[name] = target
    for directory, entries in children.items():
        uri = directory + "/index.md"
        if files.get_file_from_path(uri):
            continue
        content = "# " + directory.removeprefix("source/") + "\n\n"
        content += "\n".join(f"- [{name}]({target})" for name, target in sorted(entries.items())) + "\n"
        file = File.generated(config, uri, content=content)
        file.inclusion = InclusionLevel.NOT_IN_NAV
        files.append(file)
    return files
