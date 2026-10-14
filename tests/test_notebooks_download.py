"""Download of the notebooks and their example files (spearhead-vda-notebooks), without the network"""
import io
from urllib.error import HTTPError

import pytest

from spearhead.vda import __version__, notebooks


@pytest.fixture
def downloads(monkeypatch):
    """The URLs requested, each answered with its path as content"""
    urls = []

    def urlopen(url, timeout=None):
        urls.append(url)
        return io.BytesIO(url.rsplit("/", 1)[1].encode())

    monkeypatch.setattr(notebooks, "urlopen", urlopen)
    return urls


def test_copy_notebooks(downloads, tmp_path):
    copied = notebooks.copy_notebooks(tmp_path)
    # the files of the installed version
    assert downloads[0] == f"https://raw.githubusercontent.com/spearhead-he/VDA/v{__version__}/vda_multi_event.ipynb"
    assert len(downloads) == len(notebooks.FILES)
    assert set(copied.values()) == {"downloaded"}
    assert (tmp_path / "vda_custom_points.ipynb").read_text() == "vda_custom_points.ipynb"
    assert (tmp_path / "examples" / "custom_points_example.csv").read_text() == "custom_points_example.csv"


def test_existing_files_kept(downloads, tmp_path):
    (tmp_path / "vda_multi_event.ipynb").write_text("my notebook")
    copied = notebooks.copy_notebooks(tmp_path, version="1.2.3")
    assert copied[str(tmp_path / "vda_multi_event.ipynb")] == "kept"
    assert (tmp_path / "vda_multi_event.ipynb").read_text() == "my notebook"
    assert all("/v1.2.3/" in url for url in downloads) and len(downloads) == len(notebooks.FILES) - 1
    notebooks.copy_notebooks(tmp_path, overwrite=True)
    assert (tmp_path / "vda_multi_event.ipynb").read_text() == "vda_multi_event.ipynb"


@pytest.mark.parametrize("code", [404, 500])
def test_download_errors(monkeypatch, tmp_path, code):
    def urlopen(url, timeout=None):
        raise HTTPError(url, code, "error", None, None)

    monkeypatch.setattr(notebooks, "urlopen", urlopen)
    if code == 404:
        with pytest.raises(ValueError, match="version 0.0.1 not found on GitHub"):
            notebooks.copy_notebooks(tmp_path, version="0.0.1")
    else:
        with pytest.raises(HTTPError):
            notebooks.copy_notebooks(tmp_path)


def test_command(downloads, tmp_path, capsys):
    (tmp_path / "vda_multi_event.ipynb").write_text("my notebook")
    notebooks.main([str(tmp_path)])
    out = capsys.readouterr().out
    assert f"Kept the existing {tmp_path / 'vda_multi_event.ipynb'}" in out
    assert f"Downloaded {tmp_path / 'examples' / 'custom_points_example.csv'}" in out
    assert "Use --overwrite to replace the existing files." in out
    notebooks.main([str(tmp_path), "--overwrite"])
    assert "Use --overwrite" not in capsys.readouterr().out


def test_command_error(monkeypatch, tmp_path, capsys):
    def urlopen(url, timeout=None):
        raise HTTPError(url, 404, "Not Found", None, None)

    monkeypatch.setattr(notebooks, "urlopen", urlopen)
    with pytest.raises(SystemExit) as e:
        notebooks.main([str(tmp_path), "--version", "9.9.9"])
    assert e.value.code == 1
    assert "Error: vda_multi_event.ipynb of version 9.9.9 not found" in capsys.readouterr().err
