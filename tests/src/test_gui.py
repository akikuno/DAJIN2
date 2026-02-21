from __future__ import annotations

import subprocess

from src.DAJIN2 import gui


def test_open_browser_uses_webbrowser_outside_wsl(monkeypatch):
    opened_urls = []

    monkeypatch.setattr(gui, "is_wsl_environment", lambda: False)
    monkeypatch.setattr(gui.webbrowser, "open_new", lambda url: opened_urls.append(url))

    gui.open_browser(12345)

    assert opened_urls == ["http://localhost:12345/"]


def test_open_browser_does_not_use_linux_webbrowser_on_wsl(monkeypatch):
    opened_urls = []
    windows_open_calls = []

    monkeypatch.setattr(gui, "is_wsl_environment", lambda: True)
    monkeypatch.setattr(gui, "open_url_in_windows_default_browser", lambda url: windows_open_calls.append(url) or True)
    monkeypatch.setattr(gui.webbrowser, "open_new", lambda url: opened_urls.append(url))

    gui.open_browser(23456)

    assert windows_open_calls == ["http://localhost:23456/"]
    assert opened_urls == []


def test_open_url_in_windows_default_browser_returns_false_on_failure(monkeypatch):
    def raise_file_not_found(*args, **kwargs):
        raise FileNotFoundError

    monkeypatch.setattr(gui.subprocess, "run", raise_file_not_found)

    assert gui.open_url_in_windows_default_browser("http://localhost:34567/") is False


def test_open_url_in_windows_default_browser_uses_cmd_first(monkeypatch):
    executed_commands = []

    def mock_run(command, **kwargs):
        executed_commands.append(command)
        return None

    monkeypatch.setattr(gui.subprocess, "run", mock_run)

    result = gui.open_url_in_windows_default_browser("http://localhost:45678/")

    assert result is True
    assert executed_commands == [["cmd.exe", "/c", "start", "", "http://localhost:45678/"]]


def test_open_url_in_windows_default_browser_falls_back_to_powershell(monkeypatch):
    executed_commands = []

    def mock_run(command, **kwargs):
        executed_commands.append(command)
        if command[0] == "cmd.exe":
            raise subprocess.CalledProcessError(returncode=1, cmd=command)
        return None

    monkeypatch.setattr(gui.subprocess, "run", mock_run)

    result = gui.open_url_in_windows_default_browser("http://localhost:56789/")

    assert result is True
    assert executed_commands == [
        ["cmd.exe", "/c", "start", "", "http://localhost:56789/"],
        ["powershell.exe", "-NoProfile", "-Command", 'Start-Process "http://localhost:56789/"'],
    ]
