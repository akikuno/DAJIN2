from __future__ import annotations

import subprocess

import pytest

from src.DAJIN2 import gui


def test_build_completion_message_with_result_path():
    message = gui.build_completion_message("/tmp/DAJIN_Results")

    assert "analysis results are saved" in message
    assert "/tmp/DAJIN_Results" in message


def test_build_completion_message_for_batch():
    message = gui.build_completion_message("/tmp/DAJIN_Results", is_batch=True)

    assert "batch analysis results are saved" in message
    assert "/tmp/DAJIN_Results" in message


def test_parse_threads_returns_clamped_value():
    assert gui.parse_threads("0") == 1
    assert gui.parse_threads("64") == 32
    assert gui.parse_threads("4") == 4


def test_normalize_project_name_removes_unsafe_characters():
    assert gui.normalize_project_name("project name!!") == "project_name"


def test_normalize_project_name_raises_error_for_empty_name():
    with pytest.raises(ValueError, match="Project name must include"):
        gui.normalize_project_name("%%%")


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
