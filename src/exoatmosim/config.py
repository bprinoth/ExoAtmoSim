"""Backward-compatible reader/writer for the existing key/value config files.

Unmodified lines are preserved byte-for-byte (apart from newline normalization
when constructing a Config from a mapping). Keys are case-sensitive.
"""
from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any

_LINE = re.compile(r"^(?P<leading>\s*)(?P<key>\S+)(?P<sep>\s+)(?P<value>.*)$")


@dataclass
class _Line:
    raw: str
    key: str | None = None
    value: str | None = None
    leading: str = ""
    sep: str = "\t"


def _serialize(value: Any) -> str:
    if isinstance(value, bool):
        return "True" if value else "False"
    if value is None:
        return "nan"
    return str(value)


class Config(Mapping[str, str]):
    """A case-sensitive, ordered mapping with typed access and lossless round trips.

    ``config['P']`` returns the raw string, while ``config.get_float('P')``
    and ``config.get_bool('air')`` explicitly convert values. Original lines,
    including comments, blank lines and extra keys, remain intact on save.
    """

    def __init__(self, lines: list[_Line], *, trailing_newline: bool = True):
        self._lines = lines
        self._trailing_newline = trailing_newline
        self._index: dict[str, int] = {}
        for i, line in enumerate(lines):
            if line.key is not None:
                if line.key in self._index:
                    raise ValueError(f"Duplicate config key {line.key!r} on line {i + 1}")
                self._index[line.key] = i

    @classmethod
    def from_text(cls, text: str) -> Config:
        lines = []
        for number, raw in enumerate(text.splitlines(), 1):
            stripped = raw.strip()
            if not stripped or stripped.startswith(("#", "!")):
                lines.append(_Line(raw))
                continue
            match = _LINE.match(raw)
            if not match or not match.group("value").strip():
                raise ValueError(f"Invalid config entry on line {number}: {raw!r}")
            lines.append(_Line(raw, match.group("key"), match.group("value").strip(),
                               match.group("leading"), match.group("sep")))
        return cls(lines, trailing_newline=text.endswith(("\n", "\r")))

    @classmethod
    def from_file(cls, path: str | Path) -> Config:
        # newline='' avoids platform newline translation for exact round trips.
        with open(path, "r", encoding="utf-8", newline="") as handle:
            text = handle.read()
        return cls.from_text(text)

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> Config:
        return cls.from_text("".join(f"{key}\t{_serialize(value)}\n" for key, value in values.items()))

    def __getitem__(self, key: str) -> str:
        if key not in self._index:
            raise KeyError(key)
        value = self._lines[self._index[key]].value
        assert value is not None
        return value

    def __iter__(self) -> Iterator[str]:
        return iter(self._index)

    def __len__(self) -> int:
        return len(self._index)

    def get_float(self, key: str, default: float | None = None) -> float:
        if key not in self:
            if default is not None:
                return default
            raise KeyError(key)
        try:
            return float(self[key])
        except ValueError as exc:
            raise ValueError(f"Expected a number for {key!r}, got {self[key]!r}") from exc

    def get_int(self, key: str, default: int | None = None) -> int:
        if key not in self:
            if default is not None:
                return default
            raise KeyError(key)
        try:
            return int(self[key])
        except ValueError as exc:
            raise ValueError(f"Expected an integer for {key!r}, got {self[key]!r}") from exc

    def get_bool(self, key: str, default: bool | None = None) -> bool:
        if key not in self:
            if default is not None:
                return default
            raise KeyError(key)
        value = self[key].strip().lower()
        if value in ("true", "yes", "1"):
            return True
        if value in ("false", "no", "0"):
            return False
        raise ValueError(f"Expected a boolean for {key!r}, got {self[key]!r}")

    def set(self, key: str, value: Any) -> None:
        """Update a key without changing other lines, or append a new key."""
        if not key or any(char.isspace() for char in key) or key.startswith(("#", "!")):
            raise ValueError(f"Invalid config key: {key!r}")
        rendered = _serialize(value)
        if "\n" in rendered or "\r" in rendered:
            raise ValueError("Config values must occupy one line")
        if key in self._index:
            line = self._lines[self._index[key]]
            line.raw = f"{line.leading}{key}{line.sep}{rendered}"
            line.value = rendered
        else:
            self._index[key] = len(self._lines)
            self._lines.append(_Line(f"{key}\t{rendered}", key, rendered))
            self._trailing_newline = True

    def to_text(self) -> str:
        text = "\n".join(line.raw for line in self._lines)
        return text + ("\n" if self._trailing_newline and self._lines else "")

    def write(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8", newline="") as handle:
            handle.write(self.to_text())
