"""Read raw chromatogram data from Agilent/ChemStation output.

The main entry point is ``read_chemstation_dataframe(path)``. Pass either a
ChemStation ``.D`` directory or an individual ``.MS`` file. The return value is
a pandas DataFrame with ``time`` in the first column and one signal column
per detected m/z channel.

This module ports the ChemStation binary offsets and compression routines from
the MATLAB ImportAgilent code supplied with chromatoPy_App.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import argparse
import csv
import json
import math
import struct
from typing import TYPE_CHECKING, BinaryIO, Iterable

if TYPE_CHECKING:
    import pandas as pd


@dataclass
class ChannelData:
    """Raw data for one ChemStation channel."""

    source: str
    kind: str
    time: list[float]
    signal: list[float]
    units: str = ""
    metadata: dict[str, object] = field(default_factory=dict)
    mz: list[float] | None = None
    xic: list[list[float]] | None = None


def read_chemstation(path: str | Path, *, include_xic: bool = False, precision: int = 3) -> list[ChannelData]:
    """Read raw channels from a ChemStation ``.D`` folder or data file.

    Parameters
    ----------
    path:
        ChemStation ``.D`` directory, ``.MS`` file, or ``.CH`` file.
    include_xic:
        For mass-spec ``.MS`` files, also decode scan-by-m/z extracted ion
        intensities. This can be much larger than the TIC.
    precision:
        Decimal places used to bin m/z values when ``include_xic`` is true.
    """

    input_path = Path(path).expanduser()
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    files = _chemstation_files(input_path)
    channels: list[ChannelData] = []

    for file_path in files:
        suffix = file_path.suffix.upper()
        try:
            if suffix == ".MS":
                channels.append(_read_ms(file_path, include_xic=include_xic, precision=precision))
            elif suffix == ".CH":
                channels.append(_read_ch(file_path))
        except UnsupportedChemStationFile:
            continue

    return channels


def read_chemstation_dataframe(
    path: str | Path,
    *,
    precision: int = 3,
    time_column: str = "time",
) -> "pd.DataFrame":
    """Read ChemStation MS channels into a pandas DataFrame.

    The first column is time in minutes. Remaining columns are the detected
    m/z values, using numeric pandas column labels such as ``44.0`` or
    ``45.0``.

    If a ``.D`` folder contains more than one readable ``.MS`` file, the first
    one with decoded m/z channels is used.
    """

    try:
        import pandas as pd
    except ImportError as exc:
        raise ImportError("read_chemstation_dataframe requires pandas.") from exc

    channels = read_chemstation(path, include_xic=True, precision=precision)
    for channel in channels:
        if channel.mz is None or channel.xic is None:
            continue

        dataframe = pd.DataFrame(channel.xic, columns=channel.mz)
        dataframe.insert(0, time_column, channel.time)
        return dataframe

    raise ValueError(f"No mass/charge channels found in ChemStation data: {path}")


def write_channel_csv(channel: ChannelData, output_path: str | Path) -> None:
    """Write one channel's time/signal data to CSV."""

    with Path(output_path).open("w", newline="") as handle:
        writer = csv.writer(handle)
        unit_label = f"signal ({channel.units})" if channel.units else "signal"
        writer.writerow(["time_min", unit_label])
        writer.writerows(zip(channel.time, channel.signal))


class UnsupportedChemStationFile(ValueError):
    """Raised when a ChemStation file version is not handled by this module."""


def _chemstation_files(path: Path) -> list[Path]:
    if path.is_dir():
        if path.suffix.upper() != ".D":
            raise ValueError(f"Expected a .D folder, got directory: {path}")
        return sorted(
            child
            for child in path.iterdir()
            if child.is_file()
            and not child.name.startswith(".")
            and child.suffix.upper() in {".MS", ".CH"}
        )

    if path.suffix.upper() in {".MS", ".CH"}:
        return [path]

    raise ValueError(f"Expected a .D folder, .MS file, or .CH file: {path}")


def _read_ms(path: Path, *, include_xic: bool, precision: int) -> ChannelData:
    with path.open("rb") as handle:
        version = _read_pascal_string(handle, 0, "ascii").strip()
        if version not in {"2", "20"}:
            raise UnsupportedChemStationFile(f"Unsupported .MS version {version!r}: {path}")

        metadata = _read_ms_metadata(handle)

        scans = _read_struct_at(handle, 278, ">I")[0]
        tic_offset = _read_struct_at(handle, 260, ">i")[0] * 2 - 2

        handle.seek(tic_offset + 4)
        time = [_read_struct(handle, ">i", skip=8)[0] / 60000.0 for _ in range(scans)]

        handle.seek(tic_offset + 8)
        tic = [float(_read_struct(handle, ">i", skip=8)[0]) for _ in range(scans)]

        channel = ChannelData(
            source=str(path),
            kind="MS TIC",
            time=time,
            signal=tic,
            metadata=metadata,
        )

        if include_xic:
            xic_offsets = _read_ms_xic_offsets(handle, tic_offset, scans)
            channel.mz, channel.xic = _read_ms_xic(handle, xic_offsets, scans, len(time), precision)

        return channel


def _read_ms_metadata(handle: BinaryIO) -> dict[str, object]:
    return {
        "version": _read_pascal_string(handle, 0, "ascii").strip(),
        "sample": {
            "name": _read_pascal_string(handle, 24, "ascii"),
            "description": _read_pascal_string(handle, 86, "ascii"),
            "sequence": _read_struct_at(handle, 252, ">h")[0],
            "vial": _read_struct_at(handle, 254, ">h")[0],
            "replicate": _read_struct_at(handle, 256, ">h")[0],
        },
        "method": {
            "name": _read_pascal_string(handle, 228, "ascii"),
            "operator": _read_pascal_string(handle, 148, "ascii"),
            "date_time": _read_pascal_string(handle, 178, "ascii"),
        },
        "instrument": {
            "name": _read_pascal_string(handle, 208, "ascii"),
            "inlet": _read_pascal_string(handle, 218, "ascii"),
        },
    }


def _read_ms_xic_offsets(handle: BinaryIO, tic_offset: int, scans: int) -> list[int]:
    handle.seek(tic_offset)
    return [_read_struct(handle, ">i", skip=8)[0] * 2 - 2 for _ in range(scans)]


def _read_ms_xic(
    handle: BinaryIO,
    offsets: Iterable[int],
    scans: int,
    time_count: int,
    precision: int,
) -> tuple[list[float], list[list[float]]]:
    rows: list[list[tuple[float, float]]] = []
    mz_values: list[float] = []

    for offset in offsets:
        handle.seek(offset)
        point_count = int((_read_struct(handle, ">h")[0] - 18) / 2 + 2)
        if point_count <= 0:
            rows.append([])
            continue

        handle.seek(offset + 18)
        mz_raw = [_read_struct(handle, ">H", skip=2)[0] for _ in range(point_count)]

        handle.seek(offset + 20)
        intensity_raw = [_read_struct(handle, ">H", skip=2)[0] for _ in range(point_count)]

        row: list[tuple[float, float]] = []
        for raw_mz, raw_intensity in zip(mz_raw, intensity_raw):
            mz = round((raw_mz / 20.0), precision)
            mantissa = raw_intensity & 0x3FFF
            exponent = raw_intensity >> 14
            intensity = float(mantissa * (8**exponent))
            row.append((mz, intensity))
            mz_values.append(mz)
        rows.append(row)

    mz_axis = sorted(set(mz_values))
    mz_index = {mz: index for index, mz in enumerate(mz_axis)}
    matrix = [[0.0 for _ in mz_axis] for _ in range(time_count)]

    for scan_index in range(min(scans, time_count)):
        for mz, intensity in rows[scan_index]:
            matrix[scan_index][mz_index[mz]] = intensity

    return mz_axis, matrix


def _read_ch(path: Path) -> ChannelData:
    with path.open("rb") as handle:
        version = _read_pascal_string(handle, 0, "ascii").strip()

        if version == "8":
            offsets = _legacy_ch_offsets()
            data_offset = (_read_struct_at(handle, 264, ">i")[0] - 1) * 512
            metadata = _read_ch_metadata(handle, version, offsets)
            signal = _delta_compression(handle, data_offset)
            xmin = _read_struct_at(handle, 282, ">i")[0] / 60000.0
            xmax = _read_struct_at(handle, 286, ">i")[0] / 60000.0
            header = _read_struct_at(handle, 542, ">i")[0]
            if header in {1, 2, 3}:
                signal = [value * 1.33321110047553 for value in signal]
            else:
                intercept = _read_struct_at(handle, 636, ">d")[0]
                slope = _read_struct_at(handle, 644, ">d")[0]
                signal = [value * slope + intercept for value in signal]

        elif version == "81":
            offsets = _legacy_ch_offsets()
            data_offset = (_read_struct_at(handle, 264, ">i")[0] - 1) * 512
            metadata = _read_ch_metadata(handle, version, offsets)
            signal = _double_delta_compression(handle, data_offset)
            xmin = _read_struct_at(handle, 282, ">f")[0] / 60000.0
            xmax = _read_struct_at(handle, 286, ">f")[0] / 60000.0
            intercept = _read_struct_at(handle, 636, ">d")[0]
            slope = _read_struct_at(handle, 644, ">d")[0]
            signal = [value * slope + intercept for value in signal]

        elif version == "179":
            offsets = _modern_ch_offsets()
            data_offset = (_read_struct_at(handle, 264, ">i")[0] - 1) * 512
            metadata = _read_ch_metadata(handle, version, offsets)
            signal = _double_array(handle, data_offset)
            xmin = _read_struct_at(handle, 282, ">f")[0] / 60000.0
            xmax = _read_struct_at(handle, 286, ">f")[0] / 60000.0
            intercept = _read_struct_at(handle, 4724, ">d")[0]
            slope = _read_struct_at(handle, 4732, ">d")[0]
            signal = [value * slope + intercept for value in signal]

        elif version == "181":
            offsets = _modern_ch_offsets()
            data_offset = (_read_struct_at(handle, 264, ">i")[0] - 1) * 512
            metadata = _read_ch_metadata(handle, version, offsets)
            signal = _double_delta_compression(handle, data_offset)
            xmin = _read_struct_at(handle, 282, ">f")[0] / 60000.0
            xmax = _read_struct_at(handle, 286, ">f")[0] / 60000.0
            intercept = _read_struct_at(handle, 4724, ">d")[0]
            slope = _read_struct_at(handle, 4732, ">d")[0]
            signal = [value * slope + intercept for value in signal]

        else:
            raise UnsupportedChemStationFile(f"Unsupported .CH version {version!r}: {path}")

    return ChannelData(
        source=str(path),
        kind="Detector Channel",
        time=_linspace(xmin, xmax, len(signal)),
        signal=[float(value) for value in signal],
        units=str(metadata.get("instrument", {}).get("units", "")),
        metadata=metadata,
    )


def _legacy_ch_offsets() -> dict[str, int]:
    return {
        "sample": 24,
        "description": 86,
        "method": 228,
        "operator": 148,
        "date": 178,
        "instrument": 218,
        "inlet": 208,
        "units": 580,
    }


def _modern_ch_offsets() -> dict[str, int]:
    return {
        "sample": 858,
        "description": 1369,
        "method": 2574,
        "operator": 1880,
        "date": 2391,
        "instrument": 2533,
        "inlet": 2492,
        "units": 4172,
    }


def _read_ch_metadata(handle: BinaryIO, version: str, offsets: dict[str, int]) -> dict[str, object]:
    encoding = "ascii" if version in {"8", "81"} else "utf-16le"
    return {
        "version": version,
        "sample": {
            "name": _read_pascal_string(handle, offsets["sample"], encoding),
            "description": _read_pascal_string(handle, offsets["description"], encoding),
            "sequence": _read_struct_at(handle, 252, ">h")[0],
            "vial": _read_struct_at(handle, 254, ">h")[0],
            "replicate": _read_struct_at(handle, 256, ">h")[0],
        },
        "method": {
            "name": _read_pascal_string(handle, offsets["method"], encoding),
            "operator": _read_pascal_string(handle, offsets["operator"], encoding),
            "date_time": _read_pascal_string(handle, offsets["date"], encoding),
        },
        "instrument": {
            "name": _read_pascal_string(handle, offsets["instrument"], encoding),
            "inlet": _read_pascal_string(handle, offsets["inlet"], encoding),
            "units": _read_pascal_string(handle, offsets["units"], encoding),
        },
    }


def _delta_compression(handle: BinaryIO, offset: int) -> list[float]:
    stop = _file_size(handle)
    handle.seek(offset)
    signal: list[float] = []
    previous = 0

    while handle.tell() < stop:
        word = _read_struct(handle, ">h")[0]
        if (word << 12) == 0:
            break

        for _ in range(word & 0x0FFF):
            delta = _read_struct(handle, ">h")[0]
            if delta != -32768:
                previous += delta
            else:
                previous = _read_struct(handle, ">i")[0]
            signal.append(float(previous))

    return signal


def _double_delta_compression(handle: BinaryIO, offset: int) -> list[float]:
    stop = _file_size(handle)
    handle.seek(offset)
    signal: list[float] = []
    value = 0
    delta_sum = 0

    while handle.tell() < stop:
        delta = _read_struct(handle, ">h")[0]
        if delta != 32767:
            delta_sum += delta
            value += delta_sum
        else:
            high = _read_struct(handle, ">h")[0]
            low = _read_struct(handle, ">I")[0]
            value = high * 4294967296 + low
            delta_sum = 0
        signal.append(float(value))

    return signal


def _double_array(handle: BinaryIO, offset: int) -> list[float]:
    stop = _file_size(handle)
    count = (stop - offset) // 8
    handle.seek(offset)
    return [float(_read_struct(handle, "<d")[0]) for _ in range(count)]


def _read_pascal_string(handle: BinaryIO, offset: int, encoding: str) -> str:
    handle.seek(offset)
    length_raw = handle.read(1)
    if not length_raw:
        return ""
    length = length_raw[0]
    byte_count = length * 2 if encoding.lower().replace("_", "-") in {"utf-16le", "utf-16-le"} else length
    raw = handle.read(byte_count)
    return raw.decode(encoding, errors="ignore").strip().strip("\x00").strip()


def _read_struct_at(handle: BinaryIO, offset: int, fmt: str) -> tuple[object, ...]:
    handle.seek(offset)
    return _read_struct(handle, fmt)


def _read_struct(handle: BinaryIO, fmt: str, *, skip: int = 0) -> tuple[object, ...]:
    size = struct.calcsize(fmt)
    raw = handle.read(size)
    if len(raw) != size:
        raise EOFError("Unexpected end of ChemStation file")
    values = struct.unpack(fmt, raw)
    if skip:
        handle.seek(skip, 1)
    return values


def _file_size(handle: BinaryIO) -> int:
    current = handle.tell()
    handle.seek(0, 2)
    size = handle.tell()
    handle.seek(current)
    return size


def _linspace(start: float, stop: float, count: int) -> list[float]:
    if count <= 0:
        return []
    if count == 1:
        return [start]
    step = (stop - start) / (count - 1)
    return [start + step * index for index in range(count)]


def _channel_summary(channel: ChannelData) -> dict[str, object]:
    signal_min = min(channel.signal) if channel.signal else math.nan
    signal_max = max(channel.signal) if channel.signal else math.nan
    return {
        "source": channel.source,
        "kind": channel.kind,
        "points": len(channel.time),
        "time_min": channel.time[0] if channel.time else None,
        "time_max": channel.time[-1] if channel.time else None,
        "signal_min": signal_min,
        "signal_max": signal_max,
        "units": channel.units,
        "has_xic": channel.xic is not None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Read raw Agilent/ChemStation .D, .MS, and .CH data.")
    parser.add_argument("path", help="ChemStation .D folder or .MS/.CH file")
    parser.add_argument("--precision", type=int, default=3, help="Decimal places used to bin m/z channels")
    parser.add_argument("--include-xic", action="store_true", help="Decode .MS extracted ion matrix")
    parser.add_argument("--xic-csv", help="Write m/z channel DataFrame to this CSV file")
    parser.add_argument("--csv", help="Write the first channel's time/signal data to this CSV file")
    args = parser.parse_args()

    channels = read_chemstation(args.path, include_xic=args.include_xic, precision=args.precision)
    print(json.dumps([_channel_summary(channel) for channel in channels], indent=2))

    if args.xic_csv:
        dataframe = read_chemstation_dataframe(args.path, precision=args.precision)
        dataframe.to_csv(args.xic_csv, index=False)

    if args.csv:
        if not channels:
            raise SystemExit("No readable ChemStation channels found.")
        write_channel_csv(channels[0], args.csv)


if __name__ == "__main__":
    main()
# %%
path = '/Users/gerard/Downloads/Test_samples-selected/H2504022.D'
df = read_chemstation_dataframe(path)

import matplotlib.pyplot as plt
fig = plt.figure()
plt.plot(df['time'], df[1022.1])
plt.show()


# %%
import glob, os
path = '/Users/gerard/Downloads/Test_samples-selected/'
os.chdir(path)
for file in glob.glob("*.D"):
    
