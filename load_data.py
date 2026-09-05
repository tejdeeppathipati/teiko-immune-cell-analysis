"""Validate cell-count.csv and load it into a normalized SQLite database."""

from __future__ import annotations

import os
import sqlite3
from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
CSV_PATH = BASE_DIR / "cell-count.csv"
DATABASE_PATH = BASE_DIR / "cell_counts.db"
TEMP_DATABASE_PATH = BASE_DIR / "cell_counts.db.tmp"

METADATA_COLUMNS = [
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample_type",
]
CELL_POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
REQUIRED_COLUMNS = (
    METADATA_COLUMNS + ["sample", "time_from_treatment_start"] + CELL_POPULATIONS
)


SCHEMA_SQL = """
CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL,
    condition TEXT NOT NULL,
    age INTEGER NOT NULL CHECK (age >= 0),
    sex TEXT NOT NULL,
    treatment TEXT NOT NULL,
    response TEXT NULL,
    FOREIGN KEY (project_id) REFERENCES projects(project_id)
);

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    subject_id TEXT NOT NULL,
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL,
    FOREIGN KEY (subject_id) REFERENCES subjects(subject_id)
);

CREATE TABLE cell_counts (
    sample_id TEXT NOT NULL,
    population TEXT NOT NULL,
    cell_count INTEGER NOT NULL CHECK (cell_count >= 0),
    PRIMARY KEY (sample_id, population),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

CREATE INDEX idx_subject_analysis
ON subjects(condition, treatment, response);

CREATE INDEX idx_sample_analysis
ON samples(sample_type, time_from_treatment_start);
"""


def _validate_integral_nonnegative(data: pd.DataFrame, columns: list[str]) -> None:
    """Validate numeric count-like columns without modifying their source values."""
    for column in columns:
        numeric = pd.to_numeric(data[column], errors="coerce")
        if numeric.isna().any():
            raise ValueError(f"Column {column!r} contains missing or nonnumeric values")
        if (numeric < 0).any():
            raise ValueError(f"Column {column!r} contains negative values")
        if (numeric % 1 != 0).any():
            raise ValueError(f"Column {column!r} contains non-integral values")


def read_and_validate_csv(csv_path: Path = CSV_PATH) -> pd.DataFrame:
    """Read the source CSV and enforce assumptions required by the schema."""
    if not csv_path.exists():
        raise FileNotFoundError(f"Input file not found: {csv_path}")

    data = pd.read_csv(csv_path)
    missing_columns = sorted(set(REQUIRED_COLUMNS) - set(data.columns))
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

    if data["sample"].isna().any() or not data["sample"].is_unique:
        raise ValueError("Sample IDs must be present and unique")

    required_text = [
        "project",
        "subject",
        "condition",
        "sex",
        "treatment",
        "sample_type",
    ]
    if data[required_text].isna().any().any():
        raise ValueError("Required subject and sample metadata cannot be missing")

    _validate_integral_nonnegative(
        data, ["age", "time_from_treatment_start", *CELL_POPULATIONS]
    )

    inconsistent = []
    for column in METADATA_COLUMNS:
        counts = data.groupby("subject", dropna=False)[column].nunique(dropna=False)
        if (counts > 1).any():
            inconsistent.append(column)
    if inconsistent:
        raise ValueError("Subject-level metadata is inconsistent for: " + ", ".join(inconsistent))

    return data


def create_and_load_database(data: pd.DataFrame) -> dict[str, int]:
    """Build a fresh database and atomically replace the prior database on success."""
    # A temporary database keeps the previous valid run intact if loading fails.
    TEMP_DATABASE_PATH.unlink(missing_ok=True)
    connection = sqlite3.connect(TEMP_DATABASE_PATH)
    try:
        connection.execute("PRAGMA foreign_keys = ON")
        connection.executescript(SCHEMA_SQL)

        projects = [(value,) for value in sorted(data["project"].unique())]
        subject_rows = data.drop_duplicates("subject")
        subjects = [
            (
                row.subject,
                row.project,
                row.condition,
                int(row.age),
                row.sex,
                row.treatment,
                None if pd.isna(row.response) else row.response,
            )
            for row in subject_rows.itertuples(index=False)
        ]
        samples = [
            (
                row.sample,
                row.subject,
                row.sample_type,
                int(row.time_from_treatment_start),
            )
            for row in data.itertuples(index=False)
        ]

        # Store measurements in long form so new populations do not require schema changes.
        long_counts = data[["sample", *CELL_POPULATIONS]].melt(
            id_vars="sample", var_name="population", value_name="cell_count"
        )
        measurements = [
            (row.sample, row.population, int(row.cell_count))
            for row in long_counts.itertuples(index=False)
        ]

        with connection:
            connection.executemany("INSERT INTO projects VALUES (?)", projects)
            connection.executemany(
                "INSERT INTO subjects VALUES (?, ?, ?, ?, ?, ?, ?)", subjects
            )
            connection.executemany(
                "INSERT INTO samples VALUES (?, ?, ?, ?)", samples
            )
            connection.executemany(
                "INSERT INTO cell_counts VALUES (?, ?, ?)", measurements
            )

        table_counts = {
            table: connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            for table in ("projects", "subjects", "samples", "cell_counts")
        }
        foreign_key_errors = connection.execute("PRAGMA foreign_key_check").fetchall()
        if foreign_key_errors:
            raise RuntimeError(f"Foreign-key validation failed: {foreign_key_errors}")
    except Exception:
        connection.close()
        TEMP_DATABASE_PATH.unlink(missing_ok=True)
        raise
    else:
        connection.close()

    os.replace(TEMP_DATABASE_PATH, DATABASE_PATH)
    return table_counts


def main() -> None:
    """Validate the source file and rebuild the project database."""
    data = read_and_validate_csv()
    counts = create_and_load_database(data)
    summary = ", ".join(f"{table}={count:,}" for table, count in counts.items())
    print(f"Loaded {CSV_PATH.name} into {DATABASE_PATH.name}: {summary}")


if __name__ == "__main__":
    main()
