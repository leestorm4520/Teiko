from pathlib import Path
import pandas as pd
from sqlalchemy import (
    CheckConstraint, ForeignKey, Index, Integer, String, create_engine,
    func, insert, select)
from sqlalchemy.orm import (DeclarativeBase, Mapped, Session, mapped_column, relationship)

ROOT = Path(__file__).resolve().parent
CSV = ROOT/"cell-count.csv"
DB = ROOT/"cell_counts.db"
POPULATIONS=["b_cell", "cd8_t_cell", "cd4_t_cell","nk_cell", "monocyte"]

class Base(DeclarativeBase):
    pass

class Project(Base):
    __tablename__ = "projects"

    project_id: Mapped[str] = mapped_column(String, primary_key=True)

    subjects: Mapped[list["Subject"]] = relationship(back_populates="project")


class Subject(Base):
    """One row per patient; attributes are constant across a subject's samples."""
    __tablename__ = "subjects"
    __table_args__ = (
        CheckConstraint("sex IN ('M','F')", name="ck_subject_sex"),
        CheckConstraint("response IN ('yes','no') OR response IS NULL",
                        name="ck_subject_response"),
        Index("idx_subjects_filter", "condition", "treatment", "response", "sex"),
    )

    subject_id: Mapped[str] = mapped_column(String, primary_key=True)
    project_id: Mapped[str] = mapped_column(ForeignKey("projects.project_id"))
    condition: Mapped[str] = mapped_column(String)
    age: Mapped[int | None] = mapped_column(Integer)
    sex: Mapped[str | None] = mapped_column(String)
    treatment: Mapped[str] = mapped_column(String)
    response: Mapped[str | None] = mapped_column(String)

    project: Mapped[Project] = relationship(back_populates="subjects")
    samples: Mapped[list["Sample"]] = relationship(back_populates="subject")


class Sample(Base):
    """One row per biological sample."""
    __tablename__ = "samples"
    __table_args__ = (
        Index("idx_samples_subject", "subject_id"),
        Index("idx_samples_filter", "sample_type", "time_from_treatment_start"),
    )

    sample_id: Mapped[str] = mapped_column(String, primary_key=True)
    subject_id: Mapped[str] = mapped_column(ForeignKey("subjects.subject_id"))
    sample_type: Mapped[str] = mapped_column(String)
    time_from_treatment_start: Mapped[int] = mapped_column(Integer)

    subject: Mapped[Subject] = relationship(back_populates="samples")
    counts: Mapped[list["CellCount"]] = relationship(back_populates="sample")


class CellCount(Base):
    """Long format: one row per sample x population, so a new population
    needs no schema change."""
    __tablename__ = "cell_counts"
    __table_args__ = (
        CheckConstraint("count >= 0", name="ck_count_nonneg"),
        Index("idx_counts_pop", "population"),
    )

    sample_id: Mapped[str] = mapped_column(ForeignKey("samples.sample_id"),
                                           primary_key=True)
    population: Mapped[str] = mapped_column(String, primary_key=True)
    count: Mapped[int] = mapped_column(Integer)

    sample: Mapped[Sample] = relationship(back_populates="counts")

def records(frame: pd.DataFrame) -> list[dict]:
    """DataFrame -> list of dicts with NaN converted to None (SQL NULL)."""
    return frame.astype(object).where(frame.notna(), None).to_dict("records")

def main() -> None:
    df = pd.read_csv(CSV)
    if DB.exists():
        DB.unlink()
    engine = create_engine(f"sqlite:///{DB}")
    Base.metadata.create_all(engine)

    subjects = df[["subject", "project", "condition", "age", "sex",
                   "treatment", "response"]].drop_duplicates("subject")
    subjects = subjects.rename(columns={"subject": "subject_id", "project": "project_id"})
    samples = df[["sample", "subject", "sample_type", "time_from_treatment_start"]]
    samples = samples.rename(columns={"sample": "sample_id", "subject": "subject_id"})
    counts = df.melt(id_vars="sample", value_vars=POPULATIONS,
                     var_name="population", value_name="count")
    counts = counts.rename(columns={"sample": "sample_id"})

    with Session(engine) as session, session.begin():
        session.execute(insert(Project),
                        [{"project_id": p} for p in sorted(df["project"].unique())])
        session.execute(insert(Subject), records(subjects))
        session.execute(insert(Sample), records(samples))
        session.execute(insert(CellCount), records(counts))

    with Session(engine) as session:
        for model in (Project, Subject, Sample, CellCount):
            n = session.scalar(select(func.count()).select_from(model))
            print(f"{model.__tablename__}: {n} rows")
    print(f"Database written to {DB}")


if __name__ == "__main__":
    main()