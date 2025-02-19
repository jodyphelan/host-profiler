from .models import Pipeline, Result, ProfileResult
from pathogenprofiler.models import Variant, BamQC, VcfQC
from typing import List, Union
import argparse

def get_result(
        args: argparse.Namespace,
        variants: List[Variant],
        qc: Union[BamQC, VcfQC],
        software_version: str,
        db_version: dict = None
    ):
    return ProfileResult(
        pipeline=Pipeline(
            software_version=software_version,
            db_version=db_version,
            software=[]
        ),
        id=args.prefix,
        notes=[],
        variants=variants,
        fail_variants=[],
        qc=qc
    )

