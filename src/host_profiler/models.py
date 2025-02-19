from pydantic import BaseModel
from typing import List, Optional, Union
from pathogenprofiler.models import Variant, BamQC, VcfQC

__model_schema_version__ = "0.1.0"

class Pipeline(BaseModel):
    """
    A class to hold information about the NTM-Profiler pipeline
    
    Attributes
    ----------
    software_versio  : str
        Host-Profiler version
    db_version : dict
        Host-Profiler database version
    software : List[dict]
        Software used in the pipeline
    """
    software_version: str
    db_version: Optional[dict]
    software: List[dict]

class Result(BaseModel):
    schema_version: str = __model_schema_version__
    pipeline: Pipeline
    id: str

class ProfileResult(Result):
    result_type: str = 'Profile'    
    notes: List[str] = []
    variants: List[Variant] = []
    fail_variants: List[Variant] = []
    qc: Union[BamQC, VcfQC]

