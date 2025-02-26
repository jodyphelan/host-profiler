import argparse
import pathogenprofiler as pp
from .reformat import get_result
import host_profiler as hp
import shutil

def profile(args:argparse.Namespace):
    args.conf = pp.get_db(db_dir=args.db_dir, db_name = args.db)
    
    pp.process_args(args)

    variants_profile = pp.run_profiler(args)
    
    if args.bam:
        qc = pp.run_bam_qc(args,only_targets=True)
    elif args.vcf:
        qc = pp.run_vcf_qc(args)
    else:
        qc = None

    result = get_result(
        args=args,
        variants=variants_profile,
        qc=qc,
        software_version=hp.__version__,
        db_version=args.conf['version']
    )

    
    open(f"{args.dir}/{args.prefix}.results.json","w").write(result.model_dump_json(indent=4))

    for extension in ['.bam','.targets_for_profile.csq.vcf.gz']:
        shutil.move(f"{args.files_prefix}{extension}",f"{args.dir}/{args.prefix}{extension}")

    pp.run_cmd(f"rm {args.files_prefix}*")
