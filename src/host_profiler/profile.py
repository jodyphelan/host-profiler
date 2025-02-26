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

    files_to_keep = {
        f'{args.files_prefix}.bam':f"{args.dir}/{args.prefix}.bam",
        f'{args.files_prefix}.targets_for_profile.csq.vcf.gz':f"{args.dir}/{args.prefix}.vcf.gz"
    }
    for src,dest in files_to_keep.items():
        shutil.move(src,dest)

    pp.run_cmd(f"rm {args.files_prefix}*")
