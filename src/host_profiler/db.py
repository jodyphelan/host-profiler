import tarfile
import zipfile
import os
import pathogenprofiler as pp
import sys
import logging

def install_db(args):
    """
    Install a host-profiler database from a compressed file.

    Parameters
    ----------
    compressed_file : str
        Path to the compressed file to install
    """
    # check that the database folders are set up
    if not os.path.isdir(args.db_dir):
        os.makedirs(args.db_dir)
    pp.create_snpeff_directories(args.db_dir)

    # check the file exists
    if not os.path.isfile(args.archive):
        raise FileNotFoundError(f"File not found: {args.archive}")

    if args.archive.endswith('.tar.gz') or args.archive.endswith('.tgz'):
        # uncompress into args.db_dir
        tar = tarfile.open(args.archive, "r:gz")
        logging.info(f"Extracting {args.archive} to {args.db_dir}")
        tar.extractall(path=args.db_dir)
    elif args.archive.endswith('.zip'):
        # uncompress into args.db_dir
        with zipfile.ZipFile(args.archive, 'r') as zip_ref:
            logging.info(f"Extracting {args.archive} to {args.db_dir}")
            zip_ref.extractall(args.db_dir)
    else:
        raise ValueError(f"Unknown archive format: {args.archive}")
    
def list_db(args):
    """
    List the available host-profiler databases.
    """
    dbs = pp.list_db(args.db_dir)
    for db in dbs:
        if 'version' in db:
            d = dict(**db['version'], location=f"{args.db_dir}/{db['version']['name']}")
            print(d)
            sys.stdout.write("%(name)s\t%(commit)s\t%(author)s\t%(date)s\t%(location)s\n" % d)