import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_qc_status():
    with TemporaryDirectory() as tmpdir:
        # Modify paths to link to your test data and script
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/qc_status/data")
        expected_path = PurePosixPath(".tests/unit/qc_status/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/qc_status.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # run function
        sys.path.insert(0, workdir)
        from qc_status import main
        # import main from your script
        main(
            jsons=[
                os.path.join(workdir, '16-LI00296-0.chewiesnake.json'),
                os.path.join(workdir, '16-LI00501-0.chewiesnake.json'),
                os.path.join(workdir, '16-LI00919-0.chewiesnake.json')
            ],
            max_missing_loci=0.05,
            json_path=os.path.join(workdir, 'result.json'),
            list_path=os.path.join(workdir, 'listout.txt'),
        )

        # Insert your tests here
        # check dicts
        with open(
            os.path.join(workdir, 'result.json'), 'r'
        ) as res, open(
            os.path.join(expected_path, 'qc_status.json'), 'r'
        ) as expect:
            assert load(res) == load(expect)
        # check list
        with open(
            os.path.join(workdir, 'listout.txt'), 'r'
        ) as res, open(
            os.path.join(expected_path, 'sample_list.txt'), 'r'
        ) as expect:
            assert res.readlines() == expect.readlines()
