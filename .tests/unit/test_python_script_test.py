import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_python_script_test():
    with TemporaryDirectory() as tmpdir:
        # Modify paths to link to your test data and script
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/python_script_test/data")
        expected_path = PurePosixPath(".tests/unit/python_script_test/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/python_script_test.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(schema_path, workdir)
        shutil.copy(script_path, workdir)
        
        # run function
        sys.path.insert(0, workdir)
        from python_script_test import main # import main from your script
        main(
            arg1=os.path.join(workdir, 'data.json'),
            arg2=os.path.join(workdir, 'moredata.tsv'),
            arg3=os.path.join(workdir, 'result.json'),
        )

        # Insert your tests here
        with open(
            os.path.join(workdir, 'result.json'), 'r'
        ) as res, open(
            os.path.join(expected_path, 'expected_res.json'), 'r'
        ) as expect:
            assert res == expect
