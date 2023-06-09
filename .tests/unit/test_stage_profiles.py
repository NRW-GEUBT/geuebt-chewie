import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load


sys.path.insert(0, os.path.dirname(__file__))


def test_stage_profiles():
    with TemporaryDirectory() as tmpdir:
        # Modify paths to link to your test data and script
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/stage_profiles/data")
        expected_path = PurePosixPath(".tests/unit/stage_profiles/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/stage_profiles.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)
        
        # run function
        sys.path.insert(0, workdir)
        from stage_profiles import main # import main from your script
        main(
            jsons=[
                os.path.join(workdir, '16-LI00296-0.chewiesnake.json'),
                os.path.join(workdir, '16-LI00501-0.chewiesnake.json')
            ],
            qc=os.path.join(workdir, 'qc_status.json'),
            outdir=os.path.join(workdir),
        )

        # Insert your tests here
        # Ensure number and naming of files
        assert os.path.isfile(os.path.join(workdir, '16-LI00296-0.json'))
        assert not os.path.isfile(os.path.join(workdir, '16-LI00501-0.json'))
        # ensure dict correctness
        with open(
            os.path.join(workdir, '16-LI00296-0.json'), 'r'
        ) as res, open(
            os.path.join(expected_path, 'isolate_sheet.json'), 'r'
        ) as expect:
            assert load(res) == load(expect)
