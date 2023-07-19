import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath


sys.path.insert(0, os.path.dirname(__file__))


def test_failsafe_chewie_join_input():
    with TemporaryDirectory() as tmpdir:
        # Modify paths to link to your test data and script
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/failsafe_chewie_join_input/data")
        expected_path = PurePosixPath(".tests/unit/failsafe_chewie_join_input/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/failsafe_chewie_join_input.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # run function
        sys.path.insert(0, workdir)
        from failsafe_chewie_join_input import main
        main(
            profiles_template=os.path.join(workdir, 'profiles_filled.tsv'),
            ext_main=os.path.join(workdir, 'main.tsv'),
            ext_sub=os.path.join(workdir, 'sub.tsv'),
            ext_profiles=os.path.join(workdir, 'profiles.tsv'),
            ext_timestamps=os.path.join(workdir, 'timestamps.tsv'),
            ext_statistics=os.path.join(workdir, 'statistics.tsv'),
            ext_main_out=os.path.join(workdir, 'main_out.tsv'),
            ext_sub_out=os.path.join(workdir, 'sub_out.tsv'),
            ext_profiles_out=os.path.join(workdir, 'profiles_out.tsv'),
            ext_timestamps_out=os.path.join(workdir, 'timestamps_out.tsv'),
            ext_statistics_out=os.path.join(workdir, 'statistics_out.tsv')
        )

        # Insert your tests here
        filenames = [
            "main_out.tsv",
            "sub_out.tsv",
            "profiles_out.tsv",
            "timestamps_out.tsv",
            "statistics_out.tsv"
        ]
        for filename in filenames:
            with open(
                os.path.join(workdir, filename), 'r'
            ) as res, open(
                os.path.join(expected_path, filename), 'r'
            ) as expect:
                assert res.readlines() == expect.readlines()
