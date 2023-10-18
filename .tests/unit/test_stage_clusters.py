import os
import sys

from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath
from json import load, dumps


sys.path.insert(0, os.path.dirname(__file__))


def test_stage_clusters():
    with TemporaryDirectory() as tmpdir:
        # Modify paths to link to your test data and script
        workdir = os.path.join(Path(tmpdir), "workdir")
        data_path = PurePosixPath(".tests/unit/stage_clusters/data")
        expected_path = PurePosixPath(".tests/unit/stage_clusters/expected")
        script_path = PurePosixPath(".tests/../workflow/scripts/stage_clusters.py")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copy(script_path, workdir)

        # run function
        sys.path.insert(0, workdir)
        from stage_clusters import main
        main(
            cluster_info=os.path.join(workdir, 'main_cluster_info.tsv'),
            orphans=os.path.join(workdir, 'orphans.tsv'),
            subcluster_info=os.path.join(workdir, 'sub_cluster_info.tsv'),
            distances=os.path.join(workdir, "distance_matrix.tsv"),
            prefix="LIS",
            main_dist=10,
            sub_dist=5,
            organism="Listeria monocytogenes",
            dirout=os.path.join(workdir),
            mergedout=os.path.join(workdir, 'result.json')
        )

        # Insert your tests here
        # Ensure number and naming of files
        expected_files = [
            "LIS-2.json",
            "LIS-3.json",
            "LIS-4.json",
            "LIS-5.json",
            "LIS-6.json",
            "LIS-7.json",
            "LIS-8.json",
            "LIS-orphans.json"
        ]
        for expected_file in expected_files:
            assert os.path.isfile(os.path.join(workdir, expected_file))
        for unexpected_file in ["LIS-1.json", "LIS-9.json"]:
            assert not os.path.isfile(os.path.join(workdir, unexpected_file))
        # ensure dict correctness
        for expected_file in expected_files:
            with open(
                os.path.join(workdir, expected_file), 'r'
            ) as res, open(
                os.path.join(expected_path, expected_file), 'r'
            ) as expect:
                assert load(res) == load(expect)
        with open(
            os.path.join(workdir, 'result.json'), 'r'
        ) as res, open(
            os.path.join(expected_path, 'merged.json'), 'r'
        ) as expect:
            # make it sets because order is not conserved
            res_set = set()
            for ele in load(res):
                res_set.add(dumps(ele, sort_keys=True))
            exp_set = set()
            for ele in load(expect):
                exp_set.add(dumps(ele, sort_keys=True))
            assert res_set == exp_set
